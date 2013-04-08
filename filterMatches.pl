#!/bin/perl -w
use strict;
use diagnostics;
use Parallel::ForkManager;
$SIG{'INT'} = 'INT_handler';# manage the killing signal => remove temp files :D
my $max_process = 16;

#AUTHOR BENOIT CARRERES
#multi threads is set to 8
# for final table :"Probe.Set.Name\tSubject id\ts. start\ts. end\tmissmatch\tgaps\tPM Candidates(single multiple)".
# if no response\t the 4 last columns will be empty

my $argSize = $#ARGV;
print "there is $argSize arguments\n\n";

my $originPath = `pwd`;
chomp $originPath;
print "the origin folder is : $originPath\n";
my @folders = split("\n",`ls -la | egrep '^d' | awk '{print \$9}' | egrep -v '^\.\.?\$'`);


#separateCandidates();
#getherSequencesHomo();

blastHomoWithReads();
#estimatExpressionLevel();
#getherSequencesHetero();


###############################################################################
############Separate matched transcripts by homogeneous/non size###############
###############################################################################
sub separateCandidates {
print "==================Start filtering================================\n";
my $pmComp = new Parallel::ForkManager($max_process);
foreach my $workingFolder (@folders) {	
	#start multithreading
	my $pidComp = $pmComp->start and next;
	#now several jobs
	print "$!=>$workingFolder\n";
	my $analyseFilename = $workingFolder."vsAll";
	my $heteroSizeFilename = $workingFolder."vsAll_hom";
	my $homoSizeFilename = $workingFolder."vsAll_het";
	open (MATCHES, "<$analyseFilename") || die ("$!=>cannot open file to read results of $workingFolder\n");
	open (HETERO, ">$heteroSizeFilename") || die ("$!=>cannot open file to write results of $workingFolder\n");
	open (HOMO, ">$homoSizeFilename") || die ("$!=>cannot open file to write results of $workingFolder\n");
	while (<MATCHES>){
		my $line = $_;
		my @cells = split("\t",$line);
		my @sizes;
		foreach my $cell (@cells) {
			if (my ($start, $stop) = $cell =~ /[:](\d+)[-](\d+)$/) {
				push (@sizes,abs($stop-$start))
			}
		}
		#print "number of columns = $#cells\n";
		my $isHeterogenous = 0;
		for (my $i=1; $i<=$#sizes; $i++) {
			#print "size = $sizes[$i]\n";
			if ($sizes[0] != $sizes[$i]) {
				$isHeterogenous = 1;
			}# else {
			#	print HOMO ("$line");
			#}
		}
		if ($isHeterogenous == 1) {
			print HETERO ("$line");
		} else {
			print HOMO ("$line");
		}
	}
	close (MATCHES);
	close (HETERO);
	close (HOMO);
	#end of the multithreads
	$pmComp->finish;
}
#checks if there is nothing else running
$pmComp->wait_all_children;
}

###############################################################################
#########################gether homo sequences from db#########################
###############################################################################
sub getherSequencesHomo {
print "==================Start gethering homo seqs==========================\n";
my $pmGeth = new Parallel::ForkManager($max_process);
foreach (@folders) {
		`cp -r $_ /dev/shm/$_`;
	}
foreach my $workingFolder (@folders) {	
	#start multithreading
	my $pidGeth = $pmGeth->start and next;
	#now several jobs
	print "$pidGeth=>$workingFolder\n";

	my $homoFilename = $workingFolder."vsAll_hom";
	open (VSALLHOM, $homoFilename) || die ("$pidGeth=>pb accessing $homoFilename file from the folder: $originPath\n");
	my $heteroFilename = $workingFolder."vsAll_het";
	open (VSALLHET, $heteroFilename) || die ("$pidGeth=>pb accessing $heteroFilename file from the folder: $originPath\n");
	
	chdir($workingFolder) || die ("$pidGeth=>cannot got to directory : $_");
	my $currentPath = `pwd`;
	my $homoSeqs = "mRNA_hom.fna";
	open (HOM, ">$homoSeqs") || die ("$pidGeth=>cannot open file to write selected homogenous mRNA sequences of $workingFolder\n");
	my $heteroSeqs = "mRNA_het.fna";
	open (HET, ">$heteroSeqs") || die ("$pidGeth=>cannot open file to write selected heterogenous mRNA sequences of $workingFolder\n");
	print "$pidGeth=>the current folder is : $currentPath\n";
	
	while (my $line = <VSALLHOM>){
		chomp($line);
		#print "full line = $_\n";
		if ($line =~ /^$workingFolder[|]([-\d\w_]*?[:]\d+[-]\d+)\s.*/){
			#print "getheredID = $1\n";
			print HOM `fastacmd -d /dev/shm/$workingFolder/mRNA.fna -s $1`;
			if (`echo $?` != 0) {
				print "missing/error on fasta CMD, trying ID without start/stop\n";
				$line =~ /^$workingFolder[|]([-\d\w_]*?)[:]\d+[-]\d+\s.*/;
				print HOM `fastacmd -d /dev/shm/$workingFolder/mRNA.fna -s $1`;
			}		
		} else {
			warn "line not matching... pb in the file format:\n $line\n";
		}
	}
	close (VSALLHOM);
	close (HOM);
	while (my $line = <VSALLHET>){
		chomp($line);
		if ($line =~ /^$workingFolder[|]([-\d\w_]*?[:]\d+[-]\d+)\s.*/){
			#print "getheredID = $1\n";
			print HET `fastacmd -d /dev/shm/$workingFolder/mRNA.fna -s $1`;
		} else {
			warn "line not matching... pb in the file format\n$line\n";
		}
	}
	close (VSALLHET);
	close (HET);
	chdir("..") || die ("$pidGeth=>cannot come back to origin (previous) path\n");
	#end of the multithreads
	$pmGeth->finish;
}
#checks if there is nothing else running
$pmGeth->wait_all_children;

formatDB("hom");


}

###############################################################################
#########Build mRNA database out of the annotation (contig<-->genes)###########
###############################################################################
sub formatDB {
print "==================Building mRNA DB================================\n";
my $type = shift;
my $pmDB = new Parallel::ForkManager($max_process);
foreach (@folders) {
	#start multithreading
	my $pidDB = $pmDB->start and next;
	print "$pidDB=>$_\n";
	chdir($_) || die ("$pidDB=>cannot got to directory : $_\n");
	my $currentPath = `pwd`;
	print "$pidDB=>the current folder is : $currentPath\n";
	my $fileName = "mRNA_".$type.".fna";

	`formatdb -i $fileName -p F -o`;

	chdir("..") || die ("$pidDB=>cannot come back to origin (previous) path\n");
	#end of the multithreads
	$pmDB->finish;
}
#checks if there is nothing else running
$pmDB->wait_all_children;
}
###############################################################################


###############################################################################
########################Align/count reads on candidates########################
###############################################################################
sub blastHomoWithReads {
print "==================Expression Estimation================================\n";
my $pmEe = new Parallel::ForkManager($max_process);
foreach (@folders) {
		`mkdir /dev/shm/$_`;
		`cp $_/mRNA_hom* /dev/shm/$_`;
	}
foreach my $workingFolder (@folders) {
	#start multithreading
	my $pidEe = $pmEe->start and next;
	#now several jobs
	print "$pidEe=>$workingFolder\n";
	my $readsFilename = $workingFolder.".fna";
	open (READS, "<$readsFilename") || die ("$pidEe=>cannot open file to read reads: $readsFilename\n");
	
	chdir($workingFolder) || die ("$pidEe=>cannot got to directory : $_");
	my $currentPath = `pwd`;
	print "$pidEe=>the current folder is : $currentPath\n";
	open (RESULT, ">homMatchingReads") || die ("$pidEe=>pb writing/creating result file: $currentPath\n");
	my %resultsHash;
#	my $count = 0;
	while (<READS>){
		chomp;
		my $idQuery = $_;
		my $sequence;
		
		if ($idQuery =~ /^>.*/) {
			$idQuery =~ s/^(>\w{3}[|][-\d\w_]*?)\s.*length_(\d+)(\s.*)/$1:1-$2/;##if the id line does not have the range, add one from 1 to length
			$idQuery =~ s/^>\w{3}[|]([-\d\w_]*?[:]\d+[-]\d+)\s.*/$workingFolder|$1/;
			$sequence = <READS>; chomp($sequence);
			#print "id = $idQuery => $sequence\n";

			my $tempFile = "/dev/shm/temp$workingFolder.fa";
			open (TEMP, ">$tempFile") || die ("cannot use temp file: $tempFile; it might be used\n");
			print TEMP "$idQuery\n";
			print TEMP "$sequence\n";
			close (TEMP);

			
			#print "$pidEe=>megablast -d /dev/shm/$workingFolder/mRNA_hom.fna -i $tempFile -m 8 -P 99 -a2 -v1 | head -n1 | awk \'{print \$2}\' \n";
			my $idSubject = `megablast -d /dev/shm/$workingFolder/mRNA_hom.fna -i $tempFile -m 8 -P 99 -a2 -v1 | head -n1 | awk \'{print \$2}\'`;
			chomp $idSubject;
			if ($idSubject ne "") {
				#print "idsubject = $idSubject\n";
				if ($resultsHash{$idSubject}) {
					$resultsHash{$idSubject}++;
				} else {
					$resultsHash{$idSubject} = 1;
				}
			}
		} else {
			die "there is a problem with the format of $readsFilename\n";
		}
#		$count++;
		#print "$count\n";
#		if ($count == 1000) {
#			print "reached limit point, closing loop\n";
#			last;
#		}
	}
	close (READS);
	foreach my $key (keys %resultsHash) {
		print RESULT ("$key\t$resultsHash{$key}\n");
	}
	close (RESULT);	
	chdir("..") || die ("$pidEe=>cannot come back to origin (previous) path\n");
	#end of the multithreads
	$pmEe->finish;
}
#checks if there is nothing else running
$pmEe->wait_all_children;
}

###############################################################################





###############################################################################
###############################################################################
print "done\n";
exit 0;


sub INT_handler {
	#`rm -r *vsAll_* 2> /dev/null`;
	foreach (@folders) {
		`rm -r /dev/shm/$_ 2> /dev/null`;
#		`rm /dev/shm/temp$_.fa 2> /dev/null`;
	#	`rm -r $_/mRNA_* 2> /dev/null`;
	}
	`clear`;
	exit(0);

}

