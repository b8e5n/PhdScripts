#!/bin/perl -w
use strict;
use diagnostics;
use Parallel::ForkManager;
use Term::ProgressBar 2.00;
use Fcntl qw(:flock SEEK_END); # import LOCK_* and SEEK_END constants
$SIG{'INT'} = 'INT_handler';# manage the killing signal => remove temp files :D
my $max_process = 8;

#AUTHOR BENOIT CARRERES
#multi threads is set to 8
# for final table :"Probe.Set.Name\tSubject id\ts. start\ts. end\tmissmatch\tgaps\tPM Candidates(single multiple)".
# if no response\t the 4 last columns will be empty

my $argSize = $#ARGV+1;
print "there is $argSize arguments\n\n";

my $originPath = `pwd`;
chomp $originPath;
print "the origin folder is : $originPath\n";
my @folders = split("\n",`ls -l | egrep '^d' | awk \'{print \$9}\'`);


#separateCandidates();
#getherSequencesHomo();

#blastHomoWithReads();
#calculateExpressionLevels();
#generateComparativeExpressionFile();
#blastHomoWithNCBI_DB();#unfinished
#mRNAfastaToSingleLinedSequences(); # to be ran before mRNAsMerge
#mRNAsMerge(); #takes all hom.fa, with single lined sequences and merge fatas without replication.
#blastGoldenGenes(); #for now blast only against proteins
selectUniqueCandidates();
####getherSequencesHetero();


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
print "==================Genes/Reads matching================================\n";
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
	my $currentPath = `pwd`; chomp $currentPath;
	print "$pidEe=>the current folder is : $currentPath\n";
	
	my %resultsHash;
	my $count = 0;
	my $loop = 0;
	my $totalMatching = 0;
	my $totalReads = `wc -l ../$readsFilename | awk \'{print \$1}\'`/2;
	my $totalProcessed = 0;
	my $outputFilename = "homMatchingReads";
	while (my $idQuery = <READS>){
		chomp $idQuery;
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

			
			#print "$pidEe=>megablast -d /dev/shm/$workingFolder/mRNA_hom.fna -i $tempFile -m 8 -P 99 -a2 -v1 | head -n1 | cut -f2 \n";
			my $idSubject = `megablast -d /dev/shm/$workingFolder/mRNA_hom.fna -i $tempFile -m 8 -P 99 -a1 -v1 | head -n1 | cut -f2`;
			chomp $idSubject;
			if ($idSubject ne "") {
				#print "idsubject = $idSubject\n";
				$totalMatching++;
				if ($resultsHash{$idSubject}) {
					$resultsHash{$idSubject}++;
				} else {
					$resultsHash{$idSubject} = 1;
				}
			}
		} else {
			die "there is a problem with the format of $readsFilename\n";
		}
		$count++;
		$totalProcessed++;
		my $limit = 50000;
		#print "$count\n";
		if ($count == $limit) {
			$loop++;
			my $blastedReads = $limit*$loop;
			print "reached limit point, writing file at this level\n";
#			last;
			#intermediary result printing
			printHashFile($outputFilename, $currentPath, $totalReads, $totalMatching, $blastedReads, %resultsHash);
			$count = 0;
			my $percentage = ($blastedReads / $totalReads) * 100;
			print "$workingFolder=>".sprintf("%.3f %%", $percentage)."\n";
		}
	}
	close (READS);
	#print results
	printHashFile($outputFilename, $currentPath, $totalReads, $totalMatching, $totalProcessed, %resultsHash);
	chdir("..") || die ("$pidEe=>cannot come back to origin (previous) path\n");
	#end of the multithreads
	$pmEe->finish;
}
#checks if there is nothing else running
$pmEe->wait_all_children;
}

###############################################################################



###############################################################################
################Recalculate expressions from matching levels###################
###############################################################################
sub calculateExpressionLevels {
print "==================Expression Level Calculation===========================\n";
my $pmEe = new Parallel::ForkManager($max_process);

foreach my $workingFolder (@folders) {
	#start multithreading
	my $pidEe = $pmEe->start and next;
	#now several jobs
	print "=>$workingFolder\n";
	chdir($workingFolder) || die ("$pidEe=>cannot got to directory : $_");
	my $currentPath = `pwd`;
	print "$pidEe=>the current folder is : $currentPath\n";

	my $outputFilename = "homReadsMatchingScore";

	open (MATCHES, "<homMatchingReads") || die ("$workingFolder=>pb reading genes-reads-matches file homMatchingReads\n");
	my $totalReads;
	my $totalMatching;
	my $blastedReads;
	my %scoreHash;
	
	while (my $line =<MATCHES>) {
		chomp $line;
		my @splitedLine = split("\t", $line);
		my $id = $splitedLine[0];
		#print "line = @splitedLine\n";
		if ($id eq "TotalReads") {
			print "matching $id\n";
			$totalReads = $splitedLine[1];
		}elsif ($id eq "TotalMatchingReads") {
			print "matching $id\n";
			$totalMatching = $splitedLine[1];
		}elsif ($id eq "TotalBlastedReads") {
			print "matching $id\n";
			$blastedReads = $splitedLine[1];
		}else {
############################### testing score calculation
#			my $testNumber = 36655;
#			my $testScore = ($testNumber*$totalMatching)/($blastedReads*$blastedReads);
#			print "testScore $workingFolder = $testScore\n";
#			$scoreHash{$testNumber} = $testScore;
#			last;
###############################

			my $matchNumber = $splitedLine[1];
			my $score = ($matchNumber*$totalMatching)/($blastedReads*$blastedReads);
			$score = sprintf("%.9f", $score);
			$scoreHash{$id} = $score;
		}
	}
	close (MATCHES);
	chdir("..") || die ("$pidEe=>cannot come back to origin (previous) path\n");
	printHashFile($outputFilename, $workingFolder, $totalReads, $totalMatching, $blastedReads, %scoreHash);
	
	#end of the multithreads
	$pmEe->finish;
}
#checks if there is nothing else running
$pmEe->wait_all_children;
}
###############################################################################


###############################################################################
################Regenerate comparative file with expressions###################
###############################################################################
sub generateComparativeExpressionFile {
print "================Generate Comparative expression file======================\n";
my $pmEe = new Parallel::ForkManager($max_process);

my %expressionFiles;
my $expressionFileName = "homReadsMatchingScore";

foreach my $workingFolder (@folders) {
	my %expressionFile;
	chdir($workingFolder) || die ("$workingFolder=>cannot got to directory : $_");
	open (PEL, "<$expressionFileName") || die ("$workingFolder=>pb reading expression file $expressionFileName\n");
	while (my $line =<PEL>) {
		chomp $line;
		my @splitedLine = split("\t", $line);
		my $id = $splitedLine[0];
		#print "line = @splitedLine\n";
		if ($id ne "TotalReads" && $id ne "TotalMatchingReads" && $id ne "TotalBlastedReads") {
			$expressionFile{$id} = $splitedLine[1];
		}
	}#endwhile
	chdir("..") || die ("$workingFolder=>cannot come back to origin (previous) path\n");
	close(PEL);
	$expressionFiles{$workingFolder} = \%expressionFile;
}#endForeach

foreach my $workingFolder (@folders) {
	if ($workingFolder =~ /^\d.*/) {##################################limiting folders

	#multithreading start
	my $pidEe = $pmEe->start and next;
	#now several jobs
	print "=>$workingFolder\n";
	my $comparativeFileName = $workingFolder."vsAll_hom";
	
	my $comparativeExpressionFileName = $comparativeFileName."_PEL";

	open (COMP, "<$comparativeFileName") || die ("$workingFolder=>pb reading comparison file $comparativeFileName\n");
	open (COMPEX, ">$comparativeExpressionFileName") || die ("$workingFolder=>pb writing comparison expression file $comparativeExpressionFileName\n");
	while (my $line = <COMP>) {
		my @lineContent = split("\t",$line);
#		print "tableSize = $#lineContent\n";
		my $separation = "";
		foreach my $column (@lineContent) {
			if ($column =~ /^([\d\w]*?)[|](.*)/) {
				if ($expressionFiles{$1}{$2}) {
					#print "ID:$2\t= $column\t==> $expressionFiles{$1}{$2}, found\n--> @lineContent\n";
					print COMPEX $separation."$column,$expressionFiles{$1}{$2}";
				}# else {
				#	print "gene sample ID:$2\t=> $column was not found, verify the ID.\n";
				#}
			}
			elsif ($column eq "NA") {
				print COMPEX $separation."$column,0";
			}
			$separation = "\t";
		}
		print COMPEX "\n";
	}
	close(COMP);
	close(COMPEX);
	#end of the multithreads
	$pmEe->finish;
	}##################################limiting folders
}
#checks if there is nothing else running
$pmEe->wait_all_children;
}
###############################################################################
###############################################################################


###############################################################################
########################Align/count reads on candidates########################
###############################################################################
sub blastHomoWithNCBI_DB {
print "==================Genes/Reads matching================================\n";
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
	
	
	chdir($workingFolder) || die ("$pidEe=>cannot got to directory : $_");
	my $currentPath = `pwd`;
	print "$pidEe=>the current folder is : $currentPath\n";

	my $mRNAsFilename = "mRNA_hom.fna";
	
	my $outputFilename = "homMappingNcbi";
	open (MAPPING, ">$outputFilename") || die ("$pidEe=>cannot open file to write mapping: $outputFilename\n");
	print MAPPING "contigID\tChlamydomonas\tArabidopsis\tVitis\n";
	open (CANDIDATES, "<homMatchingReads");
	my $tempFile = "/dev/shm/temp$workingFolder.fa";
	while (my $idQuery = <CANDIDATES>){
		chomp $idQuery;
		my @line = split("\t", $idQuery);
		$idQuery = $line[0];
		if ($idQuery ne "TotalReads" && $idQuery ne "TotalMatchingReads" && $idQuery ne "TotalBlastedReads") {
			`fastacmd -d $mRNAsFilename -s $idQuery > $tempFile`;
			my $idVvi = `megablast -d ../Vvi.seq.all -i $tempFile -m 8 -P 99 -a2 -v1 | head -n1 | cut -f2`;
			my $idAt = `megablast -d ../At.seq.all -i $tempFile -m 8 -P 99 -a2 -v1 | head -n1 | cut -f2`;
			my $idCr = `megablast -d ../Cre.seq.all -i $tempFile -m 8 -P 99 -a2 -v1 | head -n1 | cut -f2`;
			chomp $idVvi;chomp $idAt;chomp $idCr;
			if ($idVvi ne "" || $idCr ne "" || $idAt ne "") {
				print MAPPING "$idQuery\t$idCr\t$idAt\t$idVvi\n";
			}
		}		
	}
	close (MAPPING);
	close (CANDIDATES);
	chdir("..") || die ("$pidEe=>cannot come back to origin (previous) path\n");
	#end of the multithreads
	$pmEe->finish;
}
#checks if there is nothing else running
$pmEe->wait_all_children;
}

###############################################################################


###############################################################################
############################Blast the golden genes#############################
###############################################################################
sub blastGoldenGenes {
my $pmEe = new Parallel::ForkManager($max_process);

my @organismsToBlast = ("Cre_uniprot.fasta", "Ath_uniprot.fasta", "Sac_uniprot.fasta", "Ptr_uniprot.fasta");

open(RNA, "<mRNA_hom_gd.fa") || die ("cannot open file to read mapping: mRNA_hom_gd.fa\n");
open(my $SIXTY, ">60.prot.mapping") || die ("cannot open file for writing\n");
open(my $FIFTY, ">50.prot.mapping") || die ("cannot open file for writing\n");
open(my $THURTYFIVE, ">35.prot.mapping") || die ("cannot open file for writing\n");

my $counter = 1;
`mkdir /dev/shm/mrnas`;
my $dbFolder = "/dev/shm/dbs";
`mkdir $dbFolder`;
`cp  /home/benoit/dbs/prots/* $dbFolder`;
`chmod -u+rw $dbFolder/*`;

while (my $id = <RNA>) {
	chomp $id;
	if ($id =~ /^>/) {
		my $sequence = <RNA>; chomp $sequence;
		my $tempFasta = ">/dev/shm/mrnas/temp$counter.fa";
		open (TEMPFILE, $tempFasta) || die ("cannot create temp fasta file in ramdisk : \"$tempFasta\"\n");
		print TEMPFILE "$id\n$sequence";
		close(TEMPFILE);
	} else {
		die("line not matching... pb in the file format, should be fasta ID:\n $id\n");
	}
	$counter++;
}
my @tempFiles = split("\n",`ls -l /dev/shm/mrnas/ | egrep \'^-\' | awk \'{print \$9}\' | egrep \'temp[[:digit:]]+.fa\$\'`);
my $progress = Term::ProgressBar->new({name  => 'BlastsX',
                                         count => $#tempFiles,
                                         ETA   => 'linear', });
$progress->max_update_rate(5);
$counter = 0;
my $next_update = 0;
foreach my $fileName (@tempFiles) {
	chomp $fileName;
	$counter++;
	#start multithreading
	my $pidEe = $pmEe->start and next;
	#now several jobs
	$fileName = "/dev/shm/mrnas/".$fileName;
#	print "blasting 60% $fileName\n";
	foreach my $organismDB (@organismsToBlast) {
		my $seqSize = `egrep -v \'^>\' $fileName | wc -m`;chomp $seqSize;
		my $protSize = $seqSize/3;
		my $protPlus = $protSize+(20*$protSize/100);
		my $protMinus = $protSize-(20*$protSize/100);
		my @blast60 = split ("\n", `blast2 -p blastx -d $dbFolder/$organismDB -i $fileName -m8 -P60 -a2`);
		foreach my $blastResult60 (@blast60) {
			my @blastResultT60 = split("\t", $blastResult60);
			if ($blastResultT60[3] >= $protMinus && $blastResultT60[3] <= $protPlus) {
				lock($SIXTY);
				print $SIXTY "$blastResult60\n";
				unlock($SIXTY);
			}
		}
		my @blast50 = split ("\n", `blast2 -p blastx -d $dbFolder/$organismDB -i $fileName -m8 -P50 -a2`);
		foreach my $blastResult50 (@blast50) {
			my @blastResultT50 = split("\t", $blastResult50);
			if ($blastResultT50[3] >= $protMinus && $blastResultT50[3] <= $protPlus) {
				lock($FIFTY);
				print $FIFTY "$blastResult50\n";
				unlock($FIFTY);
			}
		}
		my @blast35 = split ("\n", `blast2 -p blastx -d $dbFolder/$organismDB -i $fileName -m8 -P35 -a2`);
		foreach my $blastResult35 (@blast35) {
			my @blastResultT35 = split("\t", $blastResult35);
			if ($blastResultT35[3] >= $protMinus && $blastResultT35[3] <= $protPlus) {
				lock($THURTYFIVE);
				print $THURTYFIVE "$blastResult35\n";
				unlock($THURTYFIVE);
			}
		}
		
	}
	#progressbar
	$next_update = $progress->update($counter) if $counter >= $next_update;
	#end of the multithreads
	$pmEe->finish;
}
#checks if there is nothing else running
$pmEe->wait_all_children;
close($SIXTY);
close($FIFTY);
close($THURTYFIVE);
}
###############################################################################

###############################################################################
############################Blast the golden genes#############################
###############################################################################
sub selectUniqueCandidates {
open(SIXTY, "<60.prot.mapping") || die ("cannot open file for reading\n");
open(my $FIFTY, "<50.prot.mapping") || die ("cannot open file for reading\n");
open(my $THURTYFIVE, "<35.prot.mapping") || die ("cannot open file for reading\n");
my %sixty;
my %fifty;
my %thurtyfive;

while (my $line60 = <SIXTY>) {
	my @lineTab = split("\t", $line60);
	my $id = shift(@lineTab);
	print "$lineTab[0]\n";
	my $dbType = $lineTab[0];
	if ($sixty{$id}) {
		$sixty{$id} = \@lineTab;
	} else {
		my @previousTab = $sixty{$id};
		my $previousDbType = $previousTab[0];
		my $score = $lineTab[-1];
		my $previousScore = $previousTab[-1];
		if ($previousDbType =~ /^tr/ && $dbType =~ /^sp/) {
			$sixty{$id} = \@lineTab;
		} elsif ($previousDbType =~ /^tr/ && $dbType =~ /^tr/ || $previousDbType =~ /^sp/ && $dbType =~ /^sp/) {
			if ($score>$previousScore) {
				$sixty{$id} = \@lineTab;
			}
		}
	}
}
printHashFile("60.prot.mapping.uniq", "./", "", "", "", %sixty);
close(SIXTY);
close($FIFTY);
close($THURTYFIVE);
}

print "done\n";
exit 0;

sub printHashFile {
	my ($fileName, $workingFolder, $totalReads, $totalMatching, $blastedReads, %resultsHash) = @_;
	
	open (RESULT, ">$workingFolder/$fileName") || die ("$workingFolder=>pb writing-creating expression file $workingFolder/$fileName\n");
	if ($totalReads != "" && $totalMatching != "" && $blastedReads != "") {
		print RESULT ("TotalReads\t$totalReads\n");
		print RESULT ("TotalMatchingReads\t$totalMatching\n");
		print RESULT ("TotalBlastedReads\t$blastedReads\n");
	}
	foreach my $key (keys %resultsHash) {
		print RESULT ("$key\t$resultsHash{$key}\n");
	}
	close (RESULT);
}
sub mRNAfastaToSingleLinedSequences {
	`for file in \`ls -l | egrep '^d' | awk \'{print \$9}\' \`; do fasta_formatter -i \$file/mRNA_hom.fna > \$file/mRNA_hom.fa; formatdb -i \$file/mRNA_hom.fa -p F -o; done`;
}
sub mRNAsMerge {
	`cat */mRNA_hom.fa | fastx_collapser -o mRNA_hom_gd.fa`
}
######################### locking ###############################
sub lock {
	my ($fh) = @_;
	flock($fh, LOCK_EX) or die "Cannot lock file - $!\n";
	# and, in case someone appended while we were waiting...
	seek($fh, 0, SEEK_END) or die "Cannot seek - $!\n";
}
sub unlock {
	my ($fh) = @_;
	flock($fh, LOCK_UN) or die "Cannot unlock file - $!\n";
}
######################### ctrl-c ################################
sub INT_handler {
	#`rm -r *vsAll_* 2> /dev/null`;
	foreach (@folders) {
#		`rm -r /dev/shm/$_ 2> /dev/null`;
#		`rm /dev/shm/temp$_.fa 2> /dev/null`;
	#	`rm -r $_/mRNA_* 2> /dev/null`;
	}
	`clear`;
	exit(0);

}
