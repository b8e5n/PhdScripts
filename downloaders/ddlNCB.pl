use LWP::Simple;
use Parallel::ForkManager;
my $max_process = 16;

my $argSize = $#ARGV+1;
#print "there is $argSize arguments\n\n";
if ($argSize != 3) {
	print "Usage: perl ddlNCBI.pl [ncbi query] [data type] [output filename]
		The query can be previouly performed at \"http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi\" for tests.
		Query example:\t\"txid7148[Organism:exp]\"
		Data type:\t\"nucleotide\" or \"protein\"
		For outfile:\tyou can use either full path or relative path.\n";
	exit 0;
	}
my $query = @ARGV[0];
my $type = @ARGV[1];
my $outFile = @ARGV[2];
my $outFileTemp =$outFile."Temp";

#assemble the esearch URL
$base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
$url = $base . "esearch.fcgi?db=$type&term=$query&usehistory=y";


#post the esearch URL
$output = get($url);


#parse WebEnv, QueryKey and Count (# records retrieved)
$web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
$key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
$count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);


#open output file for writing
open(OUTTEMP, ">$outFileTemp") || die "Can't open file!\n";
open(OUT, ">$outFile") || die "Can't open file!\n";

#retrieve data in batches of 500
$retmax = 1000;
my $pm = new Parallel::ForkManager($max_process);
for ($retstart = 0; $retstart < $count; $retstart += $retmax) {
	my $pid = $pm->start and next;

        $efetch_url = $base ."efetch.fcgi?db=$type&WebEnv=$web";
        $efetch_url .= "&query_key=$key&retstart=$retstart";
        $efetch_url .= "&retmax=$retmax&rettype=fasta&retmode=text";
        $efetch_out = get($efetch_url);
	flock(OUTTEMP, LOCK_EX);
        print OUTTEMP "$efetch_out";
	flock(OUTTEMP, LOCK_UN);
	$pm->finish;
}
$pm->wait_all_children;
close OUTTEMP;
print "Downloading Finished, cleaning from empty sequences...\n";
`awk \'\$2\{print RS\}\$2\' FS=\'\\n\' RS=\\> ORS= $outFileTemp > $outFile`;
print "Cleaning Finished, removing temp file...\n";
`rm $outFileTemp`;

print "done\n";
exit 0;
