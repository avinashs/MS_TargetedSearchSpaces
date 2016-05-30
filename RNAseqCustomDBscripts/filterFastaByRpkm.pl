#Avinash Shanmugam, March 12th 2012
#The script will accept as input a Protein fasta database,
#a file containing rpkm values for the proteins (from rnaseq)
#and the rpkm value to be used as a threshold
#The script creates a list of protein ids above the given threshold
#and extracts the seq records of only those proteins and writes
#them into another fasta file

sub getProtList
{
	my $rpkmfile = shift @_;
	my $rpkmThreshold = shift @_;

	open(FILE,"<$rpkmfile")||die("Error.. can't open $rpkmfile!\n");
	
	#Read and discard title line
	my $line = <FILE>;

	#Declare an empty array to store the selected protIds
	my @protList;

	while($line = <FILE>)
	{
		chomp $line;

		my @fields = split(/\t/,$line);

		my $protId= $fields[0];
		my $rpkm = $fields[2];

		if( $protId ne 'null' && $rpkm > $rpkmThreshold)
		{
			push(@protList, $protId); 
		}
	}

	return @protList;

}




################################# Main ##############################

if(scalar(@ARGV) < 4)
{
	print "USAGE: filterFastaByRpkm.pl <in fasta> <rpkm file> <rpkm threshold> <out fasta>\n\n";
	exit;
}

my $infasta = shift;
my $rpkmfile = shift;
my $rpkmThreshold = shift;
my $outfasta = shift;

my @protList = getProtList($rpkmfile, $rpkmThreshold);
print "No of proteins passing threshold =".scalar(@protList)."\n\n";

#Convert the protList array into a hash for faster lookup
my %protHash = map{ $_ => 1} @protList;

##Start to look through the fasta file

#Open in and out fasta
open(IN,"<$infasta")||die("Error.. can't open $infasta!\n");
open(OUT,">$outfasta")||die("Error.. can't open $outfasta!\n");

#Find total no. of seqs in the input fasta file
$seq_tot = `grep -c ">" $infasta`;
chomp $seq_tot;

#Declare variables needed in the loop below
my $seq_no = 1;
my $seq = "";
my $protid;
my $titleLine;

while( $line = <IN>)
{
	#If line is a title line
	if( $line=~m/>/)
	{
		#From second title line onwards
		if( $seq_no > 1)
		{
			if( exists $protHash{$protid})
			{
				print OUT $titleLine.$seq;
			}

			print "$seq_no of $seq_tot processed..\r";

		}

		#Store line variable as titleLine
		$titleLine = $line;

		#Split the fields in title line and extract protid
		my @fields = split(/\s/,$titleLine);
		$protid = substr($fields[0],1);

		#Increment seq_no counter
		$seq_no = $seq_no + 1;

		#Set $seq to empty string, ready to store next seq
		$seq = "";
	}
	#If not title line, it is a seq line.
	#Keep concatenating these lines into the $seq variable
	else
	{
		$seq = $seq.$line;
	}
}

print "\n\nFiltered fasta file written to $outfasta\n\n";
