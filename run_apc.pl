use strict;
use warnings;

my $JSON = '
{
	"gff_introns": false,
	
	"cli": {
		"--min_exon": 25,
		"--min_intron": 35,
		"--max_splice": 3,
		"--flank": 99,
		"--limit": 100,
	
		"--apwm": "data/acceptor.pwm",
		"--dpwm": "data/donor.pwm",
		"--emm" : "data/exon.mm",
		"--imm" : "data/intron.mm",
		"--elen": "data/exon.len",
		"--ilen": "data/intron.len"
	},
	
	"data": [
		{
			"name": "CHROMOSOME",
			"fasta": "APCROOT/CHROMOSOME.fa",
			"gff": "APCROOT/CHROMOSOME.gff3"
		}
	]
	
}
';

die "usage: $0 <apc.genes.txt> <apc dir> <cpus>" unless @ARGV == 3;
my ($genes, $apc, $cpus) = @ARGV;

open(my $fh, $genes) or die;
my $header = <$fh>;
while (<$fh>) {
	my ($chr) = split;
	my $json = $JSON;
	$json =~ s/CHROMOSOME/$chr/g;
	$json =~ s/APCROOT/$apc/g;
	open(my $ofh, ">tmp.json") or die;
	print $ofh $json;
	close $ofh;
	print "$chr\t";
	system("./optiso --program isoformer --cpu $cpus tmp.json") == 0 or die; 
}
