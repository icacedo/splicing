use strict;
use warnings;

my $JSON = '
{
	"gff_introns": true,
	
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
			"fasta": "data/apc/CHROMOSOME.fa",
			"gff": "data/apc/CHROMOSOME.gff3"
		}
	]
	
}
';


open(my $fh, "data/719.txt") or die;
my $header = <$fh>;
while (<$fh>) {
	my ($chr) = split;
	my $json = $JSON;
	$json =~ s/CHROMOSOME/$chr/g;
	open(my $ofh, ">tmp.json") or die;
	print $ofh $json;
	close $ofh;
	system("./optiso --pop 10 --gen 10 tmp.json") == 0 or die; 
}
