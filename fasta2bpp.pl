use warnings;
use strict;

my ($indir,$list) = @ARGV;
my %seq;
open LI,$list or die "Can't open $list $!\n";
while (<LI>){
	chomp;
	my $gene = $_;
	#print"$_\n";
	open IN,"$indir/$gene.fa" or die "Can't open $gene $!\n";
	my $name;
	while (<IN>){
		chomp;
		if (/^>(.*)/){
			$name = $1;
		}
		else{
			#s/N/\?/g;
			$seq{$gene}{"$gene^$name"} .= $_;
		}
		#print"$_\n";
	}
	close IN;
}
close LI;

my $num_gene = scalar(keys %seq);
my $sum_len;
my %length;
my $loci;
foreach my $key (sort keys %seq){
	$loci .= "$key,";
	foreach my $ind (sort keys %{$seq{$key}}){
		#print"$ind\t$seq{$key}{$ind}\n";
		$sum_len += length($seq{$key}{$ind});
		$length{$key} = length($seq{$key}{$ind});
		last;
	}
}
$loci =~ s/,$//;
my $string = <<PLOT;
PLOT

foreach my $key (sort keys %seq){
	#$string .= "\n$key $length{$key}\n";
	if($length{$key} >=10000)
	{
		$string .= "\n8 10000\n";
	}
	else
	{
	$string .= "\n8 $length{$key}\n";}
	foreach my $ind (sort keys %{$seq{$key}}){
		#$ind=~/^((...).*?)(_.*)/;
		#print"$ind\t$seq{$key}{$ind}\n";
		#$ind=~/^((...).*?)_(.*)/;
		#$ind=~/^((...).*?)_(.*)/;
		#my($gene,$sp,$id)=($1,$2,$3);
		#$gene=~s/\.1//g;
		my $seq="";
		if($length{$key} >= 10000)
		{
			$seq=substr($seq{$key}{$ind},0,10000);
		}
		else
		{
			$seq=$seq{$key}{$ind};
		}
		#$string .= "$gene^$sp-$id\t$seq\n";
		$string .= "$ind\t$seq\n";
	}
}

$string .= <<PLOT;
PLOT
print "$string";


