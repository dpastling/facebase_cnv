#!/usr/bin/env perl

use strict;
use warnings;

my $genome_ref = "/vol4/home/astlingd/genomes/bowtie2.2.5_indicies/hg19/hg19.fa";
my $window_size = 1000000;
my $window_half = $window_size / 2;

my %chromosome_size;
my $chrom;
open(GENOME, $genome_ref) or die "cannot find genome reference: $genome_ref";
while(<GENOME>)
{
	chomp;
	my $line = $_;
	if ($line =~ /^>/)
	{
		$chrom = $line;
		$chrom =~ s/^>//;
		next;
	}
	$chromosome_size{$chrom} += length($line);
}
close(GENOME);

print join("\t", "Name", "Chr", "Position", "GC") . "\n";
while(<>)
{
	chomp;
	my $line = $_;
	my ($chr, $pos1, $pos2, $name) = split("\t", $line);
	my $position = $pos2;
	if ($position - $window_half < 1)
	{
		$pos1 = 1;
		$pos2 = $window_size;
	} elsif ($position + $window_half > $chromosome_size{$chr})
	{
		$pos2 = $chromosome_size{$chr};
		$pos1 = $pos2 - $window_half;
	} else 
	{
		$pos1 = $pos2 - $window_half;
		$pos2 = $pos2 + $window_half;
	}
	my @sequence = `samtools faidx $genome_ref $chr:$pos1-$pos2`;
	my $gc = 0;
	for (my $i = 1; $i <= $#sequence; $i++)
	{
		my $seq = $sequence[$i];
		chomp $seq;
		$seq =~ s/[ATN]//g;
		$gc += length($seq);
	}
	my $gc_percent = ($gc / $window_size) * 100;
	my $result = join("\t", $name, $chr, $position, $gc_percent);
	print "$result\n";
}


