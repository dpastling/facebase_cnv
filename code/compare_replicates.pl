#!/usr/bin/env perl

use strict;
use warnings;

my $folder = "results/vanilla_ice_qc";

my %samples = (
	"18508" => [ "1064691934", "1064692018" ],
	"18858" => [ "1064691581", "1064691639" ],
	"18861" => [ "1064691574", "1064691968" ],
	"18502" => [ "1054736219", "1054736256" ], 
	"18505" => [ "1054736194", "1054736220", "1064692236" ],
	"06985" => [ "1054736171", "1054736250", "1064733235" ]
);

my @method = ("vice", "penn", "dna");

foreach my $key (keys %samples)
{
	my $sample     = $key;
	my $replicate1 = $samples{$key}[0];
	my $replicate2 = $samples{$key}[1];
	my $file1      = "$folder/$sample.$replicate1";
	my $file2      = "$folder/$sample.$replicate2";
	my $output12   = "$folder/$sample\_$replicate1\_$replicate2";
	foreach my $suffix (@method)
	{
		intersect("$file1\_$suffix.bed", "$file2\_$suffix.bed", "$output12\_$suffix.bed");
	} 

	if ($#{ $samples{$key} } > 1)
	{
		my $replicate3 = $samples{$key}[2];
		my $file3 = "$folder/$sample.$replicate3";
		my $output23   = "$folder/$sample\_$replicate2\_$replicate3";
		my $output13   = "$folder/$sample\_$replicate1\_$replicate3";
		foreach my $suffix (@method)
		{
			intersect("$file2\_$suffix.bed", "$file3\_$suffix.bed", "$output23\_$suffix.bed");
			intersect("$file1\_$suffix.bed", "$file3\_$suffix.bed", "$output13\_$suffix.bed");
			intersect("$output12\_$suffix.bed", "$file3\_$suffix.bed", "$sample\_all\_$suffix.bed");
		} 
	} else {
		foreach my $suffix (@method)
                {
                        system("cp $output12\_$suffix.bed $sample\_all\_$suffix.bed");
                }
	}
}

sub intersect {
	my $file1  = $_[0];
	my $file2  = $_[1];
	my $output = $_[2];
	# find overlaps and make sure the logR Ratios are in same direction
	my $vice_command = "bedtools intersect -wo -a $file1 -b $file2 | ";
	$vice_command .= "awk '(\$11 > 0.1 && \$23 > 0.1) || (\$11 < -0.1 && \$23 < -0.1)' | ";
	$vice_command .= "cut -f 1-12 > $output";
	system($vice_command);
	# do a proper intersect, hopefully trimming regions
	#system("bedtools intersect -a temp.bed -b file2 > $output");
	#system("rm temp.bed");
}

