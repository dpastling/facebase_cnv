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

my $dup_threshold = 0.05;
my $del_threshold = -0.1;

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
		split_file("$file1\_$suffix.bed");
		split_file("$file2\_$suffix.bed");
		system("bedtools intersect -u -a $file1\_$suffix\_dup.bed -b $file2\_$suffix\_dup.bed > $output12\_$suffix\_dup.bed");
		system("bedtools intersect -u -a $file1\_$suffix\_del.bed -b $file2\_$suffix\_del.bed > $output12\_$suffix\_del.bed");
	} 

	if ($#{ $samples{$key} } > 1)
	{
		my $replicate3 = $samples{$key}[2];
		my $file3 = "$folder/$sample.$replicate3";
		my $output23   = "$folder/$sample\_$replicate2\_$replicate3";
		my $output13   = "$folder/$sample\_$replicate1\_$replicate3";
		foreach my $suffix (@method)
		{
			split_file("$file3\_$suffix.bed");
			system("bedtools intersect -u -a $output12\_$suffix\_dup.bed -b $file3\_$suffix\_dup.bed > $folder/$sample\_all_$suffix\_dup.bed");
			system("bedtools intersect -u -a $output12\_$suffix\_del.bed -b $file3\_$suffix\_del.bed > $folder/$sample\_all_$suffix\_del.bed");

			system("bedtools intersect -u -a $file1\_$suffix\_dup.bed -b $file3\_$suffix\_dup.bed > $output13\_$suffix\_dup.bed");
			system("bedtools intersect -u -a $file1\_$suffix\_del.bed -b $file3\_$suffix\_del.bed > $output13\_$suffix\_del.bed");
				
			system("bedtools intersect -u -a $file2\_$suffix\_dup.bed -b $file3\_$suffix\_dup.bed > $output23\_$suffix\_dup.bed");
			system("bedtools intersect -u -a $file2\_$suffix\_del.bed -b $file3\_$suffix\_del.bed > $output23\_$suffix\_del.bed");
		} 
	} else {
		# If there are only two samples then the intersection is all there is
		foreach my $suffix (@method)
                {
                        system("cp $output12\_$suffix\_dup.bed $folder/$sample\_all\_$suffix\_dup.bed");
                        system("cp $output12\_$suffix\_del.bed $folder/$sample\_all\_$suffix\_del.bed");
                }
	}
}

sub split_file {
	my $file = $_[0];
	my $dup = $file;
	my $del = $file;
	$dup =~ s/\.bed/_dup.bed/;
        $del =~ s/\.bed/_del.bed/;
        system("awk '\$11 > $dup_threshold' $file > $dup");
        system("awk '\$11 < $del_threshold' $file > $del");
}


