#!/usr/bin/perl -w
use strict;

my $output_sword=shift;




my %hash_junctions;
my %hash_weighted_junctions;
my $total=0;
print "#Hinge/Junction consistency:\n";

open(I,"$output_sword") or die "Cannot open file output sword \"$output_sword\" : $!\n";

WHILE:while(my $line=<I>)
{
    if ($line =~/^(\d{1,})\s{0,}\|(\d{0,})\s{0,}\|\s{0,}(.*)\|\s{0,}(\d{1,}\.\d{1,7})\|\s{0,}(.*)\|/)
    {
        my $ndomains=$1;
        my $min=$2;
        my $delineation=$3;
        my $average_kappa=$4;
        my $quality=$5;
        next WHILE if ($quality =~/n/);
        my $cnt = () =  $quality =~ /(\*)/g ;
        #print "$ndomains $min $delineation $average_kappa $quality $cnt\n";
        $total++;
        my @tab_domains=split(/\s/,$delineation);
        foreach my $domain (@tab_domains)
        {
            if ($domain=~/^(\d{1,})-/)
            {
                my $junction=$1;
                if (not exists $hash_junctions{$junction})
                {
                    $hash_junctions{$junction}=1;
                    $hash_weighted_junctions{$junction}=$cnt;
                }
                else
                {
                    $hash_junctions{$junction}++;
                    $hash_weighted_junctions{$junction}=$hash_weighted_junctions{$junction}+$cnt;
                }
            }

            if ($domain=~/-(\d{1,})$/)
            {
                my $junction=$1;
                if (not exists $hash_junctions{$junction})
                {
                    $hash_junctions{$junction}=1;
                    $hash_weighted_junctions{$junction}=$cnt;
                }
                else
                {
                    $hash_junctions{$junction}++;
                    $hash_weighted_junctions{$junction}=$hash_weighted_junctions{$junction}+$cnt;
                }
            }
        }
    }
}
close I;


printf("#%-4s %4s %4s %4s\n","Jnct","Cnt","Raw","Wei");
foreach my $junction (sort {$a <=> $b} keys %hash_junctions) 
{
    printf("%-6d %-4d %4.2f %4.2f\n",$junction,$hash_junctions{$junction},$hash_junctions{$junction}/$total,$hash_weighted_junctions{$junction}/$total);
}



