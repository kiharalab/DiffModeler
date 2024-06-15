#!/usr/bin/env perl

use strict;
use warnings;

#my @results = `cat $ARGV[0]`;
#my @results = `./SWORDdisplay -i 1jx4 -c A -d --max 9 --path foo/`;

my $cmd = $ARGV[0];

my @results = `$cmd`;

my $i       = 0 ; #JCG
my $dom_opt = 0 ; #JCG
# MODIF JCG more
FORI:for ($i =0 ; $i <= $#results ; $i++)
{
    if ($results[$i] =~ /^(\d{1,})/)
    {
        $dom_opt = $1;
        last FORI;
    }

}

my $domP1 = $dom_opt+1;
my $domP2 = $dom_opt+2;
my $domM1 = $dom_opt-1;
my $domM2 = $dom_opt-2;

my $countO = 0;
my $countP1 = 0;
my $countP2 = 0;
my $countM1 = 0;
my $countM2 = 0;

foreach my $line (@results){
    if ($line =~ m/^$dom_opt/ && $countO <= 1){
        print $line;
        $countO++;
    }
    elsif ($line =~ m/^$domP1/ && $countP1 <= 1){
        print $line;
        $countP1++;
    }
    elsif ($line =~ m/^$domP2/ && $countP2 <= 1){
        print $line;
        $countP2++;
    }
    elsif ($line =~ m/^$domM1/ && $countM1 <= 1){
        print $line;
        $countM1++;
    }
    elsif ($line =~ m/^$domM2/ && $countM2 <= 1){
        print $line;
        $countM2++;
    }
    elsif ($line !~ m/^\d/){
        print $line;
    }
}
