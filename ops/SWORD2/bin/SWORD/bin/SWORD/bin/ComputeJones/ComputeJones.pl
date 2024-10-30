#!/usr/bin/env perl

###########################################################
#CHEBREK ROMAIN
#14/04/2011
#Programme de Calcul du pourcentage de d'overlap (score Jones)

# GP 2014: originally named Compute_Jones_alter.pl
###########################################################

use strict;
use File::Basename;
use lib dirname (__FILE__);

my $delineation1	= 	shift;
my $delineation2	= 	shift;
my $dir_data		=	shift;

$delineation1 	= join(" ",split("_",$delineation1));
$delineation2	= join(" ",split("_",$delineation2));
chomp $delineation1;
chomp $delineation2;

my $pdb 		= 	"";
my $dom_peel 		= 	"";
my $dom_auth		= 	"";
my %hash_domain;

my @tab_del1	=	split(" ",$delineation1);
my @tab_del2	=	split(" ",$delineation2);

$pdb = @tab_del1[0];

if($#tab_del1 < $#tab_del2)
{
	$dom_peel 	= join(" ",@tab_del2[2..$#tab_del2]);
	$dom_auth 	= join(" ",@tab_del1[2..$#tab_del1]);
}else{
	$dom_peel 	= join(" ",@tab_del1[2..$#tab_del1]);	
	$dom_auth 	= join(" ",@tab_del2[2..$#tab_del2]);
}


my $val = 0;
my @temp = split(" ",$dom_peel);
if($temp[0] =~ /^-/)
{
	for (my $i = 0; $i <= $#temp; $i++)
	{		
		my @temp2 = split("-",$temp[$i]);
		if($temp[$i] =~ /^-/)
		{
			$val = $temp2[1];
			$temp2[1] = $val*(-1);
			for (my $j = 1; $j <= $#temp2 ;$j++)
			{
				my @temp3 = split(";",$temp2[$j]);
				if($#temp3 > 0)
				{
					for (my $k = 0; $k <= $#temp3 ;$k++)
					{
						$temp3[$k] = $temp3[$k]+$val;
					}
					$temp2[$j] = join(";",@temp3);
				}else{
					$temp2[$j] = $temp2[$j]+$val;
				}
			}
		$temp[$i] = join("-",@temp2[1..$#temp2]);
		}else{
			for (my $j = 0; $j <= $#temp2 ;$j++)
			{
				my @temp3 = split(";",$temp2[$j]);	
				if($#temp3 > 0)
				{
					for (my $k = 0; $k <= $#temp3 ;$k++)
					{
						$temp3[$k] = $temp3[$k]+$val;
					}
					$temp2[$j] = join(";",@temp3);
				}else{
					$temp2[$j] = $temp2[$j]+$val;
				}
			}
			$temp[$i] = join("-",@temp2[0..$#temp2]);
		}
	}
}
$dom_peel = join(" ",@temp);

my $val = 0;
my @temp = split(" ",$dom_auth);
if($temp[0] =~ /^-/)
{
	for (my $i = 0; $i <= $#temp; $i++)
	{		
		my @temp2 = split("-",$temp[$i]);
		if($temp[$i] =~ /^-/)
		{
			$val = $temp2[1];
			$temp2[1] = $val*(-1);
			for (my $j = 1; $j <= $#temp2 ;$j++)
			{
				my @temp3 = split(";",$temp2[$j]);
				if($#temp3 > 0)
				{
					for (my $k = 0; $k <= $#temp3 ;$k++)
					{
						$temp3[$k] = $temp3[$k]+$val;
					}
					$temp2[$j] = join(";",@temp3);
				}else{
					$temp2[$j] = $temp2[$j]+$val;
				}
			}
		$temp[$i] = join("-",@temp2[1..$#temp2]);
		}else{
			for (my $j = 0; $j <= $#temp2 ;$j++)
			{
				my @temp3 = split(";",$temp2[$j]);	
				if($#temp3 > 0)
				{
					for (my $k = 0; $k <= $#temp3 ;$k++)
					{
						$temp3[$k] = $temp3[$k]+$val;
					}
					$temp2[$j] = join(";",@temp3);
				}else{
					$temp2[$j] = $temp2[$j]+$val;
				}
			}
			$temp[$i] = join("-",@temp2[0..$#temp2]);
		}
	}
}
$dom_auth = join(" ",@temp);

#~ print "$dom_peel\n";
#~ print "$dom_auth\n";

${$hash_domain{"$pdb"}}{"cut"}=$dom_auth;

my $domain_authors	=	${$hash_domain{"$pdb"}}{"cut"};
$domain_authors		=~	s/; /;/g;

my @tab_authors=split(/\s+/,$domain_authors);
my @tab_sorted_domain_raw  = map { $_->[0] } 
sort { $a->[1] <=> $b->[1] || $a->[2] cmp $b->[2] } map { [$_, /^(\d+)/, uc($_)] } @tab_authors;

my $tmp_dom=join("-",@tab_sorted_domain_raw);

$tmp_dom=~s/;/;-/g;
my @tab_dom_authors=split(/-/,$tmp_dom);
my $start_authors=6666;
my $last_authors=-6666;
my $flag=0;

for (my $tmp=0 ; $tmp <= $#tab_dom_authors ; $tmp++)
{
	my $dom=$tab_dom_authors[$tmp];
	if ($dom =~/;/) {chop $dom;$flag=1;}
	if ($start_authors>$dom) {$start_authors=$dom;}
	if ($last_authors<$dom) {$last_authors=$dom;}
	if ($flag==1){$flag=0;$dom="$dom".";";};
}
my $length_authors=$last_authors-$start_authors;
my $domain_sorted_raw=join(" ",@tab_sorted_domain_raw);

my $domain_sorted_raw=join(" ",@tab_sorted_domain_raw);


###########################################################
###########################################################

my @tab_protein_domain_authors;
my $num_dom		=1;

my @buff_pdb	=split('',$pdb);
my $inputchain	=$buff_pdb[$#buff_pdb];

#my $temp 		= `perl $prog_parsing $pdb $dir_data../PDBs_Stand/`;

sub num_aa
{
	my $file_pdb		=shift;
    my $selectedchain	=shift;
    my $rep_tmp			="/tmp";

    #~ $selectedchain=uc($selectedchain);

    #~ my $check_chain=0;
    my $check_chain=1;
    #~ if ($selectedchain !~ /^$/) {$check_chain=1;}
        
	open (FILE_PDB,"$file_pdb") || die "Cannot open file:$file_pdb :$!\n";
	my @tab_out=<FILE_PDB>;
	close FILE_PDB;
	my $file_tmp=int(rand (1000));
	my @out;
    my @tab_num;
    FOR:foreach my $line (@tab_out)
	{
		my $record_name = substr($line,  0, 6);  # columns 1-6 = ATOM
		my $number      = substr($line,  6, 5);  # columns 7-11
		my $atom_name   = substr($line, 12, 3);  # columns 13-16
		my $alt_loc     = substr($line, 16, 1);  # columns 17
		my $aa          = substr($line, 17, 3);  # columns 18-20
		my $chain       = substr($line, 21, 1);  # columns 22
		my $no_aa       = substr($line, 22, 4);  # columns 23-26
		my $insert_code = substr($line, 26, 4);  # columns 27
		my $x           = substr($line, 30, 8);  # columns 31-38
		my $y           = substr($line, 38, 8);  # columns 39-46
		my $z           = substr($line, 46, 8);  # columns 47-54
		my $element     = substr($line, 76, 2);  # columns 77-78

        chomp $line;
		#if ($record_name =~/HEADER/) {print OUT "$line";next FOR;}
		# $number and $element may have leading spaces: strip them
		$record_name =~ s/^\s*//;
		$number      =~ s/^\s*//;
		$atom_name   =~ s/^\s*//;
		$chain       =~ s/^\s*//;
		$no_aa       =~ s/^\s*//;
		$element     =~ s/^\s*//;

        # Select the good chain
		if ($selectedchain !~ /^$/) { next FOR if  ($chain !~/$selectedchain/);}
		if ($selectedchain !~ /^$/ and  $chain =~/$selectedchain/) {$check_chain=1;}

        last FOR if ($record_name =~ /TER/   and $check_chain==1);
	    last FOR if ($record_name =~ /ENDMDL/);

        next FOR if ($record_name !~ /^ATOM/);
        #next FOR if ($atom_name   =~ /H/);
        next FOR if ($atom_name   !~ /CA/);

        if    ($aa =~/ALA/i) {}
		elsif ($aa =~/CYS/i) {}
		elsif ($aa =~/ASP/i) {}
		elsif ($aa =~/GLU/i) {}
		elsif ($aa =~/PHE/i) {}
		elsif ($aa =~/GLY/i) {}
		elsif ($aa =~/HIS/i) {}
		elsif ($aa =~/ILE/i) {}
		elsif ($aa =~/LYS/i) {}
		elsif ($aa =~/LEU/i) {}
		elsif ($aa =~/MET/i) {}
		elsif ($aa =~/ASN/i) {}
		elsif ($aa =~/PRO/i) {}
		elsif ($aa =~/GLN/i) {}
		elsif ($aa =~/ARG/i) {}
		elsif ($aa =~/SER/i) {}
		elsif ($aa =~/THR/i) {}
		elsif ($aa =~/VAL/i) {}
		elsif ($aa =~/TRP/i) {}
		elsif ($aa =~/TYR/i) {}
		elsif ($aa =~/UNK/i) {}

		else  {next FOR;}

        push @tab_num,$no_aa;
        push @out,$line;
		#print "$line";
	}
	return (\@tab_num);
}

#my $ref_tab_num_aa_authors=num_aa("$dir_data../PDBs_tmp/$pdb/$pdb.pdb","$inputchain");
my $ref_tab_num_aa_authors=num_aa("$dir_data/$pdb/$pdb.pdb","$inputchain");
#my $temp = `rm -r -f $dir_data../PDBs_tmp/$pdb/`;

my @tab_num_aa_authors=(@$ref_tab_num_aa_authors);

$length_authors=$#tab_num_aa_authors;

for (my $pos=0 ; $pos <= $length_authors ; $pos++)
{
	$tab_protein_domain_authors[$pos]=-1;
}

for (my $pos=0 ; $pos <= $length_authors ; $pos++)
{
	for (my $num_dom=0 ; $num_dom <= $#tab_sorted_domain_raw; $num_dom++)
	{
		my @tab_pu_domain=split(/;/,$tab_sorted_domain_raw[$num_dom]);

		for (my $num_pu=0; $num_pu <= $#tab_pu_domain ; $num_pu++)
		{
			my ($start_pu,$end_pu)=split(/-/,$tab_pu_domain[$num_pu]);
			if ($tab_num_aa_authors[$pos] >= $start_pu and $tab_num_aa_authors[$pos] <= $end_pu)
			{
				$tab_protein_domain_authors[$pos]=$num_dom;
			}
		}
	}
}

###########################################################
###########################################################

sub measure_jones {

    my $dom_peel=shift;
    my $domain_sorted_raw=shift;
    my $ref_tab_protein_domain_authors=shift;
    my $start_authors=shift;

    my $length_authors=$#$ref_tab_protein_domain_authors;

    my @tab_protein_domain_authors=(@$ref_tab_protein_domain_authors);

    my $debug=0;

    #If good cutoff and good_size make some measurments
    ###################################################
    my $start_peel=6666;
    my $last_peel=-6666;
    my $flag_peel=0;
    my $false_dom_peel;
    if ($dom_peel =~/NA/)
    {
        $false_dom_peel="0-$length_authors";
        $dom_peel=$false_dom_peel;
    }

    #Compute length of Peel
    #######################
    $dom_peel=~s/; /;/g;
    $dom_peel=~s/^\s{1,}//g;
    my @tab_peel=split(/\s+/,$dom_peel);
    my $tmp_dom_peel=join("-",@tab_peel);
    $tmp_dom_peel=~s/;/;-/g;
    $tmp_dom_peel=~s/\s+//g;
    my @tab_dom_peel=split(/-/,$tmp_dom_peel);
    if ($debug ==1){    printf STDERR (" LINE 270 CHECK PEEL: \@tab_dom_peel=|@tab_dom_peel|\n");foreach my $coco (@tab_dom_peel){print STDERR "    LINE 270|$coco|\n";};}
    my $flag=0;
    for (my $tmp=0 ; $tmp <= $#tab_dom_peel ; $tmp++)
    {
        my $dom=$tab_dom_peel[$tmp];
        if ($dom =~/;/) {chop $dom;$flag=1;}
        if ($debug == 1) {printf STDERR "  LINE 282: \$start_peel>\$dom:|$start_peel|>|$dom|\$last_peel<\$dom:|$last_peel|<|$dom|\n";} 
        if ($start_peel > $dom) {$start_peel=$dom;}
        if ($last_peel<$dom) {$last_peel=$dom;}
        if ($flag_peel==1){$flag_peel=0;$dom="$dom".";"}
    }
    my $length_peel=$last_peel-$start_peel;
    #my $domain_sorted_peel=join(" ",@tab_doma_peel);

    if ($debug ==1){    printf STDERR ("  LINE 278 \$length_peel=$length_peel \$last_peel=$last_peel \$start_peel=$start_peel\n");};


    #Compute difference between AUTHORS and Peel 
    ############################################
    #Lenght and starting residu
    my $diff=abs($length_peel-$length_authors);
    my $decalage=$start_authors-$start_peel;

    #Change number of residu of peel to be identical to number of residu of author
    ##############################################################################
    my $domain_sorted_peel=join(" ",@tab_dom_peel);

    # Transfert domain appartenance of PEEL at each position in a tab (eg @tab=(11111111111122222211111111111113333...))
    ################################################
    my @tab_protein_domain_peel;
    my $num_dom=1;
    FOR_POS:for (my $pos=0 ; $pos <= $length_peel ; $pos++)
    {
        for (my $num_dom=0 ; $num_dom <= $#tab_peel ; $num_dom++)
        {
            my @tab_pu_domain=split(/;/,$tab_peel[$num_dom]);
            if($debug == 1) {printf STDERR " LINE 315 \$tab_peel[$num_dom]:$tab_peel[$num_dom] \@tab_pu_domain:@tab_pu_domain \n";};

            for (my $num_pu=0; $num_pu <= $#tab_pu_domain ; $num_pu++)
            {
                my ($start_pu,$end_pu)=split(/-/,$tab_pu_domain[$num_pu]);
                if($debug == 1) {printf STDERR " LINE 325 start_pu:$start_pu \$end_pu:$end_pu \$tab_pu_domain[\$num_pu]:$tab_pu_domain[$num_pu] (\@tab_pu_domain:@tab_pu_domain)\n";};

                $start_pu=$start_pu-$start_peel;
                $end_pu=$end_pu-$start_peel;
                if($debug == 1) {printf STDERR " LINE 329: \$pos:$pos \$start_pu=$start_pu \$start_peel:$start_peel  \$end_pu=$end_pu \$start_peel:$start_peel \n";};
                if ($debug == 1) {printf STDERR " LINE 330: \$pos >= \$start_pu and \$pos <= \$end_pu: \$pos >= $start_pu and $pos <= $end_pu : ?\n";}
                if ($pos >= $start_pu and $pos <= $end_pu)
                {
                    $tab_protein_domain_peel[$pos]=$num_dom;
                    next FOR_POS;
                }
            }
        }
    }
    if($debug == 1) {printf STDERR " LINE 326 \@tab_protein_domain_peel: @tab_protein_domain_peel \n";};

#    if ($#tab_protein_domain_peel != $#tab_num_aa_authors)
#    {
#        my $diff_size=($#tab_protein_domain_peel)-($#tab_num_aa_authors);
#        print STDERR " WARNING for $pdb not same size for authors and peel: $diff_size ($#tab_protein_domain_peel)vs($#tab_num_aa_authors)";
#    }

    #Check if number of domains in peel and in authors don't differ too 
    ###################################################################
    my @tab_sorted_authors=split(/\s+/,$domain_sorted_raw);
    my $number_domains_authors=$#tab_sorted_authors+1;
    my $number_domains_peel=$#tab_peel+1; 
    if ($debug == 1) {print STDERR "  LINE 428: \$number_domains_authors=$number_domains_authors \$number_domains_peel=$number_domains_peel\n";}
    if ($number_domains_peel > 7 and ($number_domains_authors < $number_domains_peel))
    {
        return("0","0",$length_authors,$length_peel,$dom_peel,$domain_sorted_raw,$diff);
    }


    #TEST every jones model  
    #######################
    #Permutation
    #~ use Math::Combinatorics;
	require Combinatorics;
    my @tab_domain_perm;

    for (my $k=0 ; $k <=$#tab_peel ; $k++) {push @tab_domain_perm,$k;}
    #printf STDERR " LINE 353: $pdb permutation: k:$#tab_peel \@tab_dom_peel:@tab_peel \n";
    #my @tab_permutation = permute(@tab_domain_perm);
    $number_domains_peel=$#tab_domain_perm+1;


    my $combinat = Math::Combinatorics->new(count => $number_domains_peel, data => [@tab_domain_perm], );

	my $seuil_jones=85;
    my $jones_criterion=0;
    my $jones=0;
    my $old_jones=0;
    my $total=0;

    WHILE_JONES:while(my @tab_ref_tab_current = $combinat->next_permutation)
    {   
        $jones=0;
        $total=0;
        my $ref_tab_current=\@tab_ref_tab_current;
        #printf STDERR " LINE 443  \$ref_tab_current: @$ref_tab_current\n";
        FOR_POS_COMPARE:for (my $pos=0; $pos <= $#tab_protein_domain_peel ; $pos++)
        {
            last FOR_POS_COMPARE if ($pos > $#tab_protein_domain_authors);
            #Translate current domain name by anothers 
            my $actual_protein_peel_dom=$$ref_tab_current[$tab_protein_domain_peel[$pos]];
            #if($debug==1){print STDERR "  LINE 392 \$actual_protein_peel_dom:$actual_protein_peel_dom \$tab_protein_domain_authors[\$pos]:$tab_protein_domain_authors[$pos] (pos:$pos)\n";}

            if ($actual_protein_peel_dom == $tab_protein_domain_authors[$pos])
            {
                $jones++;
            }
            #Check if position in AUthors is determined
            if ($tab_protein_domain_authors[$pos] == -1) 
            {
                if ($debug==1) {print STDERR "WARNING! Undefined domain in AUTHORS for position $pos";};
            }
            else
            {
                $total++;
            }
        }
        if ($old_jones < $jones) {$old_jones=$jones;}
        if ($old_jones > $seuil_jones) {last WHILE_JONES;}
    }
    $jones=$old_jones;

    if ($total == 0) {$total=-1;}
    my $percentage_jones=$jones/$total*100; 

    if ($percentage_jones >= 85)
    {
        $jones_criterion=1;
    }
    #if ($debug==1) {print STDERR " LINE 354: \$smallest_pdp:$smallest_pdp \$smallest_size:$smallest_size";}

    return($jones_criterion,$percentage_jones,$length_authors,$length_peel,$dom_peel,$domain_sorted_raw,$diff);
}

  my($jones_criterion,$fidelity,$length_authors,$length_peel,$dom_peel,$dom_authors,$diff) = measure_jones($dom_peel,$domain_sorted_raw,\@tab_protein_domain_authors,$start_authors);

print "$jones_criterion $fidelity\n";
