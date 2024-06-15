#!/usr/bin/perl -w

use strict;
use POSIX qw(ceil floor);
#my $limit_size_dom=30;
if ($#ARGV < 1) {print STDERR "Error!\nUsage: $0 <file input matrix contact> <file \nExample: $0 file.mtx file_pu_delineation.mtx\n";exit 1;}

my $file_contact=$ARGV[0];
my $file_pu     =$ARGV[1];

print "#FILE CONTACT: $file_contact\n";
print "#FILE PU     : $file_pu\n";

my $cutoff_pdp=$ARGV[2];
my $max_number_results=500;  #OLD 200
#my $cutoff_protein_ratio_pdp =100000.0; #OLD 20

my $cutoff_nnc=0.001; #1.0

if (not defined $cutoff_pdp)
{
     $cutoff_pdp=0.0001;
}

print "# CUTOFF_PDP (merge) (cutoffratio PDPM) : $cutoff_pdp \n";
print "# MAX NUMBER OF RESULTS:                  $max_number_results\n";
#print "# CUTOFF_PROTEIN_RATIO_PDP:               $cutoff_protein_ratio_pdp\n";

# 1) OPEN AND SAVE DELINEATION OF PU
####################################

my @tab_pu;
my @tab_pu_start;
my @tab_pu_end;
my @tab_pu_size;
open(F,"$file_pu") || die "Cannot open $file_pu :$!\n";
my @tab_file_pu=<F>;
close F;
my %hash_pu;
for (my $i=0 ; $i <= $#tab_file_pu ; $i++)
{
    my $line=$tab_file_pu[$i];
    chomp $line;
    my @tab_line=split(/\s+/,$line);
    my $id_pu=$tab_line[0];
    my $start=$tab_line[1];
    my $end  =$tab_line[2];
    my $size=$end-$start+1;
    $hash_pu{$start}=$line;
}
###########################
#   Sort PUs
###########################
#Sort PU by start
my %hash_correspondance;
my $pu=0;

my %hash_pu_start_end;


my @tab_domains;
my $total_size=0;
my $cutoff_size_domain=30;

foreach my $key (sort { $a <=> $b } keys %hash_pu) 
{
    my $line=$hash_pu{$key};
    my @tab_line=split(/\s+/,$line);
    my $id_pu=$tab_line[0];
    my $start=$tab_line[1];
    my $end  =$tab_line[2];
    my $size=$end-$start+1;

    $total_size=$size+$total_size;

    #For correspondance in contact file
    $hash_correspondance{$id_pu}=$pu;

    $tab_pu[$pu]="$start"." $end"." $size";
    $tab_pu_start[$pu]="$start";
    $tab_pu_end[$pu]  ="$end";
    $tab_pu_size[$pu] ="$size";

    my $true_id_pu=$pu+1;
    $hash_pu_start_end{$true_id_pu}="$tab_pu_start[$pu]-$tab_pu_end[$pu]";
    print "# KEY $key ID_PU:$true_id_pu $tab_pu[$pu] START:$tab_pu_start[$pu] END:$tab_pu_end[$pu]  SIZE:$tab_pu_size[$pu]\n";
    push @tab_domains,$true_id_pu;
   $pu++;
}

# 2) OPEN AND SAVE CONTACTS FILE
################################
my @tab_matrix;
open(F,"$file_contact") || die "Cannot open $file_contact :$!\n";
my @tab_file_contact=<F>;
close F;
my %hash_contact;
for (my $i=0 ; $i <= $#tab_file_contact ; $i++)
{
    my $line=$tab_file_contact[$i];
    chomp $line;
    my @tab_line=split(/\s+/,$line);
    my $x=$tab_line[0];
    my $y=$tab_line[1];
    my $value=$tab_line[2];
    my $a=$hash_correspondance{$x};
    my $b=$hash_correspondance{$y};
    $tab_matrix[$a][$b]=$value;
}

my %hash_sort_pdp;
my $old_number_domain=1;
my @tab_output;

# 3) WHILE CREATE DOM WHITH PU #####
####################################

my   @tab_old_domains;
my   $number_domains=$#tab_domains;
my   $max_number_domains=$#tab_domains+1;
my   @tab_all_tab_domains; # Tab with all domains
push @tab_all_tab_domains,[(@tab_domains)];

my %hash_info_domains;

my $number_domains_1 = $#tab_domains+1;
my $print_domain_1     = print_domain(\@tab_domains,\%hash_pu_start_end);

my ($minimal_measured_size2_1,$maximal_cr_1,$mean_cr_1,$density_min_1,$mean_density_1)=measure_domain([(@tab_domains)]);

printf("#%-1s|%-3s|%50s|%10s|%10s|%10s|%10s|\n",
                "D",
                "min",
                "Domains",
                "MAX_CR",
                "MEAN_CR",
                "DENSIT_MIN",
                "MEAN_DENSITY"
                );

printf("%-2d|%-3d|%50s|%10.6f|%10.6f|%10.6f|%10.6f\n",
                    $number_domains_1,
                    $minimal_measured_size2_1,
                    $print_domain_1,
                    $maximal_cr_1,
                    $mean_cr_1,
                    $density_min_1,
                    $mean_density_1
);


while (1)
{
    #Check if exist same doms in array tab_all_domains and clean it
    #--------------------------------------------------------------
    my @tab_new_all_tab_domains; #New array
    FORSAME:for (my $i=0 ; $i <= $#tab_all_tab_domains ; $i++)
    {
        my @tab_tmp_domain1=(@{$tab_all_tab_domains[$i]});
        my $doms1=join(" ",@tab_tmp_domain1);
        #print "  LINE115: \$dom1:$doms1\n";
        for (my $j=$i+1 ; $j <= $#tab_all_tab_domains ; $j++)
        {
            my @tab_tmp_domain2=(@{$tab_all_tab_domains[$j]});
            my $doms2=join(" ",@tab_tmp_domain2);
            #print "    LINE 120:\$dom2:$doms2\n";
            if ($doms1 =~ /$doms2/) {
                #print " LINE 121: SKIP $doms1\n";
                next FORSAME;}
        } 
        push @tab_new_all_tab_domains,[@tab_tmp_domain1];
    }
    @tab_all_tab_domains=(@tab_new_all_tab_domains);

    # Try differenet merging Compute criterion
    #-----------------------------------------
    my @tab_all_results=();
    for (my $i=0 ; $i <= $#tab_all_tab_domains ; $i++)
    {
        my $ref_tab_domains=$tab_all_tab_domains[$i];

        my $ref_all_results=compute_merge_pu($ref_tab_domains);

        my @tab_output=(@$ref_all_results);
        @tab_all_results=(@tab_all_results,@tab_output);
    }

    # Sort fragment in tab
    #--------------------------------------------------------------
    for (my $i=0 ; $i <= $#tab_all_results ; $i++)
    {
        my @tab_tmp_results=(@{$tab_all_results[$i]});
        my @tab_doms=(@{$tab_tmp_results[0]});
        for (my $j=0 ; $j <= $#tab_doms ; $j++)
        {
            my $dom=$tab_doms[$j];
            my $new_dom=$dom;
            if ($dom=~/\;/)
            {
                my @tab_pu=split(/\;/,$dom);
                my @tab_sorted_pu = sort { $a <=> $b } @tab_pu;
                $new_dom=join(';',@tab_sorted_pu);
            }
            $tab_doms[$j]=$new_dom;
        }
        #my @tab_sorted_doms=sort { $a <=> $b } @tab_doms;
        my @tab_sorted_doms =  map { $_->[0] } sort { $a->[1] <=> $b->[1] } map { [$_, /^(\d{1,})/] } @tab_doms;

        ${$tab_all_results[$i]}[0]=[(@tab_sorted_doms)];
#        print "           SORT FRAGMENTS -> ||||| @{${$tab_all_results[$i]}[0]} \n";
    }


    #Sort results by one criterion
    #-----------------------------
    my %hash_all_results;
	my @tab_sorting_all_results;
    for (my $i=0 ; $i <= $#tab_all_results ; $i++)
    {
#        print " LINE 188 @{$tab_all_results[$i]}\n";
        my ($ref_tab_domains,$new_domain,$ref_old_tab_domains,$minimal_measured_size,$ratio_pdp_merge)=@{$tab_all_results[$i]};
        #my $rand=rand(0.0001);
        #my $value=$s_tot+$rand;
	#print "#JCG RATIO_PDP_MERGE:$ratio_pdp_merge RATIO_PDP_MERGE+RAND: ";
        #TMP
	#my $value=$ratio_pdp_merge+$rand;
        #$hash_all_results{$value}=[@{$tab_all_results[$i]}];
	#TMP
	#VALUE FOR SORTING RATIO PDP MERGE
	push @tab_sorting_all_results,$ratio_pdp_merge;
    } 

	#my @list_order = sort { $stuff_list[$a] cmp $stuff_list[$b] } 0 .. $#stuff_list;
	my @tab_order_sorting_all_results = sort { $tab_sorting_all_results[$b] <=> $tab_sorting_all_results[$a] } 0 .. $#tab_sorting_all_results;


    my @tab_new_all_results=();
    my $limit_results=$max_number_results;
    my $number_results=0;
    
	SORT:for (my $i=0 ; $i <= $#tab_order_sorting_all_results ; $i++)
	{
		my $sorted_indice=$tab_order_sorting_all_results[$i];
		#printf("$i SORTED INDICE $sorted_indice->%s \n", $tab_sorting_all_results[$sorted_indice]);
		push @tab_new_all_results,[@{$tab_all_results[$sorted_indice]}];
	        $number_results++;
	        last SORT if ($number_results > $limit_results);
	}

#SORT:foreach my $key (sort { $b <=> $a  } keys %hash_all_results) 
#    {
#        #last SORT if ($key < 1);
#	print "SORTHASH $number_results $key\n";
#        #push @tab_new_all_results,[@{$hash_all_results{$key}}];
##        print "LINE 226           SORT BY PDP -> |||||| $key @{${$hash_all_results{$key}}[0]} ${$hash_all_results{$key}}[1] @{${$hash_all_results{$key}}[2]} ${$hash_all_results{$key}}[3] ${$hash_all_results{$key}}[4] ${$hash_all_results{$key}}[5]\n";
#        $number_results++;
#        last SORT if ($number_results > $limit_results);
#    }
#
    @tab_all_results=(@tab_new_all_results);
    
    # Print Resutls
    #---------------
printf("#%-1s|%-3s|%50s|%10s|%10s|%10s|%10s|\n",
                "D",
                "min",
                "Domains",
                "MAX_CR",
                "MEAN_CR",
                "DENSIT_MIN",
                "MEAN_DENSITY"
                );


    @tab_all_tab_domains=();

    my %hash_print_results;
    my %hash_domains;

	my @tab_sorting_print_results;
	my @tab_sorting_domains;

	my @tab_print_results;

    for (my $i=0 ; $i <= $#tab_all_results ; $i++)
    {
        my ($ref_tab_domains,$new_domain,$ref_old_tab_domains,$minimal_measured_size,$normalized_external_contact, $ratio_pdp_merge, $mcc_criterion,$s)=@{$tab_all_results[$i]};
        @tab_domains         = ();
        @tab_domains         = (@$ref_tab_domains);

        @tab_old_domains     = ();
        @tab_old_domains     = (@$ref_old_tab_domains);

        my $print_domain     = print_domain(\@tab_domains,\%hash_pu_start_end);
        my $print_old_domain = print_domain(\@tab_old_domains,\%hash_pu_start_end);

        $print_domain        = merge_dom($print_domain);
        $print_old_domain    = merge_dom($print_old_domain);

        # Analyse domains delineation
        #-----------------------------

        my ($minimal_measured_size2,$maximal_cr,$mean_cr,$density_min,$mean_density)=measure_domain([(@tab_domains)]);

        # Check if value to cutoff
        ###########WARNING###############
        #if ($new_ratio_pdp < $cutoff_protein_ratio_pdp or $number_domains==1)

        my $pouet=1;
        if ($pouet==1)
        {
            #print "LINE 290 INSIDE\n";
            my $print_results="";
            #Hash info
            if (not defined $hash_info_domains{$print_domain})
            {
                #print "LINE 295 Define has\n";
                $hash_info_domains{$print_domain}=
                {
                    minimal_measured_size       => "$minimal_measured_size2",
                    maximal_cr                  => "$maximal_cr",
                    mean_cr                     => "$mean_cr",
                    density_min                 => "$density_min",
                    mean_density                => "$mean_density"
                };

                    $print_results=sprintf("%-2d|%-3d|%50s|%10.6f|%10.6f|%10.6f|%10.6f",
                    $number_domains,
                    $minimal_measured_size2,
                    $print_domain,
                    $maximal_cr,
                    $mean_cr,
                    $density_min,
                    $mean_density
                );
# JCG 2011  $hash_print_results{$new_ratio_pdp} =$print_results;
		
		push @tab_sorting_print_results,$mean_density;
		push @tab_sorting_domains,[(@tab_domains)];
		push @tab_print_results,$print_results;

                $hash_print_results{$mean_density}  =$print_results;
# JCG 2011  $hash_domains{$new_ratio_pdp}       =[(@tab_domains)];
                $hash_domains{$mean_density}        =[(@tab_domains)];
            }
            else
            {
                #   print "This domain exists!\n Create a new one ?\n";
            }
        }
    }

# Sort results #
#--------------#

    my $number=0;

     my @tab_order_sorting_print_results = sort { $tab_sorting_print_results[$a] <=> $tab_sorting_print_results[$b] } 0 .. $#tab_sorting_print_results;

	SORTING:for (my $j=0 ; $j <= $#tab_order_sorting_print_results ; $j++)
	{
		my $sorted_indice=$tab_order_sorting_print_results[$j];
	        push @tab_all_tab_domains,[(@{$tab_sorting_domains[$sorted_indice]})];
        #Cutoff limit of maximal number of domain
        my $limit_number_of_domains=ceil($total_size /$cutoff_size_domain);
        next SORTING if ($number_domains > $limit_number_of_domains);
        printf("$tab_print_results[$sorted_indice]\n");
        $number++;
        last SORTING if $number > $max_number_results;
    }

	
#
#
#    SORTING:foreach my $value (sort { $a <=> $b  } keys %hash_print_results)
#    {
#        #print "LINE 361 INSIDE SORTING\n";
#        push @tab_all_tab_domains,[(@{$hash_domains{$value}})];
#        #Cutoff limit of maximal number of domain
#        my $limit_number_of_domains=ceil($total_size /$cutoff_size_domain);
#        next SORTING if ($number_domains > $limit_number_of_domains);
#        printf("$hash_print_results{$value}\n");
#        $number++;
#        last SORTING if $number > $max_number_results;
#    }
#
    $number_domains--;

    last if ($number_domains==0);
    last if ($#tab_domains == 0);

}


###################################
###################################
###                             ###
###       SUB PROGRAMS          ###
###                             ###
###################################
###################################

#################
### MERGE PU ####
#################

sub compute_merge_pu
{
    my $ref_tab_domains=shift;
    my $debug=0;
    my @tab_domains=(@$ref_tab_domains);
    my @tab_old_domains=(@$ref_tab_domains);

    #1 MEASURE INTERNAL 1 and 2 and  EXTERNAL
    #----------------------------------------
    my $check_pdp=1;
    my $total_pdp=0;
    my $check_size=0;

    my %hash_results;


my @tab_results;
my @tab_sorted_results;

    # Foreach PU/DOMAINS
    FORJ:for (my $j=0 ; $j <= $#tab_domains ; $j++)
    {

        my $dom1=$tab_domains[$j];

        #SEARCH IF DOMAINS IS FRANGMENTED
        my @tab_secondary_dom1;
        if ($dom1 =~/;/)
        {
            @tab_secondary_dom1=split(/;/,$dom1);
        }
	else
	{
		push @tab_secondary_dom1,$dom1;
	}

#MEASURE INTERNAL CONTACTS
	my ($size_internal1,$contact_internal1)=measure_internal(\@tab_secondary_dom1,\@tab_matrix,\@tab_pu_size);

## JCG 2014 printf("x:$j size:$size_internal1 P:$contact_internal1 \n");

#Foreach PU/DOMAINS
FORK:for (my $k=$j ; $k <= $#tab_domains ; $k++)
     {
#Skip same PU/DOMAINS
	     next FORK if ($j==$k);

	     my $dom2=$tab_domains[$k];

	     my @tab_secondary_dom2;
	     if ($dom2 =~/;/)
	     {
		     @tab_secondary_dom2=split(/;/,$dom2);
	     }
	     else
	     {
		     push @tab_secondary_dom2,$dom2;
	     }


# COMPUTE INTERNAL SIZE AND INTERNAL CONTACTS 
	     my ($size_internal1,$contact_internal1) =measure_internal(\@tab_secondary_dom1,\@tab_matrix,\@tab_pu_size);
	     my ($size_internal2,$contact_internal2) =measure_internal(\@tab_secondary_dom2,\@tab_matrix,\@tab_pu_size);

## JCG 2014 print(" y:$k size:$size_internal2 P:$contact_internal2\n");

#External
	     my ($size_external,$contact_external)   =measure_external(\@tab_secondary_dom1,\@tab_secondary_dom2,\@tab_matrix,\@tab_pu_size);

## JCG 2014 print("   x-y:$j $k  P:$contact_external\n");

# COMPUTE PDP criterion
###########################
# EXT / SIZE : Hight value many inter(ext) contacts => merge 
#SKIP IF THERE IS NEAR NO CONTACTS ?
	     if ($contact_external <= 0.001) {$contact_external = 0.001;if ($debug ==1) {print STDERR "  Contact external 0 between ID1:$j and ID2:$k\n";};}

	     my $normalized_external_contact     = $contact_external/(($size_internal1 ** 0.43) * ($size_internal2 ** 0.43));

##JCG 2014 printf("      nnc: $normalized_external_contact (contact external:$contact_external)\n");

# JCG 2014 my $modulo                          = rand(0.0001);
#JCG 2014 $normalized_external_contact        = $normalized_external_contact+$modulo;

	     my $ratio_pdp                       = $normalized_external_contact/( ($contact_internal1+$contact_internal2) / ($size_internal1 * $size_internal2) ); # Hight external contacts imply hight ratio_pdp value
# JCG 2014 $ratio_pdp                          = $ratio_pdp+$modulo;

		my $average_number_of_contact1      = ($contact_internal1+$contact_internal2) / ($size_internal1 * $size_internal2);
	     	my $average_number_of_contact2      = ($contact_internal1+$contact_internal2) / ($size_internal1 + $size_internal2);

#NEWNEWNEW PDP CRITERION
	     my $new_new_new_pdp_criterion       = ($contact_external+$contact_internal1+$contact_internal2)/($size_internal1 + $size_internal2);
	     my $normalized_internal_contact1    = 0.5*($contact_internal1/$size_internal1);
	     my $normalized_internal_contact2    = 0.5*($contact_internal2/$size_internal2);

#Normalized contact whole domain
	     my $normalized_contact_whole_domain = ($contact_internal1+$contact_internal2+$contact_external)/($size_internal1 + $size_internal2);
	     my $total_contact_whole_domain      = ($contact_internal1+$contact_internal2+$contact_external);

#######################
	     my $normalized_internal_contact     = 0.5*($contact_internal2/$size_internal2)+($contact_internal1/$size_internal1);
	     my $s                               = $normalized_external_contact/$normalized_internal_contact;
#JCG 2014 $s                                  = $s+$modulo;

#ajouter CR
#et classer selon CR
	     push @tab_results,"$j $k $normalized_external_contact $ratio_pdp";
	     push @tab_sorted_results, $ratio_pdp;

	     $hash_results{"$ratio_pdp"}         = "$j $k $normalized_external_contact $ratio_pdp";
#$hash_results{"$mean_cr"}="$j $k $normalized_external_contact $ratio_pdp";

     }

    }
	 my $check=0;
	 my @tab_all_results;
	 my $number_of_results=0;

#Sort by specific criterion (the higher is the first the reason is that many contact imply probable merge)
#########################################################################################################
	 my @tab_order_sorted_results = sort { $tab_sorted_results[$b] <=> $tab_sorted_results[$a] } 0 .. $#tab_sorted_results;

SORT:for (my $i=0 ; $i <= $#tab_order_sorted_results ; $i++)
     {
	     my $sorted_indice=$tab_order_sorted_results[$i];

#   SORT:foreach my $value (sort { $b <=> $a  } keys %hash_results)
#   {
	my ($id1,$id2,$normalized_external_contact,$ratio_pdp)=split(/ /,$tab_results[$sorted_indice]);
#        my ($id1,$id2,$normalized_external_contact,$ratio_pdp)=split(/ /,$hash_results{$value});
	#$ratio_pdp=$value;
	my $pu1=$tab_domains[$id1];
	my $pu2=$tab_domains[$id2];

#Cutoff criterion and always keep first results
###############################################
	$check=1; #WARNING DON T KEEP THE FIRST RESULTS

####### WARNING #####
#if ($nnc > $cutoff_nnc)
		if ($number_of_results < 1 or $normalized_external_contact > 0.01)
		{
			my $new_pu="$pu1;$pu2";
			my @tab_new_domains=(@tab_domains);
			$tab_new_domains[$id1]="$new_pu";
			splice (@tab_new_domains, $id2, 1);
			my $minimal_measured_size=measure_min_size(\@tab_new_domains,\@tab_matrix,\@tab_pu_size);
			push @tab_all_results,[[@tab_new_domains],$new_pu,[@tab_old_domains],$minimal_measured_size,$normalized_external_contact];
			$number_of_results++;
		}

}
return(\@tab_all_results);
}

# SUB PROGRAM PRINT DOMAIN
##########################

sub print_domain
{

	my $ref_tab_domains=shift;
	my $ref_hash_pu_start_end=shift;

	my @tab_domains=(@$ref_tab_domains);

	my $all_domain_delineation="";
	for (my $a= 0 ; $a <= $#tab_domains ; $a++)
	{
		my @tab_sub_domains=split(/;/,$tab_domains[$a]);
		my $domain_delineation="";

		my $old_pu=$tab_sub_domains[0];

		for (my $b=0 ; $b <= $#tab_sub_domains ; $b++)
		{

			my $current_pu=$tab_sub_domains[$b];
			my $delineation=$$ref_hash_pu_start_end{$current_pu};
			my $next_pu=$old_pu+1;
			if ($current_pu==$next_pu)
			{
				my @tab_domain_delineation=split('-',$domain_delineation);

				my @tab_current_pu=split('-',$$ref_hash_pu_start_end{$current_pu});
				my $last_num=$tab_current_pu[1];
				$tab_domain_delineation[$#tab_domain_delineation]="$last_num;";
				$domain_delineation=join('-',@tab_domain_delineation);

			}
			else
			{
				$domain_delineation.="$delineation;";
			}
			$old_pu=$current_pu;
		}
		chop $domain_delineation;
		$all_domain_delineation.="$domain_delineation ";
	}
	chop $all_domain_delineation;


#my $print_domain=sprintf("$number_domain|%10s|$minimal_measured_size|$all_domain_delineation\n",$total_pdp);
	my $print_domain=sprintf("$all_domain_delineation");

	return($print_domain);
}



##################################
# SUB PROGRAM MEASURE_INTERNAL
#################################
sub measure_internal
{
	my $ref_tab_secondary_dom=shift;
	my $ref_tab_matrix       =shift;
	my $ref_tab_pu_size      =shift;

#Internal 
	my $contact_internal=0; 
	my $size_internal=0;
	for (my $l=0 ; $l <= $#$ref_tab_secondary_dom ;$l++)
	{
		my $id_pu_dom_1   = $$ref_tab_secondary_dom[$l]-1;
#my $startm=$l;
		for (my $m=0 ; $m <= $#$ref_tab_secondary_dom ;$m++)
		{
			my $id_pu_dom_2=$$ref_tab_secondary_dom[$m]-1;
			$contact_internal+=$$ref_tab_matrix[$id_pu_dom_1][$id_pu_dom_2];
			if ($id_pu_dom_1 =~/^$id_pu_dom_2$/)
			{
				$size_internal   +=$$ref_tab_pu_size[$id_pu_dom_1];
			}
		}

	}
#Divide by 2 because count twice
#$size_internal=$size_internal/2;
	$contact_internal=$contact_internal/2;
	return($size_internal,$contact_internal);
}


##########################
# SUB PROGRAM Measure external
##########################
sub measure_external
{
	my $ref_tab_secondary_dom1=shift;
	my $ref_tab_secondary_dom2=shift;
	my $ref_tab_matrix        =shift;
	my $ref_tab_pu_size       =shift;

	my $contact_external=0;
	my $size_external=0;

	for (my $l=0 ; $l <= $#$ref_tab_secondary_dom1 ;$l++)
	{
		my $id_pu_dom1   = $$ref_tab_secondary_dom1[$l]-1;
#my $startm=$l;
		$size_external   +=($$ref_tab_pu_size[$id_pu_dom1]);

		for (my $m=0 ; $m <= $#$ref_tab_secondary_dom2 ;$m++)
		{
			my $id_pu_dom2   = $$ref_tab_secondary_dom2[$m]-1;
			$contact_external+= $$ref_tab_matrix[$id_pu_dom1][$id_pu_dom2];
#MODIF 20100701
			$size_external   +=($$ref_tab_pu_size[$id_pu_dom2]);
#ENDMODIF 20100701
		}

	}
	$size_external=$size_external/2;
	$contact_external=$contact_external;
	return($size_external,$contact_external);
}


##########################
# Measure Domains partition quality
##########################
sub measure_domain
{
	my $ref_tab_domains=shift;
	my $debug=0;
	my @tab_domains=(@$ref_tab_domains);

	if ($#tab_domains==0)
	{
		my $size = &measure_min_size(\@tab_domains,\@tab_matrix,\@tab_pu_size);
		return($size,"0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0");
	}


#1 MEASURE INTERNAL 1 and 2 and  EXTERNAL
#---------------------------------------#

	my $minimal_measured_size=100000;


	my $total_contact_internal=0;
	my $total_contact_external=0;

	my $total_size_internal=0;
	my $total_size_external=0;

	my $maximal_nnc=0;

	my $total_pdp_external=0;


	my $number_of_domains=$#tab_domains+1;

	my $total_nnc=0;

	my $total_cr=0;
	my $maximal_cr=0;

	my $nnc_tot=0;

	my $min_density=10000000000;
	my $density_tot=0;

	my $external_number=0;


FORJ:for (my $j=0 ; $j <= $#tab_domains ; $j++)
     {

	     my $dom1=$tab_domains[$j];
	     my @tab_secondary_dom1;

	     if ($dom1 =~/;/)
	     {
		     @tab_secondary_dom1=split(/;/,$dom1);
	     }
	     else
	     {
		     push @tab_secondary_dom1,$dom1;
	     }

	     my ($size_internal1,$contact_internal1)=measure_internal(\@tab_secondary_dom1,\@tab_matrix,\@tab_pu_size);

	     $total_contact_internal=$total_contact_internal+$contact_internal1;
	     $total_size_internal   =$total_size_internal   +$size_internal1;

#Save the actual minimal size measured
	     if ($minimal_measured_size > $size_internal1)
	     {
		     $minimal_measured_size=$size_internal1;
	     }
	     my $density1=$contact_internal1/$size_internal1;

	     $density_tot=$density_tot+$density1;

FORK:for (my $k=$j ; $k <= $#tab_domains ; $k++)
     {

	     next FORK if ($j==$k);
#Normalization du score de contact

	     my $dom2=$tab_domains[$k];

	     my @tab_secondary_dom2;
	     if ($dom2 =~/;/)
	     {
		     @tab_secondary_dom2=split(/;/,$dom2);
	     }
	     else
	     {
		     push @tab_secondary_dom2,$dom2;
	     }


# COMPUTE INTERNAL SIZE AND INTERNAL CONTACTS 
	     my ($size_internal2,$contact_internal2)=measure_internal(\@tab_secondary_dom2,\@tab_matrix,\@tab_pu_size);
	     $total_contact_internal=$total_contact_internal+$contact_internal2;

#Current TOTAL INTERNAL
	     my $current_contact_internal=($contact_internal1+$contact_internal2)/($size_internal1+$size_internal2);;

#External
	     my ($size_external,$contact_external)=measure_external(\@tab_secondary_dom1,\@tab_secondary_dom2,\@tab_matrix,\@tab_pu_size);
	     $total_contact_external=$total_contact_external+$contact_external;
	     $total_size_external   =$total_size_external+$size_external;


#current PDP criterion
	     my $current_nnc=$contact_external/( ($size_internal1 ** 0.43) * ($size_internal2 ** 0.43) );


	     $nnc_tot=$current_nnc+$nnc_tot;

###################################

#NEW PDP CRITERION
	     my $new_new_new_pdp_criterion=($contact_external+$contact_internal1+$contact_internal2)/($size_internal1 + $size_internal2);



#Current Contact ratio
	     my $current_cr=$current_nnc/$new_new_new_pdp_criterion;


#  my $ratio_pdp                 = $normalized_external_contact/( ($contact_internal1+$contact_internal2) / ($size_internal1 * $size_internal2) ); # Hight external contacts imply hight ratio_pdp value

#~ print "708:  current_cr : $current_cr\n";

#Checkand save maximal CR
	     if ($maximal_cr <= $current_cr)
	     {
		     $maximal_cr = $current_cr;
	     }

#NEWSNEWS NEWS NEWS #New S
	     my $normalized_internal_contact=0.5*(($contact_internal2/$size_internal2)+($contact_internal1/$size_internal1));

#####################################
#DENSITY of each domain
	     my $density2=$contact_internal2/$size_internal2;

#Check Save minimal density
###########################
	     if ($min_density > $density1)
	     {
		     $min_density=$density1;
	     }
	     if ($min_density > $density2)
	     {
		     $min_density=$density2;
	     }

#CHECK SIZE
###########
	     if ($minimal_measured_size  > $size_internal2)
	     {
		     $minimal_measured_size  = $size_internal2;
	     }

#TOTAL CR (to compute mean CR)
##############################
	     $total_cr=$current_cr+$total_cr;

     }

     }
#MEAN DENSITY CRITERION (MEAN CONTACT PROBABILITY DENSITY)
##########################################################
     my $mean_density=$density_tot/$number_of_domains;

####################################################

#MEAN CR (CONTACT RATIO)
     my $mean_cr=$total_cr/$number_of_domains;
     $mean_cr=0;
     return($minimal_measured_size,$maximal_cr,$mean_cr,$min_density,$mean_density);
}

################
# merge PU in domain
#################
sub merge_dom {
	my $dom_peel=shift;
	my $new_dom_peel=$dom_peel;
	my @tab_dom_peel=split(/\s+/,$dom_peel);
	for (my $i=0 ; $i <= $#tab_dom_peel ; $i++)
	{
		my $dom=$tab_dom_peel[$i];
		if ($dom =~/;/)
		{
			my @tab_dom_pu=split(/;/,$dom);
FOR1:for (my $a=0 ; $a <= $#tab_dom_pu ; $a++)
     {
	     next FOR1 if ($tab_dom_pu[$a] eq "");
	     my ($startpu1,$endpu1)=split(/\-/,$tab_dom_pu[$a]);
	     my $startb=$a+1;
FOR2:for (my $b=$startb ; $b <= $#tab_dom_pu ; $b++)
     {
	     next FOR2 if ($tab_dom_pu[$b] eq "");
	     my ($startpu2,$endpu2)=split(/\-/,$tab_dom_pu[$b]);
	     my $jonction_1_2=$startpu2-1;
	     my $jonction_2_1=$startpu1-1;
	     if ($endpu1 == $jonction_1_2)
	     {
		     my $new_pu="$startpu1"."-"."$endpu2";
		     $tab_dom_pu[$a]=$new_pu;
		     $tab_dom_pu[$b]="";

		     my @tab_new_dom_pu=(@tab_dom_pu);
		     my $new_dom=join(";",@tab_new_dom_pu);
		     $new_dom =~s/\;{1,}$//g;
		     $new_dom =~s/\;{2,}/\;/g;
		     $tab_dom_peel[$i]=$new_dom;
		     $new_dom_peel=join(" ",@tab_dom_peel);
#print " LINE 622:MERGE 1\n";
		     $new_dom_peel=merge_dom($new_dom_peel); 
	     }
	     elsif ($endpu2 == $jonction_2_1)
	     {
		     my $new_pu="$startpu2"."-"."$endpu1";
		     $tab_dom_pu[$b]=$new_pu;
		     $tab_dom_pu[$a]="";

		     my @tab_new_dom_pu=(@tab_dom_pu);
		     my $new_dom=join(";",@tab_new_dom_pu);
		     $new_dom =~s/\;{1,}$//g;
		     $new_dom =~s/\;{2,}/\;/g;
		     $tab_dom_peel[$i]=$new_dom;
		     $new_dom_peel=join(" ",@tab_dom_peel);
		     $new_dom_peel=merge_dom($new_dom_peel); 

	     }
     }
     }

		}
	}
	return($new_dom_peel);
}

##########################
# Compute smallest domain
##########################
sub measure_min_size
{
	my $ref_tab_domains=shift;
	my @tab_domains=(@$ref_tab_domains);
	my $ref_tab_matrix       =shift;
	my $ref_tab_pu_size      =shift;

	my $minimal_measured_size=100000;

FORJ:for (my $j=0 ; $j <= $#tab_domains ; $j++)
     {

	     my $dom1=$tab_domains[$j];
	     my @tab_secondary_dom1;

	     if ($dom1 =~/;/)
	     {
		     @tab_secondary_dom1=split(/;/,$dom1);
	     }
	     else
	     {
		     push @tab_secondary_dom1,$dom1;
	     }
	     my ($size_internal1,$contact_internal1)=measure_internal(\@tab_secondary_dom1,$ref_tab_matrix,$ref_tab_pu_size);

	     if ($minimal_measured_size > $size_internal1)
	     {   
		     $minimal_measured_size=$size_internal1;
	     }


     }
     return($minimal_measured_size);
}
