
use strict;
use FindBin;
use lib "$FindBin::Bin/bin/ComputeJones";
use DistanceModel;

sub ParseMeasure {

    my ( $PU_measure, $dir_data, $pdb, $option_alt_dist, $Ndom, $alt_b, $alt_l,
        $option_alt_diff )
      = splice( @_, 0, 8 );

    $pdb = $pdb . '_1_';

    my ( @clean_measure, @temp_measure, @relevant_measure, @measure );
    my ( $flag,         $max_dom )      = (0) x 2;
    my ( $delineation1, $delineation2 ) = ('') x 2;

    my $max_dist = 0.2;

    my @PU_measure = split( "\n", $PU_measure );

    if ( $option_alt_diff == 1 ) {
        foreach my $id_line ( 0 .. $#PU_measure ) {
            if ( substr( $PU_measure[$id_line], 0, 1 ) ne '#' ) {
                if ( $max_dom == 0 ) {
                    if (
                        substr( $PU_measure[$id_line], 0, 2 ) >
                        ( $Ndom + $alt_l )
                        && substr( $PU_measure[$id_line], 0, 2 ) > 6 )
                    {
                        push( @clean_measure, $PU_measure[$id_line] );
                    }
                    else {
                        $max_dom = substr( $PU_measure[$id_line], 0, 2 );
                    }
                }

                if ( substr( $PU_measure[$id_line], 0, 2 ) == $max_dom ) {

                    if ( $option_alt_dist == 1 ) {
                        my @temp = split( '\|', $PU_measure[$id_line] );
                        my $abs = 1;
                        if ( DistanceModel( $temp[3], $temp[5], $abs ) <
                            $max_dist )
                        {
                            push( @temp_measure, $PU_measure[$id_line] );
                        }
                    }
                    else {
                        push( @temp_measure, $PU_measure[$id_line] );
                    }
                }

                if ( substr( $PU_measure[$id_line], 0, 2 ) < $max_dom ) {
                    $max_dom = substr( $PU_measure[$id_line], 0, 2 );
                    @temp_measure = reverse(@temp_measure);

                    my @last_temp;

                    if ( $#temp_measure > 0 and $alt_b > 1 ) {
                        my @temp_measure2 = ();
                        foreach my $i ( 0 .. $alt_b - 1 ) {
                            @temp_measure2 =
                              ( $temp_measure[$i], @temp_measure2 );

                            if ( $i == $#temp_measure or $i == $alt_b - 1 ) {
                                last;
                            }
                            else {
                                my $cpt = 0;
                                for (
                                    my $j = $i + 1 ;
                                    $j <= $#temp_measure ;
                                    $j++
                                  )
                                {
                                    $cpt++;

                                    last if ( $cpt > 20 );

                                    if ( $i != $#temp_measure ) {
                                        my @temp =
                                          split( '\|', $temp_measure[$i] );
                                        $delineation1 = $temp[2];
                                        $delineation1 = join( '_',
                                            split( ' ', $delineation1 ) );
                                        $delineation1 =~ s/;/\\;/g;

                                        @last_temp =
                                          split( '\|', $temp_measure[$j] );
                                        $delineation2 = $last_temp[2];
                                        $delineation2 = join( '_',
                                            split( ' ', $delineation2 ) );
                                        $delineation2 =~ s/;/\\;/g;

                                        my $temp =
`$FindBin::Bin/bin/ComputeJones/ComputeJones.pl $pdb$delineation1 $pdb$delineation2 $dir_data`;

                                        if ( substr( $temp, 0, 1 ) == 1 ) {
                                            @temp_measure = (
                                                @temp_measure[ 0 .. $j - 1 ],
                                                @temp_measure[ $j +
                                                  1 .. $#temp_measure ]
                                            );
                                            $j--;
                                        }
                                        else {
                                            last if ( $i == $alt_b - 2 );
                                        }
                                    }
                                }
                            }
                        }

                        @clean_measure = ( @clean_measure, @temp_measure2 );

                    }
                    else {
                        push( @clean_measure, $temp_measure[0] )
                          if ( $#temp_measure != -1 );
                    }

                    @temp_measure = ();
                    if ( DistanceModel( $last_temp[3], $last_temp[5], 1 ) <
                        $max_dist )
                    {
                        push( @temp_measure, $PU_measure[$id_line] );
                    }
                }
            }
        }

        @clean_measure = ( @clean_measure, $PU_measure[$#PU_measure] );

        return @clean_measure if ( $alt_b == 1 );

        $max_dom = 0;

        foreach my $id_line ( 0 .. $#clean_measure ) {
            $max_dom = substr( $clean_measure[$id_line], 0, 2 )
              if ( $max_dom == 0 );

            if ( substr( $clean_measure[$id_line], 0, 2 ) < $max_dom ) {
                for ( my $i = 1 ; $i <= $alt_b ; $i++ ) {
                    if (
                        substr( $clean_measure[ ( $id_line - $i ) ], 0, 2 ) ==
                        $max_dom )
                    {
                        push( @relevant_measure,
                            $clean_measure[ ( $id_line - $i ) ] )
                          if ( ( $id_line - $i ) >= 0 );
                    }
                }
                $max_dom = substr( $clean_measure[$id_line], 0, 2 );
            }
        }

    }
    else {
        foreach my $id_line ( 0 .. $#PU_measure ) {
            if ( substr( $PU_measure[$id_line], 0, 1 ) ne '#' ) {
                $max_dom = substr( $PU_measure[$id_line], 0, 2 )
                  if ( $max_dom == 0 );

                if ( substr( $PU_measure[$id_line], 0, 2 ) < $max_dom ) {
                    for ( my $i = 1 ; $i <= $alt_b + 1 ; $i++ ) {
                        if (
                            substr( $PU_measure[ ( $id_line - $i ) ], 0, 2 ) ==
                            $max_dom )
                        {
                            if ( $option_alt_dist == 1 ) {
                                my @temp =
                                  split( '\|',
                                    $PU_measure[ ( $id_line - $i ) ] );
                                my $abs = 1;
                                if ( DistanceModel( $temp[3], $temp[5], $abs ) <
                                    $max_dist )
                                {
                                    push( @relevant_measure,
                                        $PU_measure[ ( $id_line - $i ) ] );
                                }
                            }
                            else {
                                push( @relevant_measure,
                                    $PU_measure[ ( $id_line - $i ) ] );
                            }
                        }
                    }
                    $max_dom = substr( $PU_measure[$id_line], 0, 2 );
                }
            }
        }
    }

    push( @relevant_measure, $PU_measure[$#PU_measure] );

    return @relevant_measure;
}

1;

