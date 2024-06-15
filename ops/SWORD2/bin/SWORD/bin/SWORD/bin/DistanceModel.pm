
sub DistanceModel
{

	use List::Util qw[min max];
	my $Mx = shift;
	my $My = shift;
	my $abs= shift;

	my $minX = 0.008501;
	my $maxX = 1.099581;
	my $minY = 1.797981;
	my $maxY = 4.124218;

	#MIN/MAX normalization
	$Mx = ($Mx-$minX)/($maxX-$minX);
	$My = ($My-$minY)/($maxY-$minY);
	
	my $HORIZONTAL 			= 0.5884363 ;
	my $VERTICAL			= 0.2046999 ;
	my $DIAGONAL_intercept 	        = 0.4474396 ;
	my $DIAGONAL_slope 		= 1.680319 ;
	my $DIAGONAL_vertical	        = 0.2075470 ;
	my $DIAGONAL_horizontal	        = 0.7867766;

	#(abs(1.680319*0.1-0.8+0.4474396))/(sqrt(1+1.680319*1.680319))
	my $distance_to_diag       = abs($DIAGONAL_slope*$Mx-$My+$DIAGONAL_intercept)/sqrt(1+$DIAGONAL_slope*$DIAGONAL_slope);
	my $distance_to_horizontal = abs($My - $HORIZONTAL);
	my $distance_to_vertical   = abs($Mx - $VERTICAL);


	#Where ?
	# +----------+
	# |   3| 5 
	# | 2 /|----- 
	# |4 /1| 8
	# |----+-----
	# |6| 9|  7
	# | |  |
	# if ( $Ay < $HORIZONTAL )                     => zone 6
	# if ( $Ax > $VERTICAL   )                     => zone 5
	# if ( $Ay < $HORIZONTAL and $Ax > $VERTICAL ) => zone 7
	# if ( $Ax < $DIAGONAL_horizontal and $Ax < $VERTICAL ) => zone 5  
	# if ( $Ay > $HORIZONTAL and $Ax < $VERTICAL and $distance_to_diag < 0) => zone 1


# On which side of the line y=ax+b (A and B: dots on the line) is M(X,Y)?
# 0.5525604/1.680319
# y=ax+b A=(0,0.4474396) and B=(0.3288426,1) 
# Sign of D gives position relative to the line
# D < 0 right
# D > 0 left
# D=(xB-xA)*(yM-yA)-(yB-yA)*(xM-xA)

	my $Ax=0;
	my $Ay=$DIAGONAL_intercept;

	my $Bx=0.3288426;
	my $By=1;

	my $D=($Bx-$Ax)*($My-$Ay)-($By-$Ay)*($Mx-$Ax);
	my $sign=0;
	my $distance=0;
	#print "$D\n";	
	if    ( $My > $HORIZONTAL and $Mx < $VERTICAL and $D < 0)
	{
		# => ZONE 1
		#print "ZONE 1 (x=$Mx,y=$My)\n";
		$distance=min($distance_to_diag,$distance_to_horizontal,$distance_to_vertical);
		$sign=-1;
		return($distance*$sign);		
	}
	if    ( 
		($My < $DIAGONAL_horizontal and $My > $HORIZONTAL        and $D < 0) or 
		($Mx > $DIAGONAL_vertical   and $Mx < $VERTICAL          and $D < 0) or 
		($My < $DIAGONAL_horizontal and $Mx > $DIAGONAL_vertical and $D < 0)
		)
	{
		# => ZONE 7, 8, 9
		#print "ZONE 7, 8, 9 (x=$Mx,y=$My)\n";
		$distance=$distance_to_diag;
		$sign=-1;
		return($distance*$sign);		
	}

	elsif ( $My < $HORIZONTAL )
	{
		# => ZONE 6
		#print "ZONE 6 (x=$Mx,y=$My)\n";
		$distance=$distance_to_horizontal;
		$sign=-1;
		return($distance*$sign);

	}
	elsif ( $Mx > $VERTICAL   )
	{
		# => ZONE 5
		#print "ZONE 5 (x=$Mx,y=$My)\n";
		$distance=$distance_to_vertical;
		$sign=-1;
		return($distance*$sign);
	}
	else
	{
		# => ZONE 2 ou 3 ou 4
		#print "INSIDE (x=$Mx,y=$My)\n";
		$distance=min($distance_to_diag,$distance_to_horizontal,$distance_to_vertical);
		$sign=1;
		return($distance*$sign); 
	}



	
}

1;
