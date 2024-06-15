#! /usr/bin/perl -w
$dssp2pdb_version = "0.03" ;
# dssp2pdb
# converts dssp type output to HELIX/SHEET lines for PDB files
# distributed under gnu public license
# copyright James Stroud 2002

#version history
#0.01
#  original version
#0.02
#  added "-3" for 3(10) helices
#  added "-5" for pi helices
#  allowed for modified amino acids ("X" in dssp output)
#  fixed bug that choked on blank chain label in dssp output
#0.03
#  dispenses helix/sheet lines of original pdb

$true = 1 ;
$false = 0 ;

$helix = "HELIX" ;
$sheet = "SHEET" ;

use Getopt::Std ;
$opt_3 = $false ;
$opt_5 = $false ;
getopts ('35') ;
$helices_310 = $opt_3 ;  # flag for 3(10) helices
$helices_pi = $opt_5 ;   # flag for pi helices

sub printhelp() {
  print STDERR <<EOF ;

************************************************************************
** dssp2pdb version $dssp2pdb_version (c) James Stroud, 2002
** provided under the gnu public license
** james.stroud\@colorado.edu
************************************************************************
** Converts dssp file to HELIX/SHEET type lines for PDB files
** dssp2pdb will insert HELIX/SHEET lines above ATOM lines if
** "pdbfile.pdb "is provided
** When the "-3" flag is set, 3(10) helices will result in HELIX lines
** When the "-5" flag is set, pi helices will result in HELIX lines
************************************************************************
** Usage: dssp2pdb [-35] dsspfile.dssp [pdbfile.pdb] > new_pdbfile.pdb
************************************************************************

EOF
}

if (! $ARGV[0]) {
  printhelp() ;
  exit ;
}



%iupac = ( A => ALA,
           C => CYS,
           D => ASP,
           E => GLU,
           F => PHE,
           G => GLY,
           H => HIS,
           I => ILE,
           K => LYS,
           L => LEU,
           M => MET,
           N => ASN,
           P => PRO,
           Q => GLN,
           R => ARG,
           S => SER,
           T => THR,
           V => VAL,
           W => TRP,
           Y => TYR,
           X => UNK, ) ;  # UNK courtesy of Dirk Kostrewa

sub printss() {
  print @comments ;
  $number = 0 ;
  $helixcount = 0 ;
  $sheetcount = 0 ;
  while ($number < $elements) {
    if ($eletype[$number] eq $helix) {
      $helixcount++ ;
      $length = $resnum2[$number] - $resnum1[$number] + 1 ;
      ( defined $restype1[$number]) || ( $restype1[$number] = "UNK");
      ( defined $restype2[$number]) || ( $restype2[$number] = "UNK");
      printf "HELIX %4d %3d $restype1[$number] $reschain1[$number]%5d  ".
             "$restype2[$number] $reschain2[$number]%5d%3d".
             "                                 %3d\n",
             ($helixcount, $helixcount, $resnum1[$number], $resnum2[$number], 1, $length) ;
    }
    if ($eletype[$number] eq $sheet) {
      $sheetcount++ ;
      ( defined $restype1[$number]) || ( $restype1[$number] = "UNK");
      ( defined $restype2[$number]) || ( $restype2[$number] = "UNK");
      printf "SHEET %4d %3d 1 $restype1[$number] $reschain1[$number]%4d  ".
             "$restype2[$number] $reschain2[$number]%4d  0\n",
             ($sheetcount, $sheetcount, $resnum1[$number], $resnum2[$number]) ;
    }
    $number++ ;
  }
}


$dsspfile = $ARGV[0] ;

$comments[0] = "REMARK   1 secondary structure assigned by dssp\n" ;
$comments[1] = "REMARK   2 dssp output converted by dssp2pdb version $dssp2pdb_version\n" ;

open (DSSP, $dsspfile) or die "Could not open file \"$dsspfile\"\n" ;

$helixcount = 0 ;
$sheetcount = 0 ;

#
# loop to find start of ss info
#
while () {
  $line = <DSSP> ;
  if (! $line) { 
    printhelp() ;
    print STDERR "\n\"$dsspfile\" does not seem to be dssp output\n\n" ;
    exit ;
  }
  if ($line =~ /RESIDUE AA STRUCTURE/) {   # ready to cull ss info
    last ;
  }
}

#
# loop to cull ss info from dssp file
#
$residues = 0 ;    # will eventually become number of residues in file
while () {
    $line = <DSSP> ;
    if (! $line) { 
    last ;
    }
    $resid[$residues] = substr ($line, 5, 5) ;
    $chain[$residues] = substr ($line, 11, 1) ;
    $aa[$residues] = $iupac{substr ($line, 13, 1)} ;
    $ss[$residues] = substr ($line, 16, 1) ;
    $residues++ ;
}

close (DSSP) ;

#
# loop to parse ss info from data structures
#
$elements = 0 ;       # index to keep track of ss elements we are on
$helixcount = 0 ;    # number of helices
$sheetcount = 0 ;    # number of sheets
$number = 0 ;        # residue number, index of this loop
$last_type = "" ;
while($number < $residues) {
    $current_type = "" ; # type of ss elements, will have value only if HELIX or SHEET
    if ($ss[$number] eq "H") {
    $current_type = $helix ;
    }
    if ($ss[$number] eq "E") {
    $current_type = $sheet ;
    }
    if ($ss[$number] eq "G" && $helices_310) {
    $current_type = $helix ;
    }
    if ($ss[$number] eq "I" && $helices_pi) {
    $current_type = $helix ;
    }
    if ($current_type eq $helix) {
    if ($current_type eq $last_type) {
        # extension of helix
        $restype2[$elements-1] = $aa[$number] ;
        $resnum2[$elements-1] = $resid[$number] ;
        $reschain2[$elements-1] = $chain[$number] ;
    } else {
        # new helix
        $helixcount++ ;
        $eletype[$elements] = $helix ;
        $restype1[$elements] = $aa[$number] ;
        $resnum1[$elements] = $resid[$number] ;
        $reschain1[$elements] = $chain[$number] ;
        $restype2[$elements] = $aa[$number] ;
        $resnum2[$elements] = $resid[$number] ;
        $reschain2[$elements] = $chain[$number] ;
        $elements++ ;
    }
  }
  if ($current_type eq $sheet) {
    if ($current_type eq $last_type) {
      # extension of strand
      $restype2[$elements-1] = $aa[$number] ;
      $resnum2[$elements-1] = $resid[$number] ;
      $reschain2[$elements-1] = $chain[$number] ;
    } else {
      # new strand
      $sheetcount++ ;
      $eletype[$elements] = $sheet ;
      $restype1[$elements] = $aa[$number] ;
      $resnum1[$elements] = $resid[$number] ;
      $reschain1[$elements] = $chain[$number] ;
      $restype2[$elements] = $aa[$number] ;
      $resnum2[$elements] = $resid[$number] ;
      $reschain2[$elements] = $chain[$number] ;
      $elements++ ;
    }
  }
  $last_type = $current_type ;
  $number++ ;
}

$first_atom = $true ;
if ($ARGV[1]) {
  $pdbfile = $ARGV[1] ;
  open (PDB, $pdbfile) or die "Could not open file \"$pdbfile\"\n" ;
  #
  # loop to print pdb file lines before ATOM lines
  #
  while () {
    $line = <PDB> ;
    if (! $line) {
      if ($first_atom) {
        printhelp() ;
        print STDERR "\n\"$pdbfile\" does not seem like a pdb file\n\n" ;
        exit ;
      } else { exit }
    }
    if ($line =~ /^ATOM/ and $first_atom) {    # print out ss info
      printss() ;
      $first_atom = $false ; 
    }
    if ($line !~ /^HELIX/ and $line !~/^SHEET/) {
      print $line ;
    }
  }
} else { 
  printss() ;
}


exit ;
