package ParsePDB;

#######################################################################################################################
#
# Package:         ParsePDB
#
# Function:        Reads in a given PDB, can extract single models or chains and a lot more...
#
# Author:          Benjamin Bulheller
#
# Website:         http://comp.chem.nottingham.ac.uk/parsepdb/
#                  http://www.bulheller.com
#
# Mail address:    webmaster.-at-.bulheller.com
#
# Research Group:  Prof. Dr. Jonathan Hirst
#                  School of Physical Chemistry
#                  University of Nottingham
#
# Funded by:       EPSRC
#
# Date:            November 2005 - Feburary 2009
#
# Revision:        $Revision: 5084 $, $Date: 2011-02-13 21:04:35 +0000 (Sun, 13 Feb 2011) $
#
# Acknowledgments: Special thanks to Dr. Daniel Barthel for many, many discussions and help whenever needed!
#
# Licence:         This program is free software: you can redistribute it and/or modify
#                  it under the terms of the GNU General Public License as published by
#                  the Free Software Foundation, either version 3 of the License, or
#                  (at your option) any later version.
#
#                  This program is distributed in the hope that it will be useful,
#                  but WITHOUT ANY WARRANTY; without even the implied warranty of
#                  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#                  GNU General Public License for more details.
#
#                  You should have received a copy of the GNU General Public License
#                  along with this program.  If not, see http://www.gnu.org/licenses/.
#
#######################################################################################################################

use strict;                          # always use this!!!
use Data::Dumper;                    # for easy printout of arrays and hashes
require Error;                  # for proper exception handling

our $VERSION = "2.0";                # for version determination of the MakeFileMaker


#######################################################################################################################
# Default Values
#######################################################################################################################

our $ChainLabelAsLetterDefault = 0;     # export single chains with letter (-cA) or with numbers (-c2)
our $ChainSuffixDefault        = "-c";  # added to the file name by WriteChains
our $ModelSuffixDefault        = "-m";  # added to the file name by WriteModels
our $HeaderRemarkDefault       = 1;     # whether write should add a remark before the header
our $ThrowErrorsDefault        = 1;     # whether errors should be throws (i.e. the program dies when an error happens)
our $AtomLocationsDefault      = "All"; # if alternate atom locations are available, all will be returned
our $NoANISIGDefault           = 0;     # filter out all SIGATM, SIGUIJ and ANISOU lines on|off
our $NoHETATMDefault           = 0;     # filter all HETATM lines on|off
our $VerboseDefault            = 1;     # whether all warnings and messages are printed or not

# This is not yet implemented!
our $DetectChainBreakDefault   = 0;     # whether excessive C-C distances shall require a chain break

# defines the line after which the HeaderRemark is added and the returned lines, if MinHeader or MinFooter is requested
our $HeaderRemarkLine = "COMPND";
our $MinHeaderLines = "HEADER|TITLE|COMPND";
our $MinFooterLines = "END";

our %WARNINGS = (
	"UnknownModel"       => "The following model is not defined: ",
	"UnknownChain"       => "The following chain is not defined: ",
	"UnknownResidue"     => "The following residue is not defined: ",
	"UnknownAtom"        => "The following atom is not defined: ",
	"UnknownLabel"       => "The following label has not been found: ",
	"UnknownNumber"      => "The following number has not been found: ",
	"UnknownType"        => "The following type has not been recognized: ",
	"UnknownElement"     => "The following element has not been recognized: ",
	"ModelDefault"       => "WARNING: Model not specified. Taking model 0 as default!",
	"ChainDefault"       => "WARNING: Chain not specified. Taking chain 0 as default!",
	"NoChainLabel"       => "No chain identifier given! This can be corrected using the method \->RenumberChains!",
	"NoENDMDL"           => "No \"ENDMDL\" termination of a \"MODEL\"-block has been found.",
	"MultipleChainLabel" => "Multiple chain ID! The following chain ID has been found more than once: ",
	"ChainLabelsInvalid" => "Multiple or undefined chain labels have been found. Either correct the chain labels using ->RenumberChains or use internal identifiers.",
	"UnknownAminoAcid"   => "The following code has not been recognized: ",
	"InvalidAminoAcid"   => "Abbreviation not recognized! Only 1- or 3-letter codes can be converted!",
   "UnknownParameter"   => "Unknown parameter, please check for a typo: ",
);

our %ERRORS = (
	"NoFile"             => "No file name given!",
	"IOError"            => "Cannot write file: ",
	"CorruptFile"        => "The file is not valid, no ATOM lines have been found: ",
	"FileNotFound"       => "The file has not been found: ",
	"UnknownElement"     => "The following element has not been recognized: ",
	"BadParameter"       => "Given parameters cannot be combined: ",
	"ParamNumeric"       => "The following Parameter must be numeric: ",
);

# the keys to the elements have to match the ones in ELEMENTMASS! 
our %ELEMENTFILTERS = (
	C  => "(C|C[A-Z]|C[1-9]\*?|C[A-Z][1-9])",
	O  => "(O|O[A-Z]|O[1-9]\*?|O[A-Z][1-9]|OXT)",
	N  => "(N|N[A-Z]|N[1-9]\*?|N[A-Z][1-9])",
	H  => "(H|H[A-Z]|[1-9]H[1-9]\*?|[1-9]?H[A-Z][1-9]?)",
	S  => "(S|SD|SG)",
	P  => "(P|P[AB])",
	D  => "(D|D[1-9])",
	CU => "(CU|CU[1-9])",
	FE => "(FE|FE[1-9])",
	MG => "(MG|MG[1-9])",
	ZN => "(ZN|ZN[1-9])",
);

our %ELEMENTMASS = (
	"H"  => 1.00794,    "D"  => 2.014101,   "T"  => 3.016049,   "HE" => 4.002602,   "LI" => 6.941,
	"BE" => 9.012182,   "B"  => 10.811,     "C"  => 12.0107,    "N"  => 14.00674,   "O"  => 15.9994,
	"F"  => 18.9984032, "NE" => 20.1797,    "NA" => 22.989770,  "MG" => 24.3050,    "AL" => 26.981538,
	"SI" => 28.0855,    "P"  => 30.973761,  "S"  => 32.066,     "CL" => 35.4527,    "AR" => 39.948,
	"K"  => 39.0983,    "CA" => 40.078,     "SC" => 44.955910,  "TI" => 47.867,     "V"  => 50.9415,
	"CR" => 51.9961,    "MN" => 54.938049,  "FE" => 55.845,     "CO" => 58.933200,  "NI" => 58.6934,
	"CU" => 63.546,     "ZN" => 65.39,      "GA" => 69.723,     "GE" => 72.61,      "AS" => 74.92160,
	"SE" => 78.96,      "BR" => 79.904,     "KR" => 83.80,      "RB" => 85.4678,    "SR" => 87.62,
	"Y"  => 88.90585,   "ZR" => 91.224,     "NB" => 92.90638,   "MO" => 95.94,      "TC" => 98.0,
	"RU" => 101.07,     "RH" => 102.90550,  "PD" => 106.42,     "AG" => 107.8682,   "CD" => 112.411,
	"IN" => 114.818,    "SN" => 118.710,    "SB" => 121.760,    "TE" => 127.60,     "I"  => 126.90447,
	"XE" => 131.29,     "CS" => 132.90545,  "BA" => 137.327,    "LA" => 138.9055,   "CE" => 140.116,
	"PR" => 140.90765,  "ND" => 144.24,     "PM" => 145.0,      "SM" => 150.36,     "EU" => 151.964,
	"GD" => 157.25,     "TB" => 158.92534,  "DY" => 162.50,     "HO" => 164.93032,  "ER" => 167.26,
	"TM" => 168.93421,  "YB" => 173.04,     "LU" => 174.967,    "HF" => 178.49,     "TA" => 180.9479,
	"W"  => 183.84,     "RE" => 186.207,    "OS" => 190.23,     "IR" => 192.217,    "PT" => 195.078,
	"AU" => 196.96655,  "HG" => 200.59,     "TL" => 204.3833,   "PB" => 207.2,      "BI" => 208.98038,
	"PO" => 209.0,      "AT" => 210.0,      "RN" => 222.0,      "FR" => 223.0,      "RA" => 226.0,
	"AC" => 227.0,      "TH" => 232.038,    "PA" => 231.03588,  "U"  => 238.0289,   "NP" => 237.0,
	"PU" => 244.0,      "AM" => 243.0,      "CM" => 247.0,      "BK" => 247.0,      "CF" => 251.0,
	"ES" => 252.0,      "FM" => 257.0,      "MD" => 258.0,      "NO" => 259.0,      "LR" => 262.0,
	"RF" => 261.0,      "DD" => 262.0,      "SG" => 266.0,      "BH" => 264.0,      "HS" => 269.0,
	"MT" => 268.0,
);


#######################################################################################################################
# Object constructor
#######################################################################################################################

sub new { # create a new object
	my ($invocant, $class, $self);
	
	$invocant = shift; # either class name or object
	$class = ref($invocant) || $invocant; # object or class name
	
	$self = {};
	bless ($self, $class); # bless $self into class $class
	$self->init (@_);
	
	return $self;
} # of sub new


sub init { # initialize the global variables of the object
	my $self = shift;
	my (%args, $FileName, $BaseName);
	
	%args = @_;
	
	$FileName = $args{FileName};
	
	
	####################################################################################################################
	# initialize the object variables
	####################################################################################################################
	
	$self->{FileName}     = $FileName;        # the file name of the PDB
	$self->{Header}       = [];               # complete header until first MODEL or ATOM
	$self->{Footer}       = [];               # complete footer from the last ENDMDL or TER till the end
	$self->{Content}      = [];               # complete content of the PDB file including header and footer (but AtomLocations, NoANISIG and NoHETATM do affect it)
	$self->{Models}       = ();               # hash for the separated groups of the PDB (models, chains, atoms)
	
	$self->{AllAtoms}           = undef;      # contains the RegEx which line starts are atoms (ATOM|HETATM|...)
	$self->{ChainLabelAsLetter} = undef;      # export single chains with letter (-cA) or with numbers (-c2)
	$self->{ChainSuffix}        = undef;      # added to the file name by WriteChains
	$self->{ModelSuffix}        = undef;      # added to the file name by WriteModels
	$self->{HeaderRemark}       = undef;      # whether write should add a remark before the header
	$self->{AtomLocations}      = undef;      # if alternate atom locations are available, all will be returned
	$self->{NoANISIG}           = undef;      # filter out all SIGATM, SIGUIJ and ANISOU lines on|off
	$self->{NoHETATM}           = undef;      # filter all HETATM lines on|off
	$self->{Verbose}            = undef;      # whether all warnings and messages are printed or not
	
	$self->{WarningMsg} = [];                 # warnings that errors have been found which have been corrected
	$self->{Warning} = ();                    # hash with true/false values for "if" checks
	$self->{Warning}{Warning} = 0;            # whether any warning has been issued
	$self->{Warning}{UnknownModel} = 0;       # if the given model is not defined
	$self->{Warning}{UnknownChain} = 0;       # if the given chain is not defined
	$self->{Warning}{UnknownResidue} = 0;     # if the given residue is not defined
	$self->{Warning}{UnknownAtom} = 0;        # if the given atom is not defined
	$self->{Warning}{UnknownLabel} = 0;       # if the given label could not be found in the PDB
	$self->{Warning}{NoENDMDL} = 0;           # if a MODEL without a corresponding ENDMDL has been found
	$self->{Warning}{NoChainLabel} = 0;       # if no chain ID is given at all
	$self->{Warning}{MultipleChainLabel} = 0; # if a chain ID has been found more than once in a model
	$self->{Warning}{ChainLabelsInvalid} = 0; # if it was tried to access chains with external identifieres and invalid chain labels
	$self->{Warning}{UnknownAminoAcid} = 0;   # if a 1- or 3-letter-code given to AminoAcidConvert has not been recognized
	$self->{Warning}{InvalidAminoAcid} = 0;   # if a code given to AminoAcidConvert has not 1 or 3 letters
	$self->{ReportNoWarning} = 0;             # this is only used internally, if controlled actions are taken which might produce warnings
	
	$self->{ErrorMsg} = undef;                # error message a of fatal error that caused the process to be aborted
	
	
	####################################################################################################################
	# read parameters or set to default values
	####################################################################################################################
	
	if (defined $args{ChainLabelAsLetter}) {
		$self->{ChainLabelAsLetter} = $args{ChainLabelAsLetter};
	}
	else {
		$self->{ChainLabelAsLetter} = $ChainLabelAsLetterDefault;
	}
	
	if (defined $args{ChainSuffix}) {
		$self->{ChainSuffix} = $args{ChainSuffix};
	}
	else {
		$self->{ChainSuffix} = $ChainSuffixDefault;
	}
	
	if (defined $args{ModelSuffix}) {
		$self->{ModelSuffix} = $args{ModelSuffix};
	}
	else {
		$self->{ModelSuffix} = $ModelSuffixDefault;
	}
	
	if (defined $args{HeaderRemark}) {
		$self->{HeaderRemark} = $args{HeaderRemark};
	}
	else {
		$self->{HeaderRemark} = $HeaderRemarkDefault;
	}
	
	if (defined $args{AtomLocations}) {
		$self->{AtomLocations} = uc ($args{AtomLocations});
	}
	else {
		$self->{AtomLocations} = uc ($AtomLocationsDefault);
	}
	
	if (defined $args{NoANISIG}) {
		$self->{NoANISIG} = $args{NoANISIG};
	}
	else {
		$self->{NoANISIG} = $NoANISIGDefault;
	}
	
	if (defined $args{NoHETATM}) {
		$self->{NoHETATM} = $args{NoHETATM};
	}
	else {
		$self->{NoHETATM} = $NoHETATMDefault;
	}
	
	if (defined $args{Verbose}) {
		$self->{Verbose} = $args{Verbose};
	}
	else {
		$self->{Verbose} = $VerboseDefault;
	}
	
	if    ( $self->{NoHETATM} and $self->{NoANISIG}) {
		$self->{AllAtoms} = "ATOM";
	}
	elsif ( not $self->{NoHETATM} and $self->{NoANISIG} ) {
		$self->{AllAtoms} = "ATOM|HETATM";
	}
	elsif ( $self->{NoHETATM} and not $self->{NoANISIG} ) {
		$self->{AllAtoms} = "ATOM|SIGATM|SIGUIJ|ANISOU";
	}
	elsif ( not $self->{NoHETATM} and not $self->{NoANISIG} ) {
		$self->{AllAtoms} = "ATOM|HETATM|SIGATM|SIGUIJ|ANISOU";
	}

	# This is not yet implemented
	if (defined $args{DetectChainBreak}) {
		$self->{DetectChainBreak} = $args{DetectChainBreak}
	}
	else {
		$self->{DetectChainBreak} = $DetectChainBreakDefault;
	}
	
	$self->{ANISIG} = "ANISOU|SIGATM|SIGUIJ";
	
	####################################################################################################################
	# check for file
	####################################################################################################################
	
	if (not defined $FileName) { # if no file name was given
		$self->_ReportError ("NoFile");
	}
	
	if (not -f $FileName) {       # if the given file is not found
		if (-f "$FileName.pdb") {  # try with additional extension
			$self->{FileName} = "$FileName.pdb";
		}
		else {
			$self->_ReportError ("FileNotFound", $FileName);
		}
	}
	
	$BaseName = $self->_GetBaseName ($FileName);
	$self->{BaseName} = $BaseName;
} # sub init


#######################################################################################################################
#######################################################################################################################
## "Public" Methods
#######################################################################################################################
#######################################################################################################################

sub Read { # reads the whole PDB file into the content array of the object and checks for several possible problems
	my $self = shift;
	my %args = @_;
	my ($ParsePDB, $LineCount, $AtomFound);
	my ($Line, $FileName, @Content, @ANISIG, @HETATM);
	
	# if Read was called by Parse then ParsePDB will be defined and
	# therefore be set false
	# if Read was called on its own, then ParsePDB will be set true
	# and therefore the read content be reparsed
	if (defined $args{ParsePDB}) { $ParsePDB = 0 }
	                        else { $ParsePDB = 1 }
	
	$FileName = $self->{FileName};
	
	$AtomFound = 0;
	
	open (PDB, "<$FileName") or $self->_ReportError ("IOError", $FileName);
	while ( $Line = <PDB> ) {
		# make sure that the line is 80 charcters in length (more or less standard
		# for PDB files and to avoid "substr outside of string"-errors)
		if (length $Line < 80) {
			# remove the newline character at the end of the line
			chomp $Line;
			
			# Remove a carriage return. If this is omitted, ^M characters are left over if later on
			# \n is added when the file is written out
			$Line =~ s/\r//;
			
			until (length $Line > 79) { $Line = $Line . " " }
			$Line = $Line . "\n";
		}
		
		# see, whether there is at least one atom in the file
		if ( (not $AtomFound) and ($Line =~ m/^($self->{AllAtoms})/) ) { $AtomFound = 1 }
		if ($self->{NoHETATM} and ($Line =~ m/^HETATM/) ) { next }
		if ($self->{NoANISIG} and ($Line =~ m/^(SIGATM|SIGUIJ|ANISOU)/) ) { next }
		
		push @Content, $Line;
	}
	close (PDB);
	
	if ( (not @Content) or (not $AtomFound) ) {
		$self->_ReportError ("CorruptFile", $FileName);
	}
	
	$self->{Content} = \@Content; # save the main content
	if ($ParsePDB) { $self->Parse }
} # of sub Read


sub Reset { # reads the file again (after having done things one likes to make unhappened)
	my $self = shift;
	
	$self->Read;
} # of sub Reset


sub Parse { # parses the PDB and sorts all information into the {Models} hash
	my $self = shift;
	
	my ($Content, @Header, @Footer);
	my (%Models, @AllModels, @Chains, $FirstModel);
	my ($LineCount, $ModelCount, $ModelEntries, $ModelStart, $ModelEnd);
	my ($ChainCount, $CurChainLabel, $LastChainLabel, $NextChainLabel);
	
	
	####################################################################################################################
	# reset all warning messages
	####################################################################################################################
	
	$self->{Warning} = ();                     # hash with true/false values for "if" checks
	$self->{Warning}{Warning} = 0;             # whether any warning has been issued
	$self->{Warning}{UnknownModel} = 0;        # if the given model is not defined
	$self->{Warning}{UnknownChain} = 0;        # if the given chain is not defined
	$self->{Warning}{UnknownResidue} = 0;      # if the given residue is not defined
	$self->{Warning}{UnknownAtom} = 0;         # if the given atom is not defined
	$self->{Warning}{UnknownLabel} = 0;        # if the given label could not be found in the PDB
	$self->{Warning}{NoENDMDL} = 0;            # if a MODEL without a corresponding ENDMDL has been found
	$self->{Warning}{NoChainLabel} = 0;        # if no chain ID is given at all
	$self->{Warning}{MultipleChainLabel} = 0;  # if a chain ID has been found more than once in a model
	$self->{Warning}{ChainLabelsInvalid} = 0; # if it was tried to access chains with external identifieres and invalid chain labels
	$self->{Warning}{UnknownAminoAcid} = 0;    # if a 1- or 3-letter-code given to AminoAcidConvert has not been recognized
	$self->{Warning}{InvalidAminoAcid} = 0;    # if a code given to AminoAcidConvert has not 1 or 3 letters
	
	####################################################################################################################
	
	# if the content array is empty, read the file
	#if (not defined @{$self->{Content}}) { $self->Read (ParsePDB => 0) }
	if (!@{$self->{Content}}) { $self->Read (ParsePDB => 0) } # G.Postic (2014): use of defined on aggregates is deprecated
	
	$Content = $self->{Content}; # take the complete content of the file
	
	
	####################################################################################################################
	# read the header of the file and save it
	####################################################################################################################
	
	$LineCount = 0;
	
	# until the first line that begins with MODEL or ATOM	
	until ($Content->[$LineCount] =~ m/^(MODEL|ATOM)/) {
		$Header[$LineCount] = $Content->[$LineCount];     # save the header into @Header
		++$LineCount;
		
		if ($LineCount > $#{$Content}) {
			$self->_ReportError ("CorruptFile", $self->{FileName});
		}
	}
	
	$self->{Header} = \@Header; # save the header into the object
	
	
	####################################################################################################################
	# sort the atom section into $self->{Models}
	####################################################################################################################
	
	$LastChainLabel = "something"; # to have it different from the first chain label
	$ChainCount = -1;              # this is increased when the first chain is started
	
	# The current line of the PDB contains either MODEL or ATOM. If it is MODEL, then the $ModelCount will be increase
	# at first and hence has to be set to -1 to reach number 0 in the first loop iteration. If the first line is ATOM,
	# then there are no MODEL tags and the $ModelCount is set 0 at first. 
	if ($Content->[$LineCount] =~ m/^MODEL/) { $ModelCount = -1 }
	                                    else { $ModelCount =  0 }

	# The following loop splits the chains up and saves the content of each chain and the chain label
	# until the end of the array or the first CONECT, MASTER or END line
	until ( ($LineCount > $#{$Content}) or ($Content->[$LineCount] =~ m/^(CONECT|MASTER|END |END$)/) ) {
		if ($Content->[$LineCount] =~ m/^MODEL/) {
			++$ModelCount;

			$self->{Models}{$ModelCount}{MODEL} = $Content->[$LineCount];
			$ChainCount = -1;
			$LastChainLabel = "something";
		}
		
		elsif ($Content->[$LineCount] =~ m/^ENDMDL/) {
			$self->{Models}{$ModelCount}{ENDMDL} = $Content->[$LineCount];
			$ChainCount = -1;
			$LastChainLabel = "something";
		}
		
		elsif ($Content->[$LineCount] =~ m/^($self->{AllAtoms}|TER)/) {
			# The chain ID of TER lines is ignored because some/many PDBs don't add the current ID to the TER
			# line (e.g. the line is only "TER\n". In such a case the TER could end up in a new chain, thus for
			# the determination whether a new chain has to be started, the IDs before and after the TER are taken.
			if ($Content->[$LineCount] !~ m/^TER/) {
				$CurChainLabel  = substr ($Content->[$LineCount], 21, 1);
			}
			
			# a changed chain ID forces the beginning of a new chain (this works also without a TER)
			if ($CurChainLabel ne $LastChainLabel) {
				++$ChainCount;
				
				$self->{Models}{$ModelCount}{$ChainCount}{ChainLabel} = $CurChainLabel;
				$self->{Models}{$ModelCount}{$ChainCount}{ChainLabel} =~ s/\s//g; # return an empty string (false) if no chain ID is found
				
				$LastChainLabel = $CurChainLabel;
			}
			
			# A TER provokes a chain change, as long as the chain ID (if defined) changes.
			# Some subroutines rely on the fact that a TER can *only* be at the very end of a chain, changing
			# this would require looping over a chain to check for a TER and thus slow down parsing.
			if ($Content->[$LineCount] =~ m/^TER/) {
				if ($Content->[$LineCount+1] and $Content->[$LineCount+1] =~ m/^($self->{AllAtoms})/) { # if an atom follows the TER
					$NextChainLabel = substr ($Content->[$LineCount+1], 21, 1);
					
					# If there are chain labels (i.e. the label is not a blank) and the chain label after the TER
					# is same as before it, then the file format is crap. Real crap.
					# The TER is DISCARDED in that case
					if ( ($CurChainLabel ne " ") and ($NextChainLabel eq $CurChainLabel) ) {
						# do nothing here, the TER is discarded
					}
					else {
						# If there are no chain labels or the next one differs from the current, a new chain is started.
						# This works also, if the current label is defined and the next one is a blank (e.g. HETATMs)
						$LastChainLabel = "something";
						push @{$self->{Models}{$ModelCount}{$ChainCount}{Content}}, $Content->[$LineCount];
					}
				}
				else { # if no atom is following the TER, e.g. an ENDMDL
					$LastChainLabel = "something";
					push @{$self->{Models}{$ModelCount}{$ChainCount}{Content}}, $Content->[$LineCount];
				}
			}
			else { # if the line is not a TER
				push @{$self->{Models}{$ModelCount}{$ChainCount}{Content}}, $Content->[$LineCount];
			}
		}
		else { # if the line does not begin with MODEL, $self->{AllAtoms}, TER or ENDMDL
			chomp $Content->[$LineCount];
			
			# read the line, remove all blanks to "allow" empty lines without an error message
			my $Line = $Content->[$LineCount];
			$Line =~ s/\s+$//g;
			if ($Line ne "") {
				$self->_ReportWarning ("Parse", "Misc", "Line \"$Content->[$LineCount]\" not recognized!");
			}
		}
		
		++$LineCount;
	} # of until ( ($LineCount == $#{$Content} or ($Content->[$LineCount] =~ m/^(CONECT|MASTER|END)/) )
	
	
	####################################################################################################################
	# read the footer of the file and save it
	####################################################################################################################
	
	# save everything from the last ENDMDL entry until the end
	until ($LineCount > $#{$Content}) {
		push @Footer, $Content->[$LineCount];
		++$LineCount;
	}
	
	$self->{Footer} = \@Footer;
	
	$self->_IndexChains;
#	$self->_IndexAtoms;
	
	$self->_DetermineAltLocIndicator ($Content);
	
	return 1;
} # of sub Parse


sub Get { # returns a specified part of the PDB or the whole PDB, if no part is specified
	my $self = shift;
	my %args = @_;
	my ($SetChainLabel, $AtomLocations, $KeepInsertions, $GetReference, $CurLine);
	my ($Race, $ModelStart, $ChainStart, $ResidueStart, $AtomStart, $Header, $Footer, $AtomIndex, $ResidueIndex);
	my ($Residue, $CurResidueNumber, $OldResidueNumber, $MinHeader, $MinFooter, $IgnoreTER);
	my ($ResidueNumber, $InsResidue);
	
	$self->_CheckParameters (%args); # track typos in the given paramters
	
	if (not defined $self->{Models}) { $self->Parse }
	
	
	####################################################################################################################
	# Read the configuration parameters
	####################################################################################################################
	
	if (defined $args{ModelStart})     { $ModelStart = $args{ModelStart} }
	                          else     { $ModelStart = undef }
	
	if (defined $args{ChainStart})     { $ChainStart = $args{ChainStart} }
	                          else     { $ChainStart = undef }
	
	if (defined $args{ResidueStart})   { $ResidueStart = $args{ResidueStart} }
	                            else   { $ResidueStart = undef }
	
	if (defined $args{AtomStart})      { $AtomStart = $args{AtomStart} }
	                         else      { $AtomStart = undef }
	
	if (defined $args{SetChainLabel})  { $SetChainLabel = $args{SetChainLabel} }
	                             else  { $SetChainLabel = undef }
	
	if (defined $args{KeepInsertions}) { $KeepInsertions = $args{KeepInsertions} }
	                              else { $KeepInsertions = 1 }
	
	if (defined $args{Header})         { $Header = $args{Header} }
	                      else         { $Header = 0 }
	
	if (defined $args{Footer})         { $Footer = $args{Footer} }
	                      else         { $Footer = 0 }
	
	if (defined $args{MinHeader})      { $MinHeader = $args{MinHeader} }
	                         else      { $MinHeader = 0 }
	
	if (defined $args{MinFooter})      { $MinFooter = $args{MinFooter} }
	                         else      { $MinFooter = 0 }
	
	if (defined $args{GetReference})   { $GetReference = $args{GetReference} }
	                            else   { $GetReference = 0 }
	
	if (defined $args{AtomIndex})      { $AtomIndex = $args{AtomIndex} }
	                         else      { $AtomIndex = 0 }
	
	if (defined $args{ResidueIndex})   { $ResidueIndex = $args{ResidueIndex} }
	                            else   { $ResidueIndex = 0 }
	
	if (defined $args{Race})           { $Race = $args{Race} }
	                    else           { $Race = "$self->{AllAtoms}" }
	
	if (defined $args{IgnoreTER})      { $IgnoreTER = $args{IgnoreTER} }
	                         else      { $IgnoreTER = 0 }
	
	if (defined $args{AtomLocations})  { $args{AtomLocations} = uc ($args{AtomLocations}) }
	                             else  { $args{AtomLocations} = uc ($self->{AtomLocations}) }
	
	
	####################################################################################################
	# check for contradicting arguments
	####################################################################################################
	
	if ( (defined $args{Model}) and (defined $args{ModelNumber}) ) {
		$self->_ReportError ("BadParameter", "Model and ModelNumber");
	}
	
	if ( (defined $args{Chain}) and (defined $args{ChainLabel}) ) {
		$self->_ReportError ("BadParameter", "Chain and ChainLabel");
	}
	
	if ( (defined $args{Residue}) and (defined $args{ResidueNumber}) ) {
		$self->_ReportError ("BadParameter", "Residue and ResidueNumber");
	}
	
	if ( (defined $args{Atom}) and (defined $args{AtomNumber}) ) {
		$self->_ReportError ("BadParameter", "Atom and AtomNumber");
	}
	
	
	####################################################################################################################
	# Retrieve the Content
	####################################################################################################################
	
	# initialize an array where the requested content is saved to
	my @Content;
	$args{Content} = \@Content;
	
	my @AtomIndex;
	$args{AtomIndex} = \@AtomIndex;
	
	my @ResidueIndex;
	if ($ResidueIndex) {
		$args{ResidueIndex} = \@ResidueIndex;
	}
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (%args)) {
		$self->_ReportWarning ("Get", "ChainLabelsInvalid");
		return undef;
	}
	
	if ( (defined $args{Atom}) or (defined $args{AtomNumber}) ) {
		$self->_RetrieveAtom (\%args);
	}
	elsif ( (defined $args{Residue}) or (defined $args{ResidueNumber}) ) {
		$self->_RetrieveResidue (\%args);
	}
	elsif ( (defined $args{Chain}) or (defined $args{ChainLabel}) ) {
		$self->_RetrieveChain (\%args);
		if ($#Content < 0) { return undef }
		
		# if the last line in the residue is TER, remove it
		if ($AtomIndex[$#AtomIndex]{Race} =~ m/TER/) {
			pop @AtomIndex;
			pop @Content;
		}
	}
	elsif ( (defined $args{Model}) or (defined $args{ModelNumber}) ) {
		$self->_RetrieveModel (\%args);
	}
	else {
		$self->_RetrieveAllModels (\%args);
	}
	
	if ($#Content < 0) { return undef }
	
	####################################################################################################################
	# Post-processing of the retrieved content
	####################################################################################################################
	
	
	if (defined $ModelStart)    { $self->RenumberModels   (Content      => \@Content,
	                                                       AtomIndex    => \@AtomIndex,
	                                                       ModelStart   => $ModelStart,   ) }
	
	if (defined $ChainStart)    { $self->RenumberChains   (Content      => \@Content,
	                                                       AtomIndex    => \@AtomIndex,
	                                                       ChainStart   => $ChainStart    ) }
	
	if (defined $ResidueStart)  { $self->RenumberResidues (Content        => \@Content,
	                                                       AtomIndex      => \@AtomIndex,
	                                                       ResidueStart   => $ResidueStart,
	                                                       KeepInsertions => $KeepInsertions) }
	
	if (defined $AtomStart)     { $self->RenumberAtoms    (Content      => \@Content,
	                                                       AtomIndex    => \@AtomIndex,
	                                                       AtomStart    => $AtomStart,
	                                                       IgnoreTER    => $IgnoreTER     ) }
	
	if ($SetChainLabel) { $self->SetChainLabel    (Content      => \@Content,
	                                               AtomIndex    => \@AtomIndex,
	                                               ChainLabel   => $SetChainLabel ) }
	
	if ($Header or $MinHeader) {
		$self->GetHeader (Content => \@Content, AtomIndex => \@AtomIndex, MinHeader => $MinHeader)
	}
	
	if ($Footer or $MinFooter) {
		$self->GetFooter (Content => \@Content, AtomIndex => \@AtomIndex, MinFooter => $MinFooter)
	}
	
	if ($args{CharmmFormat}) { $self->_PDB2CHARMM (Content => \@Content, AtomIndex => \@AtomIndex) }
	if ($args{PDB2CHARMM})   { $self->_PDB2CHARMM (Content => \@Content, AtomIndex => \@AtomIndex) }
	if ($args{CHARMM2PDB})   { $self->_CHARMM2PDB (Content => \@Content, AtomIndex => \@AtomIndex) }
	
	# $GetReference is checked first, so if $ResidueIndex oder $AtomIndex is given to any method which uses Get,
	# it will not interfere with that. However, ResidueIndex only works with Get at the moment
	if ($GetReference) {
		return \@AtomIndex, \@Content;
	}
	elsif ($AtomIndex) {
		return @AtomIndex;
	}
	elsif ($ResidueIndex) {
		# The atom lines have to be divided up now. Doing that at once, when the content is
		# filtered would be too extensive and slow down processing too much, since it would
		# have to be done in all renumber routines. This way, it is really only done if it is
		# requested
		
		$Residue = -1;                # will be increased for the first residue
		$OldResidueNumber = "-1000";  # to catch the first residue
		
		for $CurLine (0 .. $#AtomIndex) {
			# MODEL, TER and ENDMDL entries do not belong to any residue and have
			# to be ignored here
			if ($AtomIndex[$CurLine]->{Race} =~ m/(MODEL|TER|ENDMDL)/) { next }
			
			if (not defined $AtomIndex[$CurLine]{ResidueNumber}) {
				$ResidueNumber = "";
			}
			else {
				$ResidueNumber = $AtomIndex[$CurLine]{ResidueNumber};
			}
			
			if (not defined $AtomIndex[$CurLine]{InsResidue}) {
				$InsResidue = "";
			}
			else {
				$InsResidue = $AtomIndex[$CurLine]{InsResidue};
			}
			
			$CurResidueNumber = $ResidueNumber . $InsResidue;
			
			if ($OldResidueNumber ne $CurResidueNumber) {
				++$Residue;
				
				$OldResidueNumber = $CurResidueNumber;
				my $Atoms = []; # create a new anonymous array
				$ResidueIndex[$Residue]{Atoms} = $Atoms;
			}
			
			if ($AtomIndex[$CurLine]->{Race} ne "TER") {
				my %Atom = %{$AtomIndex[$CurLine]};
				push @{$ResidueIndex[$Residue]->{Atoms}}, \%Atom;
			}
		}
		
		return @ResidueIndex;
	}
	else {
		return @Content;
	}
} # of sub Get


sub Write { # writes a specified part of the PDB or the whole PDB, if no part is specified
	my $self = shift;
	my %args = @_;
	my ($Content, $AtomIndex, $OutFile, $BaseName, $Header);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	if (not defined $args{Header}) { $args{Header} = 1 } # by default header is included, MinHeader still possible
	if (not defined $args{Footer}) { $args{Footer} = 1 } # by default footer is included, MinFooter still possible
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (%args)) {
		$self->_ReportWarning ("Write", "ChainLabelsInvalid");
		return undef;
	}
	
	# forward all given parameters to Get and retrieve its return value
	($AtomIndex, $Content) = $self->Get (%args, GetReference => 1);
	$BaseName = $self->_GetBaseName ($self->{FileName});
	
	if (defined $args{FileName}) { $OutFile = $args{FileName}                      }
	                        else { $self->_ReportError ("NoFile", $args{FileName}) }
	
	open (OutFile, ">$OutFile") or $self->_ReportError ("IOError", $OutFile);
	print OutFile @{$Content};
	close OutFile;
} # of sub Write


sub WriteModels { # writes the single models and returns a list with the file names
	my $self = shift;
	my %args = @_;
	my (@Models, $Model, $BaseName, $Suffix, $CurFile, @FileList);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	if (defined $args{FileName}) {
		$BaseName = $args{FileName};
		$BaseName = $self->_GetBaseName ($BaseName);
	}
	else {
		$BaseName = $self->_GetBaseName ($self->{FileName})
	}
	
	if (defined $args{ModelSuffix}) { $Suffix = $args{ModelSuffix} } # if Suffix is defined, take the specified value
	                           else { $Suffix = $self->{ModelSuffix} } # else the default value
	
	@Models = $self->IdentifyModels;
	
	foreach $Model (@Models) {
		$CurFile = $BaseName . $Suffix . $Model . ".pdb";
		push @FileList, $CurFile;
		$args{FileName} = $CurFile;
		$args{Model} = $Model;
		$self->Write (%args); # all other paramters like header, footer, atom, type, ... are "forwarded"
	} # of foreach $Model (@Models)
	
	return @FileList; # return the list with all the produced files
} # of sub WriteModels


sub WriteChains { # writes the single chains and returns a list with the file names
	my $self = shift;
	my %args = @_;
	my ($Model, @Chains, $Chain, $ChainLabel, $FileName, $Suffix, $CurFile, @FileList);
	if (not defined $self->{Models}) { $self->Parse }
	
	if (defined $args{FileName}) {
		$FileName = $args{FileName};
		$FileName = $self->_GetBaseName ($FileName);
	}
	else {
		$FileName = $self->_GetBaseName ($self->{FileName})
	}

	if (defined $args{ChainSuffix}) { $Suffix = $args{ChainSuffix}   } # if Suffix is defined, take the specified value
	                           else { $Suffix = $self->{ChainSuffix} } # else take the default value
	
	$Model = $self->_GetModelOrDefault (\%args); # if Model is defined, take the specified value, otherwise take the default value
	if (not defined $Model) { return }
	
	@Chains = $self->IdentifyChains (Model => $Model);
	
	foreach $Chain (@Chains) {
		if ($self->{ChainLabelAsLetter} and ($self->ChainLabelsValid (Model => $Model)) ) {
			$ChainLabel = $self->GetChainLabel (Chain => $Chain);
		}
		else { $ChainLabel = $Chain }
		
		$CurFile = $FileName . $Suffix . $ChainLabel . ".pdb";
		push @FileList, $CurFile;
		$args{FileName} = $CurFile;
		$args{Chain} = $Chain;
		$self->Write (%args); # all other parameters like header, footer, atom, type, ... are "forwarded"
	} # of foreach $Chain (@Chains)
	
	return @FileList; # return the list with all the produced files
} # of sub WriteChains


sub WriteFASTA { # writes a FASTA file
	my $self = shift;
	my %args = @_;
	my ($BaseName, @FASTA);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	if (defined $args{FileName}) {
		$BaseName = $self->_GetBaseName ($args{BaseName});
	}
	else {
		$BaseName = $self->{BaseName};
	}
	
	@FASTA = $self->GetFASTA;
	
	if (not @FASTA) { return undef }
	
	open OUTPUT, ">$BaseName.fasta";
	print OUTPUT @FASTA;
	close OUTPUT;
	
	return @FASTA; # just to return anything
} # of sub WriteFASTA


#######################################################################################################################
## Get-Methods
#######################################################################################################################

# These methods are used to interconvert identifiers, that is to convert one *single* identifier.
# GetInternalIdentifier returns the external identifier, and vice versa, e.g. GetModel (ModelNumber => ...) returns
# the internal Model number.

sub GetModel { # returns the internal model number of an external ModelNumber
	my $self = shift;
	my %args = @_;
	my ($ModelNumber, $ModelIndex, $Model, $CurModelNumber);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$ModelNumber = $args{ModelNumber};
	
	if (not defined $ModelNumber) { return undef }
	
	$ModelIndex = $self->{Models}{ModelIndex};
	
	for $Model (0 .. $#{$ModelIndex}) {
		if ($ModelNumber eq $ModelIndex->[$Model]) { return $Model }
	}
	
	$self->_ReportWarning ("GetModel", "UnknownNumber", "ModelNumber $ModelNumber");
	return undef; # if the model has not been found
} # of sub GetModel


sub GetModelNumber { # returns the ModelNumber of an internal model number
	my $self = shift;
	my %args = @_;
	my $Model;
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$Model = $self->_GetModelOrDefault (\%args); # if Model is defined, take the specified value, otherwise take the default value
	
	if (not defined $self->{Models}{ModelIndex}[$Model]) {
		$self->_ReportWarning ("GetModelNumber", "UnknownModel", "Model $Model");
		return undef;
	}
	else {
		return $self->{Models}{ModelIndex}[$Model];
	}
} # of sub GetModelNumber


sub GetChain { # returns the internal chain number of a given ChainLabel
	my $self = shift;
	my %args = @_;
	my ($ChainLabel, $ChainIndex, $Model, @Chains, $Chain, $ChainCount, $CurLabel);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$ChainLabel = $args{ChainLabel};
	
	if (not defined $ChainLabel) { # if ChainLabel is not defined, return undef
		$self->_ReportWarning ("GetChain", "Misc", "ChainLabel must be defined for GetChain!");
		return undef;
	}
	
	$Model = $self->_GetModelOrDefault (\%args); # if Model is defined, take the specified value, otherwise take the default value
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (Model => $Model)) {
		$self->_ReportWarning ("GetChain", "ChainLabelsInvalid");
		return undef;
	}
	
	$ChainIndex = $self->{Models}{$Model}{ChainIndex};
	
	for $Chain (0 .. $#{$ChainIndex}) {
		if ($ChainLabel eq $ChainIndex->[$Chain]) { return $Chain }
	}
	
	# if the ChainLabel has not been found (there would have been a "return $Chain" otherwise)
	$self->_ReportWarning ("GetChain", "UnknownLabel", "ChainLabel $ChainLabel");
	return undef;
} # of sub GetChain


sub GetChainLabel { # returns the "real" chain ID of a chain
	my $self = shift;
	my %args = @_;
	my ($Model, $Chain, $ChainLabel);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$Model = $self->_GetModelOrDefault (\%args); # if Model is defined, take the specified value, otherwise take the default value
	if (not defined $Model) { return undef }
	
	$Chain = $self->_GetArgsInternal (\%args, "Chain");
	
	if (not defined $Chain) {
		if (not defined $args{Chain}) {
			$self->_ReportWarning ("GetChainLabel", "Misc", "Chain must be defined for GetChainLabel!");
		}
		else {
			$self->_ReportWarning ("GetChainLabel", "UnknownChain", "Chain $Chain");
		}
		
		return undef;
	}
	
	$ChainLabel = $self->{Models}{$Model}{$Chain}{ChainLabel};
	
	return $ChainLabel;
} # of sub GetChainLabel


sub GetResidue { # returns the internal residue number of a given ResidueNumber (to be correct: the FIRST matching one)
	my $self = shift;
	my %args = @_;
	
	my ($Model, $Chain, $Residue, $ResidueNumber, $ResidueIndex);
	my (@Chains, $ResidueCount, $InsResidue);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$Model = $self->_GetModelOrDefault (\%args);
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (Model => $Model)) {
		$self->_ReportWarning ("GetResidue", "ChainLabelsInvalid");
		return undef;
	}
	
	$Chain = $self->_GetArgsInternal (\%args, "Chain");
	
	# if Atoms and Residues have not been indexed in this model yet, do it now
	if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
	
	if (not defined $args{ResidueNumber}) {
		$self->_ReportWarning ("GetResidue", "Misc", "ResidueNumber needs to be defined!");
		return undef;
	}
	
	$ResidueNumber = $args{ResidueNumber};
	
	# if the ResidueNumber ends with a character, split it up into number and InsResidue identifier
	if ($ResidueNumber =~ m/\d+[A-Za-z]+/) {
		$InsResidue = $ResidueNumber;
		$InsResidue =~ s/^\d+//;    # remove all leading numbers
		$ResidueNumber =~ s/[A-Za-z]+$//; # remove all trailing letters
	}
	else {
		$InsResidue = " ";  # if no InsResidue letter is given, a blank is always the default
	}
	
	if (defined $Chain) {
		$ResidueIndex = $self->{Models}{$Model}{$Chain}{ResidueIndex};
		
		for $Residue (0 .. $#{$ResidueIndex}) {
			if ($ResidueNumber . $InsResidue eq 
				     $ResidueIndex->[$Residue]{ResidueNumber} . $ResidueIndex->[$Residue]{InsResidue}) {
					return $Residue;
			}
		}
	}
	else { # if no chain was given
		$ResidueCount = -1; # gets increased before return
		@Chains = $self->_GetNumericKeys ($self->{Models}{$Model});
		
		foreach $Chain (@Chains) {
			$ResidueIndex = $self->{Models}{$Model}{$Chain}{ResidueIndex};
			
			for $Residue (0 .. $#{$ResidueIndex}) {
				++$ResidueCount;
				
			if ($ResidueNumber . $InsResidue eq 
				     $ResidueIndex->[$Residue]{ResidueNumber} . $ResidueIndex->[$Residue]{InsResidue}) {
						return $ResidueCount;
				}
			}
		} # of foreach $Chain (@Chains)
	} # of else
	
	$self->_ReportWarning ("GetResidue", "UnknownResidue", "ResidueNumber $ResidueNumber");
	return undef; # if the residue has not been found
} # of sub GetResidue


sub GetResidueNumber { # returns the external residue number of a given Residue
	my $self = shift;
	my %args = @_;
	
	my ($Model, $Chain, $Residue, $ResidueNumber, $ResidueIndex);
	my (@Chains, $ResidueCount);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$Model = $self->_GetModelOrDefault (\%args);
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (Model => $Model)) {
		$self->_ReportWarning ("GetResidueNumber", "ChainLabelsInvalid");
		return undef;
	}
	
	$Chain = $self->_GetArgsInternal (\%args, "Chain");
	
	# if Atoms and Residues have not been indexed in this model yet, do it now
	if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
	
	if (not defined $args{Residue}) {
		$self->_ReportWarning ("GetResidueNumber", "Misc", "Residue needs to be defined!");
		return undef;
	}
	
	$Residue = $args{Residue};
	
	if (defined $Chain) {
		if ($Residue > $#{$self->{Models}{$Model}{$Chain}{ResidueIndex}}) {
			$self->_ReportWarning ("GetResidueNumber", "UnknownResidue", "Residue $Residue");
		}
		else {
			if ($args{GetReference}) {
				return $self->{Models}{$Model}{$Chain}{AtomIndex},
				       $self->{Models}{$Model}{$Chain}{ResidueIndex}[$Residue];
			}
			else {
				return $self->{Models}{$Model}{$Chain}{ResidueIndex}[$Residue]{ResidueNumber};
			}
		}
	}
	else { # if no chain was given
		$ResidueCount = 0;
		@Chains = $self->_GetNumericKeys ($self->{Models}{$Model});
		
		foreach $Chain (@Chains) {
			$ResidueIndex = $self->{Models}{$Model}{$Chain}{ResidueIndex};
			
			if ($Residue > $ResidueCount + $#{$ResidueIndex}) {
				$ResidueCount = $ResidueCount + $#{$ResidueIndex} + 1;
				next;
			}
			
			if ($args{GetReference}) {
				return $self->{Models}{$Model}{$Chain}{AtomIndex},
						 $ResidueIndex->[$Residue - $ResidueCount];
			}
			else {
				return $ResidueIndex->[$Residue - $ResidueCount]{ResidueNumber};
			}
		} # of foreach $Chain (@Chains)
	} # of else
	
	$self->_ReportWarning ("GetResidueNumber", "UnknownResidue", "Residue $Residue");
	return undef; # if the residue has not been found
} # of sub GetResidueNumber


sub GetResidueLabel { # returns the ResidueLabel of a given Residue or ResidueNumber
	my $self = shift;
	my %args = @_;
	
	my ($OneLetterCode, $Residue, $ResidueNumber, $Model, $Chain, @Chains);
	my ($ResidueCount, $ResidueIndex, $ResidueLabel);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$Model = $self->_GetModelOrDefault (\%args);
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (Model => $Model)) {
		$self->_ReportWarning ("GetResidueLabel", "ChainLabelsInvalid");
		return undef;
	}
	
	$OneLetterCode = $args{OneLetterCode};
	
	$Chain = $self->_GetArgsInternal (\%args, "Chain");
	
	if (defined $args{Residue}) {
		$ResidueIndex = $self->GetResidueNumber (%args, GetReference => 1);
		if (not defined $ResidueIndex) { return undef }
		
		$ResidueLabel = $ResidueIndex->{ResidueLabel};
		
		if ($OneLetterCode) {
			$ResidueLabel = $self->AminoAcidConvert ($ResidueLabel);
		}
		
		return $ResidueLabel;
	}
	elsif (defined $args{ResidueNumber}) {
		$ResidueIndex = $self->GetResidue (%args, GetReference => 1);
		if (not defined $ResidueIndex) { return undef }
		
		$ResidueLabel = $ResidueIndex->{ResidueLabel};
		
		if ($OneLetterCode) {
			$ResidueLabel = $self->AminoAcidConvert ($ResidueLabel);
		}
		
		return $ResidueLabel;
	}
	else {
		$self->_ReportWarning ("GetResidueLabel", "Misc", "Residue or ResidueNumber needs to be defined!");
		return undef;
	}
} # of sub GetResidueLabel


sub GetAtom { # returns the internal atom number of a given external atom number
	my $self = shift;
	my %args = @_;
	my ($Model, $Chain, $Residue, $Atom, $AtomNumber, $AtomCount);
	my (@Chains, $ResidueCount, $ResidueIndex, $AtomIndex, $CurNumber);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$Model = $self->_GetModelOrDefault (\%args);
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (Model => $Model)) {
		$self->_ReportWarning ("GetAtom", "ChainLabelsInvalid");
		return undef;
	}
	
	$Chain = $self->_GetArgsInternal (\%args, "Chain");
	$Residue = $self->_GetArgsInternal (\%args, "Residue");
	
	# if Atoms and Residues have not been indexed in this model yet, do it now
	if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
	
	if (not defined $args{AtomNumber}) {
		$self->_ReportWarning ("GetAtom", "Misc", "AtomNumber needs to be defined!");
		return undef;
	}
	
	$AtomNumber = $args{AtomNumber};
	
	if (defined $Residue) { # find an atom relative to a residue
		# That's admitedly a bit heavy now. Since we got $Residue, we need
		# GetResidueNumber (since $Residue is the input for it). With the
		# Switch GetReference it returns the index of that residue along with
		# the AtomIndex of that chain
		($AtomIndex, $ResidueIndex) = $self->GetResidueNumber (%args, GetReference => 1);
		
		$AtomCount = 0; # count the atoms in that very residue
		
		foreach $Atom ( @{$ResidueIndex->{Atoms}} ) {
			if ($AtomNumber eq $AtomIndex->[$Atom]{AtomNumber}) {
				if ($args{GetReference}) { return ($AtomIndex, $Atom) }
				                    else { return $AtomCount          }
			}
			else {
				++$AtomCount;
			}
		}
		
		$self->_ReportWarning ("GetAtom", "UnknownAtom", "AtomNumber $AtomNumber");
		return undef; # if the atom has not been found
	}
	elsif (defined $Chain) { # find an atom relative to a chain
		$AtomIndex = $self->{Models}{$Model}{$Chain}{AtomIndex};
		$AtomCount = 0;
		
		for $Atom (0 .. $#{$AtomIndex}) {
			# a TER does not count as an "Atom", so skip it
			if ($AtomIndex->[$Atom]{Race} eq "TER") { next }
			
			if ($AtomNumber eq $AtomIndex->[$Atom]{AtomNumber}) {
				if ($args{GetReference}) { return ($AtomIndex, $Atom) }
				                    else { return $AtomCount          }
			}
			else {
				++$AtomCount;
			}
		}
		
		$self->_ReportWarning ("GetAtom", "UnknownAtom", "AtomNumber $AtomNumber");
		return undef; # if the atom has not been found
	}
	else { # find an atom relative to a model
		@Chains = $self->_GetNumericKeys ($self->{Models}{$Model});
		$AtomCount = 0;
		
		foreach $Chain ( @Chains ) {
			$AtomIndex = $self->{Models}{$Model}{$Chain}{AtomIndex};
			
			for $Atom (0 .. $#{$AtomIndex}) {
				# a TER does not count as an "Atom", so skip it
				if ($AtomIndex->[$Atom]{Race} eq "TER") { next }
				
				if ($AtomNumber eq $AtomIndex->[$Atom]{AtomNumber}) {
					if ($args{GetReference}) { return ($AtomIndex, $Atom) }
					                    else { return $AtomCount          }
				}
				else {
					++$AtomCount;
				}
			}
		}
		
		$self->_ReportWarning ("GetAtom", "UnknownAtom", "AtomNumber $AtomNumber");
		return undef; # if the atom has not been found
	}
} # of sub GetAtom


sub GetAtomNumber { # returns the external atom number of a given internal atom number
	my $self = shift;
	my %args = @_;
	my ($Model, $Chain, $Residue, $Atom, $AtomNumber, $AtomCount);
	my (@Chains, $ResidueCount, $ResidueIndex, $AtomIndex, $CurNumber);
	my ($ResidueAtoms, $CurAtom);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$Model   = $self->_GetModelOrDefault (\%args);
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (Model => $Model)) {
		$self->_ReportWarning ("GetAtomNumber", "ChainLabelsInvalid");
		return undef;
	}
	
	if (not defined $args{Atom}) {
		$self->_ReportWarning ("GetAtomNumber", "Misc", "Atom needs to be defined for GetAtomNumber!");
		return undef;
	}
	
	# if Atoms and Residues have not been indexed in this model yet, do it now
	if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
	
	$Chain   = $self->_GetArgsInternal (\%args, "Chain");
	$Residue = $self->_GetArgsInternal (\%args, "Residue");
	
	$Atom = $args{Atom};
	
	if (defined $Residue) { # find an atom relative to a residue
		($AtomIndex, $ResidueIndex) = $self->GetResidueNumber (%args, GetReference => 1);
		
		# read all atom indeces in the requested residue
		$ResidueAtoms = $ResidueIndex->{Atoms};
		
		if ($#{$ResidueAtoms} < $Atom) {
			$self->_ReportWarning ("GetAtomNumber", "UnknownAtom", "Atom $Atom");
			return undef; # if the atom has not been found
		}
		
		# read the index of the requested atom in that residue
		$CurAtom = $ResidueAtoms->[$Atom];
		
		$AtomNumber = $AtomIndex->[$CurAtom]{AtomNumber};
		
		if ($args{GetReference}) { return ($AtomIndex, $CurAtom) }
		                    else { return $AtomNumber            }
	}
	elsif (defined $Chain) { # find an atom relative to a chain
		$AtomIndex = $self->{Models}{$Model}{$Chain}{AtomIndex};
		
		if ($#{$AtomIndex} < $Atom) {
			$self->_ReportWarning ("GetAtomNumber", "UnknownAtom", "Atom $Atom");
			return undef; # if the atom has not been found
		}
		
		$AtomNumber = $AtomIndex->[$Atom]{AtomNumber};
		
		if ($args{GetReference}) { return ($AtomIndex, $CurAtom) }
		                    else { return $AtomNumber              }
	}
	else { # find an atom relative to a model
		@Chains = $self->_GetNumericKeys ($self->{Models}{$Model});
		$AtomCount = 0;
		
		foreach $Chain ( @Chains ) {
			$AtomIndex = $self->{Models}{$Model}{$Chain}{AtomIndex};
			
			if ($Atom > $AtomCount + $#{$AtomIndex}) {
				$AtomCount = $AtomCount + $#{$AtomIndex};
				next;
			}
			
			$AtomNumber = $AtomIndex->[$Atom - $AtomCount]{AtomNumber};
			
			if ($args{GetReference}) {
				return ($AtomIndex, $Atom - $AtomCount);
			}
			else {
				return $AtomNumber;
			}
		}
		
		$self->_ReportWarning ("GetAtomNumber", "UnknownAtom", "Atom $Atom");
		return undef; # if the atom has not been found
	}
	
} # of sub GetAtomNumber


#######################################################################################################################
## Identify-Methods
#######################################################################################################################

# These Methods return a list of the requested identifier (Model, ChainLabel, AtomType, ...) in the given domain
# (the whole model, a single chain or just one residue). Note that only the domain can be narrowed, there is no
# filtering (e.g. for an AtomType) done.
# Implementation of filtering would make it possible to retrieve all AtomNumbers of all C-Alpha atoms, but first of
# all it would complicate things in the code and second of all are such queries also possible via Get (retrieve all
# C-Alphas and read their numbers).
# That said, the definition of the identify methods is that they return an array with all the respective
# identifiers, which then can be fed into Get. Example:
# @Chains = IndentifyChains;
# foreach $Chain (@Chains) { @Chain = Get (Chain = $Chain) };


sub IdentifyModels { # returns an array with the hash keys to all models
	my $self = shift;
	
	if (not defined $self->{Models}) { $self->Parse }
	
	return $self->_GetNumericKeys ($self->{Models});
} # of sub IdentifyModels


sub IdentifyModelNumbers { # returns an array with the hash keys to all models
	my $self = shift;
	
	if (not defined $self->{Models}) { $self->Parse }
	
	return @{$self->{Models}{ModelIndex}};
} # of sub IdentifyModelNumbers


sub IdentifyChains { # returns an array with the chain indeces
	my $self = shift;
	my %args = @_;
	my ($Model, @Chains);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$Model = $self->_GetModelOrDefault (\%args); # if Model is defined, take the specified value, otherwise take the default value
	if (not defined $Model) { return undef }
	
	@Chains = $self->_GetNumericKeys ($self->{Models}{$Model});
	
	return @Chains;
} # of sub IdentifyChains


sub IdentifyChainLabels { # returns an array with the letters of the chains (rather use IdentifyChains for accessing the chains!))
	my $self = shift;
	my %args = @_;
	my ($Model, @ChainLabels);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	# if Model is defined, take the specified value, otherwise take the default value
	$Model = $self->_GetModelOrDefault (\%args);
	if (not defined $Model) { return undef }
	
	@ChainLabels = @{$self->{Models}{$Model}{ChainIndex}};
	
	return @ChainLabels;
} # of sub IdentifyChainLabels


sub IdentifyResidues { # returns an array with the residue indeces
	my $self = shift;
	my %args = @_;
	my ($Model, $Chain, @Chains, @Residues, $ResidueCount);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$Model = $self->_GetModelOrDefault (\%args);
	if (not defined $Model) { return undef }
	
	# if Atoms and Residues have not been indexed in this model yet, do it now
	if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (Model => $Model)) {
		$self->_ReportWarning ("IdentifyResidues", "ChainLabelsInvalid");
		return undef;
	}
	
	$Chain = $self->_GetArgsInternal (\%args, "Chain");
	$ResidueCount = 0;
	
	if (defined $args{Chain}) {
		push @Chains, $Chain;
	}
	else {
		@Chains = $self->_GetNumericKeys ($self->{Models}{$Model});
	}
	
	foreach $Chain (@Chains) {
		$ResidueCount = $ResidueCount + scalar @{$self->{Models}{$Model}{$Chain}{ResidueIndex}}
	}
	
	@Residues = (0 .. $ResidueCount - 1);
	
	return @Residues;
} # of IdentifyResidues


sub IdentifyResidueLabels { # returns an array with the residue labels (the actual number in the PDB)
	my $self = shift;
	my %args = @_;
	my ($Model, @Chains, $Chain, $Residue, @ResidueLabels, $ResidueIndex, $CurResidue, $OneLetterCode);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$Model = $self->_GetModelOrDefault (\%args);
	
	if (not defined $Model) { return undef }
	
	# if Atoms and Residues have not been indexed in this model yet, do it now
	if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (Model => $Model)) {
		$self->_ReportWarning ("IdentifyResidueLabels", "ChainLabelsInvalid");
		return undef;
	}
	
	$Chain = $self->_GetArgsInternal (\%args, "Chain");
	$OneLetterCode = $args{OneLetterCode};
	
	if (defined $args{Chain}) {
		push @Chains, $Chain;
	}
	else {
		@Chains = $self->_GetNumericKeys ($self->{Models}{$Model});
	}
	
	foreach $Chain (@Chains) {
		$ResidueIndex = $self->{Models}{$Model}{$Chain}{ResidueIndex};
		
		for $Residue (0 .. $#{$ResidueIndex}) {
			$CurResidue = $ResidueIndex->[$Residue]{ResidueLabel};
			
			if ($OneLetterCode) { $CurResidue = $self->AminoAcidConvert ($CurResidue) }
			push @ResidueLabels, $CurResidue;
		}
	}
	
	return @ResidueLabels;
} # of sub IdentifyResidueLabels


sub IdentifyResidueNumbers { # returns an array with the external residue numbers
	my $self = shift;
	my %args = @_;
	my ($Model, @Chains, $Chain, $Residue, @ResidueNumbers, $ResidueIndex);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$Model = $self->_GetModelOrDefault (\%args);
	if (not defined $Model) { return undef }
	
	# if Atoms and Residues have not been indexed in this model yet, do it now
	if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (Model => $Model)) {
		$self->_ReportWarning ("IdentifyResidueNumbers", "ChainLabelsInvalid");
		return undef;
	}
	
	$Chain = $self->_GetArgsInternal (\%args, "Chain");
	
	if (defined $args{Chain}) {
		push @Chains, $Chain;
	}
	else {
		@Chains = $self->_GetNumericKeys ($self->{Models}{$Model});
	}
	
	foreach $Chain (@Chains) {
		$ResidueIndex = $self->{Models}{$Model}{$Chain}{ResidueIndex};
		
		for $Residue (0 .. $#{$ResidueIndex}) {
			if ($ResidueIndex->[$Residue]{InsResidue} ne " ") {
				push @ResidueNumbers, $ResidueIndex->[$Residue]{ResidueNumber} . $ResidueIndex->[$Residue]{InsResidue};
			}
			else {
				push @ResidueNumbers, $ResidueIndex->[$Residue]{ResidueNumber};
			}
		}
	}
	
	return @ResidueNumbers;
} # of sub IdentifyResidueNumbers


sub IdentifyAtoms { # returns an array with the atom indeces
	my $self = shift;
	my %args = @_;
	my ($Model, $Chain, @Chains, $Residue, $Atom, @Atoms);
	my ($ChainHash, $AtomIndex, $AtomCount);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$Model = $self->_GetModelOrDefault (\%args);
	if (not defined $Model) { return undef }
	
	# if Atoms and Residues have not been indexed in this model yet, do it now
	if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (Model => $Model)) {
		$self->_ReportWarning ("IdentifyAtoms", "ChainLabelsInvalid");
		return undef;
	}
	
	$Chain   = $self->_GetArgsInternal (\%args, "Chain");
	$Residue = $self->_GetArgsInternal (\%args, "Residue");
	
	if (defined $Residue) {
		$ChainHash = $self->_RetrieveResidue (\%args, 1);
		$Residue = $args{Residue}; # update $Residue, might have been changed by RetrieveResidue
		$AtomCount = $#{$ChainHash->{ResidueIndex}[$Residue]{Atoms}};
		return (0 .. $AtomCount);
	}
	
	if (defined $args{Chain}) {
		push @Chains, $Chain;
	}
	else {
		@Chains = $self->_GetNumericKeys ($self->{Models}{$Model});
	}
	
	$AtomCount = -1; # gets increased before push
	
	foreach $Chain (@Chains) {
		$AtomIndex = $self->{Models}{$Model}{$Chain}{AtomIndex};
		
		foreach $Atom (@{$AtomIndex}) {
			if ($Atom->{Race} ne "TER") {
				++$AtomCount;
				push @Atoms, $AtomCount;
			}
		}
	}
	
	return @Atoms;
} # of sub IdentifyAtoms


sub IdentifyAtomNumbers { # returns an array with the atom numbers
	my $self = shift;
	my %args = @_;
	my ($Model, $Chain, @Chains, $Residue, $Atom);
	my ($ChainHash, $AtomIndex, @AtomNumbers);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$Model = $self->_GetModelOrDefault (\%args);
	if (not defined $Model) { return undef }
	
	# if Atoms and Residues have not been indexed in this model yet, do it now
	if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (Model => $Model)) {
		$self->_ReportWarning ("IndentifyAtomNumbers", "ChainLabelsInvalid");
		return undef;
	}
	
	$Chain   = $self->_GetArgsInternal (\%args, "Chain");
	$Residue = $self->_GetArgsInternal (\%args, "Residue");
	
	if (defined $Residue) {
		$ChainHash = $self->_RetrieveResidue (\%args, 1);
		$Residue = $args{Residue}; # update $Residue, might have been changed by RetrieveResidue
		
		foreach $Atom (@{$ChainHash->{ResidueIndex}[$Residue]{Atoms}}) {
			push @AtomNumbers, $ChainHash->{AtomIndex}[$Atom]{AtomNumber};
		}
	}
	else {
		if (defined $args{Chain}) {
			push @Chains, $Chain;
		}
		else {
			@Chains = $self->_GetNumericKeys ($self->{Models}{$Model});
		}
		
		foreach $Chain (@Chains) {
			$AtomIndex = $self->{Models}{$Model}{$Chain}{AtomIndex};
			
			foreach $Atom (@{$AtomIndex}) {
				if ($Atom->{Race} ne "TER") {
					push @AtomNumbers, $Atom->{AtomNumber};
				}
			}
		}
	}
	
	return @AtomNumbers;
} # of sub IdentifyAtomNumbers


sub IdentifyAtomTypes { # returns an array with the atom types
	my $self = shift;
	my %args = @_;
	my ($Model, $Chain, @Chains, $Residue, $Atom, @AtomTypes);
	my ($ChainHash, $AtomIndex, $AtomCount);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$Model = $self->_GetModelOrDefault (\%args);
	if (not defined $Model) { return undef }
	
	# if Atoms and Residues have not been indexed in this model yet, do it now
	if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (Model => $Model)) {
		$self->_ReportWarning ("IdentifyAtomTypes", "ChainLabelsInvalid");
		return undef;
	}
	
	$Chain   = $self->_GetArgsInternal (\%args, "Chain");
	$Residue = $self->_GetArgsInternal (\%args, "Residue");
	
	if (defined $Residue) {
		$ChainHash = $self->_RetrieveResidue (\%args, 1);
		$Residue = $args{Residue}; # update $Residue, might have been changed by RetrieveResidue
		
		foreach $Atom (@{$ChainHash->{ResidueIndex}[$Residue]{Atoms}}) {
			push @AtomTypes, $ChainHash->{AtomIndex}[$Atom]{AtomType};
		}
	}
	else {
		if (defined $args{Chain}) {
			push @Chains, $Chain;
		}
		else {
			@Chains = $self->_GetNumericKeys ($self->{Models}{$Model});
		}
		
		$AtomCount = -1; # gets increased before push
		foreach $Chain (@Chains) {
			$AtomIndex = $self->{Models}{$Model}{$Chain}{AtomIndex};
			
			foreach $Atom (@{$AtomIndex}) {
				++$AtomCount;
				
				if ($Atom->{Race} !~ m/^TER/) {
					push @AtomTypes, $Atom->{AtomType};
				}
			}
		}
	}
	
	return @AtomTypes;
} # of sub IdentifyAtomTypes


sub IdentifyElements { # returns an array with the elements
	my $self = shift;
	my %args = @_;
	my ($Model, $Chain, @Chains, $Residue, $Atom, $Element, @Elements);
	my ($ChainHash, $AtomIndex, $AtomCount);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$Model = $self->_GetModelOrDefault (\%args);
	if (not defined $Model) { return undef }
	
	# if Atoms and Residues have not been indexed in this model yet, do it now
	if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (Model => $Model)) {
		$self->_ReportWarning ("IdentifyElements", "ChainLabelsInvalid");
		return undef;
	}
	
	$Chain   = $self->_GetArgsInternal (\%args, "Chain");
	$Residue = $self->_GetArgsInternal (\%args, "Residue");
	
	if (defined $Residue) {
		$ChainHash = $self->_RetrieveResidue (\%args, 1);
		$Residue = $args{Residue}; # update $Residue, might have been changed by RetrieveResidue
		
		foreach $Atom ( @{$ChainHash->{ResidueIndex}[$Residue]{Atoms}} ) {
			$Element = $self->_GetSingleElement ($ChainHash->{AtomIndex}[$Atom]);
			push @Elements, $Element; 
		}
	}
	else {
		if (defined $args{Chain}) {
			push @Chains, $Chain;
		}
		else {
			@Chains = $self->_GetNumericKeys ($self->{Models}{$Model});
		}
		
		$AtomCount = -1; # gets increased before push
		foreach $Chain (@Chains) {
			$AtomIndex = $self->{Models}{$Model}{$Chain}{AtomIndex};
			
			foreach $Atom (@{$AtomIndex}) {
				++$AtomCount;
				
				if ($Atom->{Race} !~ m/^TER/) {
					$Element = $self->_GetSingleElement ($Atom);
					push @Elements, $Element;
				}
			}
		}
	}
	
#		else  {$self->_ReportWarning ("IdentifyElements", "UnknownType", $CurType)}
	
	return @Elements;
} # of sub IdentifyElements


#######################################################################################################################
## Several methods to fetch information from the PDB hash
#######################################################################################################################

sub GetAtomType { # returns the AtomType of a given Atom or AtomNumber
	my $self = shift;
	my %args = @_;
	my ($AtomType, $AtomIndex, $Atom);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	if ( (not defined $args{Atom}) and (not defined $args{AtomNumber}) ) {
		$self->_ReportWarning ("GetAtomType", "Misc", "Atom or AtomNumber needs to be defined!");
		return undef;
	}
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (%args)) {
		$self->_ReportWarning ("GetAtomType", "ChainLabelsInvalid");
		return undef;
	}
	
	if (defined $args{Atom}) {
		($AtomIndex, $Atom) = $self->GetAtomNumber (%args, GetReference => 1);
		
		if ($args{GetReference}) {
			return ($AtomIndex, $Atom);
		}
		else {
	 		return $AtomIndex->[$Atom]{AtomType};
		}
	}
	elsif (defined $args{AtomNumber}) {
		($AtomIndex, $Atom) = $self->GetAtom (%args, GetReference => 1);
		
		if ($args{GetReference}) {
			return ($AtomIndex, $Atom);
		}
		else {
			return $AtomIndex->[$Atom]{AtomType};
		}
	}
} # of sub GetAtomType


sub GetElement { # returns the element of a given Atom or AtomNumber
	my $self = shift;
	my %args = @_;
	my (@Models, $Model, @Chains, $Chain, $Residue, $AtomIndex);
	my ($LineCount, $CurAtom, $Atom, $AtomNumber, $AtomType, $Key, $Element, $Content);
		
	if (not defined $self->{Models}) { $self->Parse }
	
	$Model = $self->_GetModelOrDefault (\%args);
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (Model => $Model)) {
		$self->_ReportWarning ("GetElement", "ChainLabelsInvalid");
		return undef;
	}
	
	$Chain   = $self->_GetArgsInternal (\%args, "Chain");
	$Residue = $self->_GetArgsInternal (\%args, "Residue");
	$Atom    = $self->_GetArgsInternal (\%args, "Atom");
	
	# if Atoms and Residues have not been indexed in this model yet, do it now
	if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
	
	# if no parameters are given, the whole Model is processed and the results saved in the main hash
	if (not defined $args{Chain} and not defined $args{Residue} and not defined $args{Atom}) {
		@Chains = $self->IdentifyChains (Model => $Model);
		
		foreach $Chain (@Chains) {
			$AtomIndex = $self->{Models}{$Model}{$Chain}{AtomIndex};
			
			for $Atom (0 .. $#{$AtomIndex} ) {
				$Element = $self->_GetSingleElement ($AtomIndex->[$Atom]);
				$AtomIndex->[$Atom]{Element} = $Element;
			}
		} # of foreach $Chain (@Chain)
	}
	else { # if parameters are given, the request is processed
		($AtomIndex, $Atom) = $self->GetAtomType (%args, GetReference => 1);
		$AtomType = $AtomIndex->[$Atom]{AtomType};
		
		if (not defined $AtomType) { return undef }
		else {
			if ($AtomIndex->[$Atom]{AtomType} eq "CA") {
				
				if ($AtomIndex->[$Atom]{Race} eq "ATOM") { return "C" }
				                                    else { return "CA" }
			}
			else {
				foreach $Key ( keys %ELEMENTFILTERS ) {
					if ($AtomType =~ m/^$ELEMENTFILTERS{$Key}$/) {
						return $Key;
					}
				}
			}
		}
		
		$self->_ReportWarning ("GetElement", "Misc", "The Element $AtomType has not been recognized.");
		return undef;
	}
} # of sub GetElement


sub GetMass { # returns the mass of a given set of atoms
	my $self = shift;
	my %args = @_;
	my (@Models, $Model, @Chains, $Chain, $Residue, $AtomIndex);
	my ($AtomType, $Element, $Mass, $CurMass, $Content, $Atom);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$Model   = $self->_GetModelOrDefault (\%args);
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (Model => $Model)) {
		$self->_ReportWarning ("GetMass", "ChainLabelsInvalid");
		return undef;
	}
	
	$Chain   = $self->_GetArgsInternal (\%args, "Chain");
	$Residue = $self->_GetArgsInternal (\%args, "Residue");
	$Atom    = $self->_GetArgsInternal (\%args, "Atom");
	
	# if Atoms and Residues have not been indexed in this model yet, do it now
	if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
	
	# if no parameters are given, the whole model is processed and the results saved in the main hash
	if (not defined $Chain and not defined $Residue and not defined $Atom) {
		@Chains = $self->IdentifyChains (Model => $Model);
		$Mass = 0;
		
		foreach $Chain (@Chains) {
			$AtomIndex = $self->{Models}{$Model}{$Chain}{AtomIndex};
			
			for $Atom (0 .. $#{$AtomIndex} ) {
				if (not defined $AtomIndex->[$Atom]{Element}) {
					$Element = $self->_GetSingleElement ($AtomIndex->[$Atom]);
					$AtomIndex->[$Atom]{Element} = $Element;
				}
				
				if ( ($AtomIndex->[$Atom]{Race} eq "ATOM") or ($AtomIndex->[$Atom]{Race} eq "HETATM") ) {
					if (defined $Element) {
						$CurMass = $ELEMENTMASS{$Element};
						$AtomIndex->[$Atom]{Mass} = $CurMass;
						$Mass = $Mass + $CurMass;
					}
				}
			}
		} # of foreach $Chain (@Chain)
		
		$Mass = sprintf ("%.5f", $Mass);
		return $Mass;
	}
	else { # if parameters are given, the request is processed
		($AtomIndex, $Content) = $self->Get (%args, GetReference => 1);
		$Mass = 0;
		
		for $Atom (0 .. $#{$AtomIndex}) {
			if ( ($AtomIndex->[$Atom]{Race} eq "ATOM") or ($AtomIndex->[$Atom]{Race} eq "HETATM") ) {
				$Element = $self->_GetSingleElement ($AtomIndex->[$Atom]);
				
				if (not defined $Element) { return 0 } # warning is issued by _GetSingleElement
				
				$Mass = $Mass + $ELEMENTMASS{$Element};
			}
		} # of for $Atom (0 .. $#{$AtomIndex})
		
		$Mass = sprintf ("%.5f", $Mass);
		return $Mass;
	}
} # of sub GetMass


sub GetCenterOfMass { # returns the center of mass of a given domain (chain, residue, ...)
	my $self = shift;
	my %args = @_;
	my ($Atom, $AtomIndex, $Content, %CenterOfMass, $CurMass, $Element);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (%args)) {
		$self->_ReportWarning ("GetCenterOfMass", "ChainLabelsInvalid");
		return undef;
	}
	
	($AtomIndex, $Content) = $self->Get (%args, GetReference => 1);
	
	%CenterOfMass = ( Mass => 0, x => 0, y => 0, z => 0 );
	
	for $Atom (0 .. $#{$AtomIndex}) {
		if ( ($AtomIndex->[$Atom]{Race} eq "ATOM") or ($AtomIndex->[$Atom]{Race} eq "HETATM") ) {
			$Element = $self->_GetSingleElement ($AtomIndex->[$Atom]);
			
			if (not defined $Element) { return undef } # warning is issued by _GetSingleElement
			
			$CurMass = $ELEMENTMASS{$Element};
			if (not defined $CurMass) { return }
			
			$CenterOfMass{Mass} = $CenterOfMass{Mass} + $CurMass;
			$CenterOfMass{x} = $CenterOfMass{x} + $AtomIndex->[$Atom]{x} * $CurMass;
			$CenterOfMass{y} = $CenterOfMass{y} + $AtomIndex->[$Atom]{y} * $CurMass;
			$CenterOfMass{z} = $CenterOfMass{z} + $AtomIndex->[$Atom]{z} * $CurMass;
		}
	}
	
	if ($CenterOfMass{Mass} > 0) { # if at least one atom was found in the residue
		$CenterOfMass{x} = $CenterOfMass{x} / $CenterOfMass{Mass};
		$CenterOfMass{y} = $CenterOfMass{y} / $CenterOfMass{Mass};
		$CenterOfMass{z} = $CenterOfMass{z} / $CenterOfMass{Mass};
	}
	
	return %CenterOfMass;
} # of sub GetCenterOfMass


sub GetCentreOfMass { # to catch the british speakers
	my $self = shift;
	my %CentreOfMass = $self->GetCenterOfMass (@_);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	return %CentreOfMass;
}


sub GetCoordinates { # returns the coordinates of the requested atoms
	my $self = shift;
	my %args = @_;
	my (@Coordinates, $AtomIndex, $Content, $xCoord, $yCoord, $zCoord, $LineCount);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (%args)) {
		$self->_ReportWarning ("GetCoordinates", "ChainLabelsInvalid");
		return undef;
	}
		
	($AtomIndex, $Content) = $self->Get (@_, GetReference => 1);
	
	for $LineCount (0 .. $#{$AtomIndex}) {
		if ($AtomIndex->[$LineCount]{Race} !~ m/^($self->{AllAtoms})/) { next }
		
		$Coordinates[$LineCount] = {
			x => $AtomIndex->[$LineCount]{x},
			y => $AtomIndex->[$LineCount]{y},
			z => $AtomIndex->[$LineCount]{z}
		};
	} # of for $LineCount (0 .. $#{$AtomIndex})
	
	return @Coordinates;
} # of sub GetCoordinates


sub GetAngles { # returns the phi and psi angles of a given residue
	my $self = shift;
	my %args = @_;
	my (@Models, $Model, @Chains, $Chain, $Residue, @Residues, @AllAngles);
	my (@AllResidues, $ReturnResidue, $Content, $Distance);
	my ($AtomC, $AtomN, $AtomCA, $Temp, $CurChain);
	my ($AtomIndex, $Line, $CurResidueNumber, $ResidueCount);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (%args)) {
		$self->_ReportWarning ("GetAngles", "ChainLabelsInvalid");
		return undef;
	}
	
	$Model   = $self->_GetModelOrDefault (\%args);
	$Chain   = $self->_GetArgsInternal (\%args, "Chain");
	$Residue = $self->_GetArgsInternal (\%args, "Residue");
		
	if (defined $Residue) { # if a single residue has to be processed
		$ReturnResidue = -1; # this is the index of the residue which is to be returned
		
		$args{Residue} = $Residue;
		($Temp, $Content) = $self->Get (%args, GetReference => 1);
		
		if (defined $Temp) {
			push @{$AtomIndex}, @{$Temp};
			++$ReturnResidue;
		}
		
		# Turn off warnings, in case the residue doesn't exist
		$self->{ReportNoWarnings} = 1;
		
		if ($Residue > 0) {
			$args{Residue} = $Residue - 1;
			($Temp, $Content) = $self->Get (%args, GetReference => 1);
			
			if (defined $Temp) {
				unshift @{$AtomIndex}, (@{$Temp}); # add it before the requested residue
				++$ReturnResidue;
			}
		}
		
		$args{Residue} = $Residue + 1;
		($Temp, $Content) = $self->Get (%args, GetReference => 1);
		
		if (defined $Temp) {
			push @{$AtomIndex}, @{$Temp};
		}
		
		# Turn on the warnings again
		$self->{ReportNoWarnings} = 0;
	} # of if (defined $Residue)
	elsif (defined $Chain) { # if a whole chain has to be processed
		($AtomIndex, $Content) = $self->Get (%args, GetReference => 1);
	}
	else { # if nothing was specified, the whole model will be calculated and the data saved in the main hash
		my @ModelAngles;
		
		@Chains = $self->IdentifyChains (Model => $Model);
		
		foreach $Chain (@Chains) {
			my @Angles = $self->GetAngles (Model => $Model, Chain => $Chain);
			push @ModelAngles, @Angles;
			
			for $Residue (0 .. $#Angles) { # for each line
				$self->{Models}{$Model}{$Chain}{ResidueIndex}[$Residue]{Phi} = $Angles[$Residue]{Phi};
				$self->{Models}{$Model}{$Chain}{ResidueIndex}[$Residue]{Psi} = $Angles[$Residue]{Psi};
			}	
		}
		
		return @ModelAngles;
	}
	
	if ($#{$AtomIndex} < 1) {
		$self->_ReportWarning ("GetAngles", "Misc", "At least two residues are needed to calculate angles!");
		return undef;
	}

	$CurResidueNumber = $AtomIndex->[0]{ResidueNumber} . $AtomIndex->[0]{InsResidue};
	
	for $Line (0 .. $#{$AtomIndex} + 1) {
		
		# after the last residue, add it to the array
		if ($Line > $#{$AtomIndex}) {
			push @AllResidues, {N => $AtomN, CA => $AtomCA, C => $AtomC};
			last;
		}
		
		# if the residue hast changed, add the last one to the array
		if ( (defined $CurResidueNumber)
				and
			  ($AtomIndex->[$Line]{ResidueNumber}.$AtomIndex->[$Line]{InsResidue}  ne $CurResidueNumber) ) {
			$CurResidueNumber = $AtomIndex->[$Line]{ResidueNumber} . $AtomIndex->[$Line]{InsResidue};
			push @AllResidues, {N => $AtomN, CA => $AtomCA, C => $AtomC};
		}
		
		if (defined $AtomIndex->[$Line]{AtomType}) {
			
			if ($AtomIndex->[$Line]{AtomType} eq "N") {
				$AtomN = $AtomIndex->[$Line];
			}
			
			if ($AtomIndex->[$Line]{AtomType} eq "CA") {
				$AtomCA = $AtomIndex->[$Line];
			}
			
			if ($AtomIndex->[$Line]{AtomType} eq "C") {
				$AtomC = $AtomIndex->[$Line];
			}
		}
	} # of for $Line (0 .. $#{$AtomIndex)}
	
	for $Residue (0 .. $#AllResidues) {
		my %Angles;
		
		if ($Residue > 0) { # the first residue of a chain has no Phi angle
			# check for a gap in the chain
			$Distance = $self->_Distance ($AllResidues[$Residue-1]->{C}, $AllResidues[$Residue]->{N});
			
			if ( (not defined $Distance) or ($Distance > 2) ) {
				$Angles{Phi} = "360.000";
			}
			else {
				if ( (not defined $AllResidues[$Residue-1]->{C})  or
				     (not defined $AllResidues[$Residue]->{N})    or
				     (not defined $AllResidues[$Residue]->{CA})   or
				     (not defined $AllResidues[$Residue]->{C})      ) {
					$Angles{Phi} = "364.000";
				}
				else {
					$Angles{Phi} = $self->_DihedralAngle ($AllResidues[$Residue-1]->{C},
					                                      $AllResidues[$Residue]->{N},
					                                      $AllResidues[$Residue]->{CA},
					                                      $AllResidues[$Residue]->{C});
				}
			}
		}
		else { $Angles{Phi} = "360.000" }
		
		if ($Residue < $#AllResidues) { # the last residue of a chain has no Psi angle
			# check for a gap in the chain
			$Distance = $self->_Distance ($AllResidues[$Residue]->{C}, $AllResidues[$Residue+1]->{N});
			
			if ( (not defined $Distance) or ($Distance > 2) ) {
				$Angles{Psi} = "360.000";
			}
			else {
				if ( (not defined $AllResidues[$Residue]->{N})   or
				     (not defined $AllResidues[$Residue]->{CA})  or
				     (not defined $AllResidues[$Residue]->{C})   or
				     (not defined $AllResidues[$Residue+1]->{N}) ) {
					$Angles{Psi} = "360.000";
				}
				else {
					$Angles{Psi} = $self->_DihedralAngle ($AllResidues[$Residue]->{N},
					                                      $AllResidues[$Residue]->{CA},
					                                      $AllResidues[$Residue]->{C},
					                                      $AllResidues[$Residue+1]->{N});
				}
			}
		}
		else { $Angles{Psi} = "360.000" }
		
		push @AllAngles, \%Angles;
	}
	
	if (defined $ReturnResidue) { # if only one single residue is to be returned
		return ($AllAngles[$ReturnResidue]);
	}
	else { # return all
		return @AllAngles;
	}
} # of sub GetAngles


sub GetHeader { # returns the header
	my $self = shift;
	my %args = @_;
	my ($Content, $AtomIndex, $Header, $HeaderIndex);
	my ($LineCount, $Race, $Rest);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	if ($args{Content} and $args{AtomIndex}) {
		$Content = $args{Content};
		$AtomIndex = $args{AtomIndex};
		
		# if the header has not been indexed yet
		if (not defined $self->{HeaderIndex}) {
			$Header = $self->{Header};
			
			for $LineCount (0 .. $#{$Header}) {
				$Race = substr ($Header->[$LineCount], 0, 6);
				$Race =~ s/\s//g;
				
				# at the moment, just the first six characters and the rest is splitted
				$Rest = substr ($Header->[$LineCount], 6);
				
				push @{$HeaderIndex}, { Race => $Race, Rest => $Rest };
			}
			
			$self->{HeaderIndex} = $HeaderIndex;
		}
		
		if (not $args{MinHeader}) {
			# add Header and HeaderIndex to the current dataset in the arguments
			if (defined $self->{Header}) {
				unshift @{$Content}, @{$self->{Header}};
			}
			
			if (defined $self->{HeaderIndex}) {
				unshift @{$AtomIndex}, @{$self->{HeaderIndex}};
			}
			
			if ( $self->{HeaderRemark} ) {
				$self->_HeaderRemark (Content => $Content, AtomIndex => $AtomIndex);
			}
		}
		else { # if only the minimal header was requested
			# two seperate arrays are required, unshifting it line by line would sort
			# it in the wrong order
			my (@MinHeader, @MinHeaderIndex);
			
			for $LineCount (0 .. $#{$self->{HeaderIndex}}) {
				if ($self->{HeaderIndex}[$LineCount]{Race} =~ m/^($MinHeaderLines)/) {
					push @MinHeaderIndex, $self->{HeaderIndex}[$LineCount];
					push @MinHeader, $self->{Header}[$LineCount];
				}
			}
			
			unshift @{$Content}, @MinHeader;
			unshift @{$AtomIndex}, @MinHeaderIndex;
		}
	}
	else {
		if (not $args{MinHeader}) {
			return @{$self->{Header}};
		}
		else {
			my @MinHeader;
			for $LineCount (0 .. $#{$self->{Header}}) {
				if ($self->{Header}[$LineCount] =~ m/^($MinHeaderLines)/) {
					push @MinHeader, $self->{Header}[$LineCount];
				}
			}
			
			return @MinHeader;
		}
	}
} # of sub GetHeader


sub GetFooter { # returns the footer
	my $self = shift;
	my %args = @_;
	my ($Content, $AtomIndex, $Footer, $FooterIndex);
	my ($LineCount, $Race, $Rest);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	if ($args{Content} and $args{AtomIndex}) {
		$Content = $args{Content};
		$AtomIndex = $args{AtomIndex};
		
		# if the footer has not been indexed yet
		if (defined $self->{Footer} and not defined $self->{FooterIndex}) {
			$Footer = $self->{Footer};
			
			for $LineCount (0 .. $#{$Footer}) {
				$Race = substr ($Footer->[$LineCount], 0, 6);
				$Race =~ s/\s//g;
				
				# at the moment, just the first six characters and the rest is splitted
				$Rest = substr ($Footer->[$LineCount], 6);
				
				push @{$FooterIndex}, { Race => $Race, Rest => $Rest };
			}
			
			$self->{FooterIndex} = $FooterIndex;
		}
		
		if (not $args{MinFooter}) {
			if (defined $self->{Footer} and defined $self->{FooterIndex}) {
				# add Footer and FooterIndex to the current dataset in the arguments
				push @{$Content}, @{$self->{Footer}};
				push @{$AtomIndex}, @{$self->{FooterIndex}};
			}
		}
		else {
			for $LineCount (0 .. $#{$self->{FooterIndex}}) {
				if ($self->{FooterIndex}[$LineCount]{Race} =~ m/^($MinFooterLines)/) {
					push @{$Content}, $self->{Footer}[$LineCount];
					push @{$AtomIndex}, $self->{FooterIndex}[$LineCount];
				}
			}
		}
	}
	else {
		if (not $args{MinFooter}) {
			return @{$self->{Footer}};
		}
		else {
			my @MinFooter;
			for $LineCount (0 .. $#{$self->{Footer}}) {
				if ($self->{Footer}[$LineCount] =~ m/^($MinFooterLines)/) {
					push @MinFooter, $self->{Footer}[$LineCount];
				}
			}
			
			return @MinFooter;
		}
	}
} # of sub GetFooter


sub GetMinHeader { # return only a minimal header
	my $self = shift;
	
	if (not defined $self->{Models}) { $self->Parse }
	
	return $self->GetHeader (@_, MinHeader => 1);
} # of sub GetMinHeader


sub GetMinFooter { # return only a minimal footer
	my $self = shift;
	
	if (not defined $self->{Models}) { $self->Parse }
	
	return $self->GetFooter (@_, MinFooter => 1);
} # of sub GetMinFooter


sub GetSection { # returns lines starting with a certain word
	my $self = shift;
	my $Pattern = shift;
	my (%args, @Content);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$args{Pattern} = $Pattern;
	
	@Content = $self->GetHeader;
	push @{$args{Content}}, @Content;
	
	$self->_RetrieveAllModels (\%args);
	
	@Content = $self->GetFooter;
	push @{$args{Content}}, @Content;
	
	$self->_FilterLines (\%args);
	
	return @{$args{Content}};
} # of sub GetSection


sub GetResolution { # returns the resolution
	my $self = shift;
	my ($Resolution, @Header, %args, $Line);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	@Header = $self->GetHeader;
	
	$args{Content} = \@Header;
	$args{Pattern} = "REMARK   2";
	
	$self->_FilterLines (\%args);
	
	for $Line (0 .. $#{$args{Content}}) {
		if ($args{Content}->[$Line] =~ m/RESOLUTION/) {
			
			if ($args{Content}->[$Line] =~ m/NOT APPLICABLE/) {
				return undef; # returns the element of a given atom number
			}
			
			$Resolution = substr ($args{Content}->[$Line], 22, 5);
			$Resolution =~ s/\s//g; # remove all blanks
			return $Resolution;
		}
	}
} # of sub GetResolution


sub GetFASTA { # returns an array with FASTA format line
	my $self = shift;
	my %args = @_;
	my ($Model, @Chains, $Chain, @Residues, $Residue);
	my (@FASTA, $Line);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$Model = $self->_GetModelOrDefault (\%args);
	if (not defined $Model) { return undef }
	
	# Check whether an external chain label was specified and block this at once if the
	# chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (Model => $Model)) {
		$self->_ReportWarning ("GetFASTA", "ChainLabelsInvalid");
		return undef;
	}
	
	$Chain = $self->_GetArgsInternal (\%args, "Chain");
	
	if ($Chain) { push @Chains, $Chain }
	       else { @Chains = $self->IdentifyChains (Model => $Model) }
	
	foreach $Chain (@Chains) {
		if ($self->CountAtoms (Model => $Model, Chain => $Chain, Race => "ATOM") == 0) {
			# skip chains containing only HETATMs
			next;
		}
		
		push @FASTA, ">$self->{BaseName}:_|PDBID|CHAIN|SEQUENCE\n";
		
		@Residues = $self->IdentifyResidueLabels (Model => $Model, Chain => $Chain, OneLetterCode => 1);
		$Line = "";
		
		foreach $Residue ( @Residues ) {
			# to ignore for example water and other HETATMs
			if (not defined $Residue) { next }
			
			$Line = $Line . $Residue;
			
			if (length $Line == 80) {
				$Line = $Line . "\n";
				push @FASTA, $Line;
				$Line = "";
			}
		}
	}
	
	# add the remaining data to the array
	if ($Line) { push @FASTA, $Line; }
	
	return @FASTA
} # of sub GetFASTA


sub RemoveInsertedResidues { # removes inserted residues in case they are superpositions
	my $self = shift;
	my %args = @_;
	
	my (@Models, $Model, @Chains, $Chain, $Residue, $Atom, $Remove, $SkipIt);
	my ($Content, $ResidueIndex, $AtomIndex, $Atoms, $Warning, $Intensive);
	my ($CurLine, %Res1, %Res2, %Distance, $Type, $RepeatSearch);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	# if Intensive is true, only residues with a non-blank InsResidue tag and
	# the one before them is checked
	if (defined $args{Intensive}) { $Intensive = $args{Intensive} }
	                         else { $Intensive = 0                }
	
	# @Models = $self->IdentifyModels;
	
	$Model = $self->_GetModelOrDefault (\%args);
	
	# if Atoms and Residues have not been indexed in this model yet, do it now
	if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
	
	@Chains = $self->IdentifyChains (Model => $Model);
	
	# if two superimposed residues are found, the search has to be repeated to catch
	# e.g. three successive residues
	$RepeatSearch = 0;
	
	foreach $Chain (@Chains) {
		$Content      = $self->{Models}{$Model}{$Chain}{Content};
		$AtomIndex    = $self->{Models}{$Model}{$Chain}{AtomIndex};
		$ResidueIndex = $self->{Models}{$Model}{$Chain}{ResidueIndex};
		
		$Residue = 0;
		
		# run over all residues and determine distance to the next residue in the chain
		while ($Residue < $#{$ResidueIndex} - 1) {
			$SkipIt = 0;
			
			if ( (not $Intensive) and
					($ResidueIndex->[$Residue]{InsResidue} eq " ") and
					($ResidueIndex->[$Residue+1]{InsResidue} eq " ") )
			{
				++$Residue;
				next;
			}
			
			TYPE: foreach $Type ( "C", "N", "O", "CA" ) {
				# Determine C,N,O,CA atoms of the current and the next residue
				foreach $Atom ( @{$ResidueIndex->[$Residue]{Atoms}} ) {
					if (not defined $AtomIndex->[$Atom]{AtomType}) {
						$SkipIt = 1;
						last TYPE;
					}
					
					if ($AtomIndex->[$Atom]{AtomType} eq $Type) {
						$Res1{$Type} = $AtomIndex->[$Atom];
					}
				}
				
				if (not defined $Res1{$Type}) {
					$SkipIt = 1;
					last;
				}
				
				# Determine C,N,O,CA atoms of the next residue
				foreach $Atom ( @{$ResidueIndex->[$Residue+1]{Atoms}} ) {
					if (not defined $Atom or not defined $AtomIndex->[$Atom]{AtomType}) {
						last;
					}
					
					if ($AtomIndex->[$Atom]{AtomType} eq $Type) {
						$Res2{$Type} = $AtomIndex->[$Atom];
					}
				}
				
				if (not defined $Res2{$Type}) {
					$SkipIt = 1;
					last;
				}
			}
			
			if ($SkipIt) {
				$self->_ReportWarning ("RemoveInsertedResidues", "Misc", "Skipped residue $Residue in chain $Chain");
				++$Residue;
				next;
			}
			
			# calculate the distances of the atoms
			foreach $Type ( "C", "N", "O", "CA" ) {
				$Distance{$Type} = sqrt (
				                          ($Res1{$Type}->{x} - $Res2{$Type}->{x})**2
				                        + ($Res1{$Type}->{y} - $Res2{$Type}->{y})**2
				                        + ($Res1{$Type}->{z} - $Res2{$Type}->{z})**2 );
			}
			
			# if all distances a smaller than one Angstrom, remove one of them
			if ($Distance{C} < 1 and $Distance{N} < 1 and
					$Distance{O} < 1 and $Distance{CA} < 1) {
				
				$Warning = "Residues " .
				           $ResidueIndex->[$Residue]{ResidueNumber} . $ResidueIndex->[$Residue]{InsResidue} . " and " .
				           $ResidueIndex->[$Residue+1]{ResidueNumber} . $ResidueIndex->[$Residue+1]{InsResidue} . " are superimposed. ";
				
				$RepeatSearch = 1;
				$Remove = -1;
				
				# check whether the second residue has an InsResidue Tag
				if ($ResidueIndex->[$Residue+1]{InsResidue} ne " ") {
					$Warning = $Warning . " The second one had an InsResidue mark and was removed.";
					$Remove = $Residue + 1;
				}
				# check whether the first residue has an InsResidue Tag
				elsif ($ResidueIndex->[$Residue]{InsResidue} ne " ") {
					$Warning = $Warning . " The first one had an InsResidue mark and was removed.";
					$Remove = $Residue;
				}
				else { # if none of the residues has an InsResidue tag
					$Warning = $Warning . " No InsResidue mark was found, thus the second one was removed.";
					$Remove = $Residue + 1;
				}
				
				if ($Remove > -1) { # Remove this residue
					$Atoms = $ResidueIndex->[$Remove]{Atoms};
					splice (@{$AtomIndex}, $Atoms->[0], scalar @{$Atoms});
					splice (@{$Content},   $Atoms->[0], scalar @{$Atoms});
					splice (@{$ResidueIndex}, $Remove, 1);
					
					$self->_ReportWarning ("RemoveInsertedResidues", "Misc", $Warning);
				}
			}
			
			++$Residue;
		} # of for $Residue (1 .. $#{$ResidueIndex})
	} # of foreach $Chain (@Chains)
} # of sub RemoveInsertedResidues


sub RemoveAtomLocations { # removes alternative atom locations
	my $self = shift;
	my %args = @_;
	
	my (@Models, $Model, @Chains, $Chain, $ChainHash);
	my ($AtomLocations);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	if (not defined $args{AtomLocations}) { $args{AtomLocations} = $self->{AtomLocations} }
	                                 else { $args{AtomLocations} = $args{AtomLocations}   }
	
	$Model = $self->_GetModelOrDefault (\%args);
	
	# if Atoms and Residues have not been indexed in this model yet, do it now
	if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
	
	@Chains = $self->IdentifyChains (Model => $Model);
	
	foreach $Chain (@Chains) {
		my (@Content, @AtomIndex);
		
		$args{Content} = \@Content;
		$args{AtomIndex} = \@AtomIndex;
		
		$ChainHash = $self->{Models}{$Model}{$Chain};
		
		$self->_FilterAllLines (\%args, $ChainHash);
		
		$self->{Models}{$Model}{$Chain}{Content} = \@Content;
		$self->{Models}{$Model}{$Chain}{AtomIndex} = \@AtomIndex;
	}
	
	$self->_IndexResidues ($Model);
} # of sub RemoveAtomLocations


#######################################################################################################################
## Count Methods
#######################################################################################################################

sub CountModels { # returns the number of models
	my $self = shift;
	my @Models = $self->IdentifyModels;
	
	if (not defined $self->{Models}) { $self->Parse }
	
	return scalar @Models;
} # of sub CountModels


sub CountChains { # returns the number of chains in a model
	my $self = shift;
	my %args = @_;
	
	if (not defined $self->{Models}) { $self->Parse }
	
	my $Model = $self->_GetModelOrDefault (\%args); # if Model is defined, take the specified value, otherwise take the default value
	if (not defined $Model) { return undef }
	
	my @Chains = $self->IdentifyChains (Model => $Model);
	return scalar @Chains; 
} # of sub CountChains


sub CountResidues { # returns the number of residues
	my $self = shift;
	my %args = @_;
	
	if (not defined $self->{Models}) { $self->Parse }
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (%args)) {
		$self->_ReportWarning ("CountResidues", "ChainLabelsInvalid");
		return undef;
	}
	
	my @Residues = $self->Get (%args, ResidueIndex => 1);
	
	if (not $Residues[0]) { return 0                }
	                 else { return scalar @Residues }
} # of sub CountResidues


sub CountAtoms { # counts all ATOM and HETATM lines
	my $self = shift;
	my %args = @_;
	my ($Pattern, $AtomIndex, $Content, $AtomCount, $Line);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	if ($args{Pattern}) { $Pattern = $args{Pattern}  }
	               else { $Pattern = "($self->{AllAtoms})" }
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	if (defined $args{ChainLabel} and not $self->ChainLabelsValid (%args)) {
		$self->_ReportWarning ("CountAtoms", "ChainLabelsInvalid");
		return undef;
	}
	
	($AtomIndex, $Content) = $self->Get (%args, GetReference => 1);
	
	$AtomCount = 0;
	
	for $Line (0 .. $#{$AtomIndex}) {
		if ($AtomIndex->[$Line]{Race} =~ m/$Pattern/) {
			++$AtomCount;
		}
	}
	
	return $AtomCount;
} # of sub CountAtoms


sub CountATOMs { # counts all ATOM lines
	my $self = shift;
	my %args = @_;
	
	if (not defined $self->{Models}) { $self->Parse }
	
	unless ($args{Pattern}) { $args{Pattern} = "ATOM" }
	return $self->CountAtoms (%args);
} # of sub CountATOMs


sub CountHETATMs { # counts all HETATM lines
	my $self = shift;
	my %args = @_;
	
	if (not defined $self->{Models}) { $self->Parse }
	
	unless ($args{Pattern}) { $args{Pattern} = "HETATM" }
	return $self->CountAtoms (%args);
} # of sub CountHETATMs


#######################################################################################################################
## BEGIN Renumber Methods
#######################################################################################################################

sub RenumberModels { # renumbers the models starting at 1 or the given value
	my $self = shift;
	my %args = @_;
	
	my ($Content, $AtomIndex, $CurModel, $CurLine, @Models, $Model);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	if (defined $args{ModelStart}) { $CurModel = $args{ModelStart} } # if a start value is given, use it
	                          else { $CurModel = 1 } # otherwise start counting at 0
	
	if (defined $args{Content} and defined $args{AtomIndex}) {
		# if Content and AtomIndex are defined, the routine has been called by Get
		
		$Content = $args{Content};
		$AtomIndex = $args{AtomIndex};
		
		for $CurLine (0 .. $#{$AtomIndex}) {
			if ($AtomIndex->[$CurLine]{Race} eq "MODEL") {
				$AtomIndex->[$CurLine]{ModelNumber} = $CurModel;
				
				substr ($Content->[$CurLine], 10, 4) = sprintf ("%4d", $CurModel);
				
				++$CurModel;
			}
		}
	}
	else { # renumber the whole PDB data
		@Models = $self->IdentifyModels;
		
		foreach $Model (@Models) {
			if (not defined $self->{Models}{$Model}{MODEL}) {
				$self->{Models}{$Model}{MODEL} = " " x 80;
			}
			
			$CurLine = \$self->{Models}{$Model}{MODEL};
			
			$self->{Models}{ModelIndex}[$Model] = $CurModel;
			substr (${$CurLine}, 10, 4) = sprintf ("%4d", $CurModel); # put in the current model number
			++$CurModel;
		}
	}
} # of sub RenumberModels


sub RenumberChains { # renumbers the chains starting at A or the given letter
	my $self = shift;
	my %args = @_;
	my (@Models, $Model, @Chains, $Chain, $AtomIndex, $Content);
	my ($ChainStart, $ChainLabel, $CurChainLabel, $OldChainLabel, $CurLine);
	my ($ModelString, $CurModel, $NewChain, $Pattern);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	if (defined $args{ChainStart}) { $ChainStart = $args{ChainStart} } # if a start letter is given, use it
	                          else { $ChainStart = 'A' } # otherwise start counting at 'A'
	
	if (defined $args{Model}) { push @Models, $args{Model} }
	                     else { @Models = $self->IdentifyModels }
	
	if (length $ChainStart > 1) {
		$self->_ReportWarning ("RenumberChains", "Misc", "Invalid chain ID given to start renumbering. A chain ID must be only one letter.");
		return undef;
	}
	
	$Pattern = "^($self->{AllAtoms}|TER)";
	
	if (defined $args{Content} and defined $args{AtomIndex}) {
		# if Content and AtomIndex are defined, the routine has been called by Get
		
		$Content = $args{Content};
		$AtomIndex = $args{AtomIndex};
		$CurChainLabel = $ChainStart;
		$NewChain = 1;
		
		for $CurLine (0 .. $#{$AtomIndex}) {
			if ( ($AtomIndex->[$CurLine]{Race} eq "MODEL") and ($CurChainLabel ne $ChainStart) ) {
				# if a new model begins and it is not the first one
				$CurChainLabel = $ChainStart; # reset the ChainLabel
				$NewChain = 1;
				next;
			}
			
			if ($AtomIndex->[$CurLine]{Race} =~ m/$Pattern/) {
				if ( (not $NewChain) and ($OldChainLabel ne $AtomIndex->[$CurLine]{ChainLabel}) ) { ++$CurChainLabel }
				
				$OldChainLabel = $AtomIndex->[$CurLine]{ChainLabel};
				$NewChain = 0;
				
				$AtomIndex->[$CurLine]{ChainLabel} = $CurChainLabel;
				substr ($Content->[$CurLine], 21, 1) = $CurChainLabel; # put in the current chain ID
				
				if ($AtomIndex->[$CurLine]{Race} eq "TER") {
					++$CurChainLabel;
					$NewChain = 1;;
					if ($CurChainLabel eq "AA") { $CurChainLabel = "a" }
					if ($CurChainLabel eq "aa") { $CurChainLabel = "0" }
					if ($CurChainLabel eq "10") { $CurChainLabel = "A" }
					next;
				}
			}
		}
	}
	else {
		$Model = $self->_GetModelOrDefault (\%args);
		
		# if Atoms and Residues have not been indexed in this model yet, do it now
		if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
		
		$CurChainLabel = $ChainStart;
		
		@Chains = $self->IdentifyChains (Model => $Model);
		
		foreach $Chain (@Chains) {
			$Content   = $self->{Models}{$Model}{$Chain}{Content};
			$AtomIndex = $self->{Models}{$Model}{$Chain}{AtomIndex};
			
			$self->{Models}{$Model}{$Chain}{ChainLabel} = $CurChainLabel;
			$self->{Models}{$Model}{ChainIndex}[$Chain] = $CurChainLabel;
			
			for $CurLine (0 .. $#{$AtomIndex}) {
				$AtomIndex->[$CurLine]{ChainLabel} = $CurChainLabel;   # put in the current chain ID
				substr ($Content->[$CurLine], 21, 1) = $CurChainLabel; # put in the current chain ID
			}
			
			++$CurChainLabel;
			if ($CurChainLabel eq "AA") { $CurChainLabel = "a" }
			if ($CurChainLabel eq "aa") { $CurChainLabel = "0" }
			if ($CurChainLabel eq "10") { $CurChainLabel = "A" }
		}
	}
	
	# check the chain labels again to correct the ChainLabelValid key if possible now
	$self->_IndexChains;
} # of sub RenumberChains


sub RenumberResidues { # renumbers the residues starting at 1 or the given number
	my $self = shift;
	my %args = @_;
	my (@Models, $Model, @Chains, $Chain, $AtomIndex, $ResidueIndex, $Content, $ChainLabel);
	my ($CurResidueNumber, $OldResidueNumber, $ResidueStart, $Residue, $CurLine, $Pattern);
	my ($CurModel, $FieldLength, $OldInsResidue, $KeepInsertions, $Atoms, $Atom, $NumberString);
	
	if (not defined $self->{Models}) { $self->Parse }

	if (defined $args{ResidueStart}) { $ResidueStart = $args{ResidueStart} } # if a start number is given, use it
	                            else { $ResidueStart = 1 } # otherwise start counting at 0
	
	if (defined $args{Model}) { push @Models, $args{Model} }
	                     else { @Models = $self->IdentifyModels }
	
	if (defined $args{KeepInsertions}) { $KeepInsertions = $args{KeepInsertions} }
	                              else { $KeepInsertions = 1 }
	
	if ($KeepInsertions) { $FieldLength = 4 } # the insertion codes are ignored, the same number persists with the code
	                else { $FieldLength = 5 } # the insertion codes are removed and all residues have different numbers
	
	$Pattern = "$self->{AllAtoms}|TER";
	
	if (defined $args{Content} and defined $args{AtomIndex}) {
		# if Content and AtomIndex are defined, the routine has been called by Get
		
		$Content = $args{Content};
		$AtomIndex = $args{AtomIndex};
		$CurResidueNumber = $ResidueStart - 1; # it gets increased at the first residue
		$OldResidueNumber = "-1000";
		$OldInsResidue = "-1000";
		
		for $CurLine (0 .. $#{$AtomIndex}) {
			if ($AtomIndex->[$CurLine]{Race} eq "MODEL") {
				$CurResidueNumber = $ResidueStart - 1; # reset the ResidueNumber
				$OldResidueNumber = "-1000";
				$OldInsResidue = "-1000";
				next;
			}
			
			if ($AtomIndex->[$CurLine]{Race} =~ m/$Pattern/) {
				if ($KeepInsertions) {
					if ( (defined $AtomIndex->[$CurLine]{ResidueNumber}) and
						  ($AtomIndex->[$CurLine]{ResidueNumber} != $OldResidueNumber) ) {
						++$CurResidueNumber;
						$OldResidueNumber = $AtomIndex->[$CurLine]{ResidueNumber};
					}
				}
				else { # if the insertion codes are to be removed
					if ( ( ($AtomIndex->[$CurLine]{ResidueNumber} != $OldResidueNumber)    # if the residue number is different
					        or                                                             # or
					        ( ($AtomIndex->[$CurLine]{ResidueNumber} == $OldResidueNumber) # the residuenumber is equal
					           and                                                         # but
					          ($AtomIndex->[$CurLine]{InsResidue} ne $OldInsResidue)       # the insertion code differs
					        )
					     )
					   )
					{
						++$CurResidueNumber;                                        # regard this residue as a new one
						$OldInsResidue = $AtomIndex->[$CurLine]{InsResidue};        # save the old marker for the next comparison
						$AtomIndex->[$CurLine]{InsResidue} = " ";                   # remove the old InsResidue mark
						$OldResidueNumber = $AtomIndex->[$CurLine]{ResidueNumber};
					}
				}
				
				$AtomIndex->[$CurLine]{ResidueNumber} = $CurResidueNumber;
				
				$NumberString = sprintf ("%4d", $CurResidueNumber); # make it 5 characters long
				if (not $KeepInsertions) {$NumberString = $NumberString . " " }
				
				substr ($Content->[$CurLine], 22, $FieldLength) = $NumberString; # put in the current ResidueNumber
			}
		}
	}
	else { # renumber the whole PDB data
		$Model = $self->_GetModelOrDefault (\%args);
		
		# if Atoms and Residues have not been indexed in this model yet, do it now
		if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
		
		$CurResidueNumber = $ResidueStart - 1; # reset the ResidueNumber
		$OldResidueNumber = "-1000";
		$OldInsResidue = "-1000";
		
		@Chains = $self->IdentifyChains (Model => $Model);
			
		foreach $Chain (@Chains) {
			$Content      = $self->{Models}{$Model}{$Chain}{Content};
			$AtomIndex    = $self->{Models}{$Model}{$Chain}{AtomIndex};
			$ResidueIndex = $self->{Models}{$Model}{$Chain}{ResidueIndex};
			
			foreach $Residue (@{$ResidueIndex}) {
				# check whether CurResidueNumber needs to be increased
				if ($KeepInsertions) {
					
					if ($Residue->{ResidueNumber} != $OldResidueNumber) {
						++$CurResidueNumber;
						$OldResidueNumber = $Residue->{ResidueNumber};
						$OldInsResidue = $Residue->{InsResidue};
						
						$Residue->{ResidueNumber} = $CurResidueNumber;
					}
					else { # if the ResidueNumber stayed the same, it is an inserted residue
						$Residue->{ResidueNumber} = $CurResidueNumber;
					}
				}
				else { # if the insertion codes are to be removed
					if ( ( ($Residue->{ResidueNumber} != $OldResidueNumber)    # if the residue number is different
							or                                                             # or
							( ($Residue->{ResidueNumber} == $OldResidueNumber) # the residuenumber is equal
								and                                                         # but
								($Residue->{InsResidue} ne $OldInsResidue)       # the insertion code differs
							)
						)
						)
					{
						++$CurResidueNumber;
						$OldResidueNumber = $Residue->{ResidueNumber};
						$OldInsResidue = $Residue->{InsResidue};
						
						$Residue->{ResidueNumber} = $CurResidueNumber;
						$Residue->{InsResidue} = " ";
					}
				}
				
				# read the reference to the array with the atom indeces
				$Atoms = $Residue->{Atoms};
				
				$NumberString = $CurResidueNumber;
				$NumberString = sprintf ("%4d", $CurResidueNumber); # make it 5 characters long
				if (not $KeepInsertions) { $NumberString = $NumberString . " " }
				
				# renumber all atom lines of this residue
				foreach $Atom (@{$Atoms}) { # $Atom becomes the index of the atom in the big hash
					$AtomIndex->[$Atom]{ResidueNumber} = $CurResidueNumber;
					substr ($Content->[$Atom], 22, $FieldLength) = $NumberString; # put in the current ResidueNumber
				}
			} # of for $Residue (0 .. $#{$ResidueIndex})
			
		}
	}
} # of sub RenumberResidues


sub RenumberAtoms { # renumbers the atoms starting at 1 or the given number
	my $self = shift;
	my %args = @_;
	my (@Models, $Model, @Chains, $Chain, $AtomIndex, $Content, $IgnoreTER);
	my ($AtomStart, $ChainLabel, $CurAtomNumber, $CurLine, $Pattern);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	if (defined $args{AtomStart}) { $AtomStart = $args{AtomStart} } # if a start letter is given, use it
	                         else { $AtomStart = 1 } # otherwise start counting at 'A'
	
	if (defined $args{IgnoreTER}) { $IgnoreTER = $args{IgnoreTER} }
	                         else { $IgnoreTER = 0 }
	
	$Model = $self->_GetModelOrDefault (\%args);
	$Pattern = "^($self->{AllAtoms}|TER)";
	
	if (defined $args{Content} and defined $args{AtomIndex}) {
		# if Content and AtomIndex are defined, the routine has been called by Get
		
		$Content = $args{Content};
		$AtomIndex = $args{AtomIndex};
		$CurAtomNumber = $AtomStart;
		
		for $CurLine (0 .. $#{$AtomIndex}) {
			
			if ($AtomIndex->[$CurLine]{Race} =~ "MODEL") {
				# if a new model begins
				$CurAtomNumber = $AtomStart - 1; # reset the AtomNumber
				next;
			}
			
			if ($AtomIndex->[$CurLine]{Race} =~ m/$Pattern/) {
				# if $IgnoreTER is true and a TER is found, $CurAtomNumber is NOT increased
				if ( $IgnoreTER and ($AtomIndex->[$CurLine]{Race} eq "TER") ) {
					# do nothing
				}
				elsif ($AtomIndex->[$CurLine]{Race} !~ m/^($self->{ANISIG})/) {
					++$CurAtomNumber;
				}
				
				$AtomIndex->[$CurLine]{AtomNumber} = $CurAtomNumber;
				substr ($Content->[$CurLine], 6, 5) = sprintf ("%5d",$CurAtomNumber); # put in the current AtomNumber
			}
			
		} # of while ($CurLine < $#{$AtomIndex})
	}
	else {
		# if Atoms and Residues have not been indexed in this model yet, do it now
		if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
		
		$CurAtomNumber = $AtomStart - 1;
		
		@Chains = $self->IdentifyChains (Model => $Model);
		
		foreach $Chain (@Chains) {
			$Content   = $self->{Models}{$Model}{$Chain}{Content};
			$AtomIndex = $self->{Models}{$Model}{$Chain}{AtomIndex};
			
			for $CurLine (0 .. $#{$AtomIndex}) {
				if ($AtomIndex->[$CurLine]{Race} =~ "MODEL") {
					# if a new model begins
					$CurAtomNumber = $AtomStart - 1; # reset the AtomNumber
					++$CurLine;
					next;
				}
				
				if ($AtomIndex->[$CurLine]{Race} =~ m/$Pattern/) {
					# if $IgnoreTER is true and a TER is found, $CurAtomNumber is NOT increased
					if ( $IgnoreTER and ($AtomIndex->[$CurLine]{Race} eq "TER") ) {
						# do nothing
					}
					elsif ($AtomIndex->[$CurLine]{Race} !~ m/^($self->{ANISIG})/) {
						++$CurAtomNumber;
					}
					
					$AtomIndex->[$CurLine]{AtomNumber} = $CurAtomNumber;
					substr ($Content->[$CurLine], 6, 5) = sprintf ("%5d", $CurAtomNumber); # put in the current AtomNumber
				}
			} # of for $CurLine (0 .. $#{$AtomIndex})
		}
	}
} # of sub RenumberAtoms


#######################################################################################################################
## END Renumber Methods
#######################################################################################################################


#######################################################################################################################
## Some other public methods
#######################################################################################################################

sub SetChainLabel { # changes the chain ID in all lines to the given value
	my $self = shift;
	my %args = @_;
	my ($AtomIndex, $Content, $ChainLabel, $LineCount);
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$Content = $args{Content};
	$AtomIndex = $args{AtomIndex};
	$ChainLabel = $args{ChainLabel};
	$LineCount = 0;
	
	if ( (not $Content) or (not $ChainLabel) )  { return undef }
	
	for $LineCount (0 .. $#{$AtomIndex}) { # loop over the whole file
		if ($AtomIndex->[$LineCount]{Race} =~ m/^($self->{AllAtoms}|TER)/) { # search for certain lines
			substr ($Content->[$LineCount], 21, 1) = $ChainLabel;
			$AtomIndex->[$LineCount]{ChainLabel} = $ChainLabel;
		}
	}
} # of sub SetChainLabel


sub AminoAcidConvert { # converts a 3-letter-code to a 1-letter-abbreviation
	my $self = shift;
	my $Code = shift;
	
	if (not defined $self->{Models}) { $self->Parse }
	
	$Code = uc ($Code);
	$Code =~ s/\s//g;  # remove all blanks
	
	if (length ($Code) == 3) {
		if    ($Code eq 'ALA') {$Code = 'A'} # Alanine
		elsif ($Code eq 'ASX') {$Code = 'B'} # Aspartic acid or Asparagine
		elsif ($Code eq 'CYS') {$Code = 'C'} # Cysteine
		elsif ($Code eq 'ASP') {$Code = 'D'} # Aspartic acid
		elsif ($Code eq 'GLU') {$Code = 'E'} # Glutamic acid
		elsif ($Code eq 'PHE') {$Code = 'F'} # Phenylalanine
		elsif ($Code eq 'GLY') {$Code = 'G'} # Glycine
		elsif ($Code eq 'HIS') {$Code = 'H'} # Histidine
		elsif ($Code eq 'ILE') {$Code = 'I'} # Isoleucine
		elsif ($Code eq 'LYS') {$Code = 'K'} # Lysine
		elsif ($Code eq 'LEU') {$Code = 'L'} # Leucine
		elsif ($Code eq 'MET') {$Code = 'M'} # Methionine
		elsif ($Code eq 'ASN') {$Code = 'N'} # Asparagine
		elsif ($Code eq 'PRO') {$Code = 'P'} # Proline
		elsif ($Code eq 'GLN') {$Code = 'Q'} # Glutamine
		elsif ($Code eq 'ARG') {$Code = 'R'} # Arginine
		elsif ($Code eq 'SER') {$Code = 'S'} # Serine
		elsif ($Code eq 'THR') {$Code = 'T'} # Threonine
		elsif ($Code eq 'SEC') {$Code = 'U'} # Selenocysteine
		elsif ($Code eq 'VAL') {$Code = 'V'} # Valine
		elsif ($Code eq 'TRP') {$Code = 'W'} # Tryptophan
		elsif ($Code eq 'UNK') {$Code = 'X'} # Unknown amino acid
		elsif ($Code eq 'TYR') {$Code = 'Y'} # Tyrosine
		elsif ($Code eq 'GLX') {$Code = 'Z'} # Glutamic acid or Glutamine
		else {
			$self->_ReportWarning ("AminoAcidConvert", "UnknownAminoAcid", $Code);
			return (undef);
		}
	} # of if (length ($Code) == 3)
	elsif (length ($Code) == 1) {
		if    ($Code eq 'A') {$Code = 'ALA'} # Alanine
		elsif ($Code eq 'B') {$Code = 'ASX'} # Aspartic acid or Asparagine
		elsif ($Code eq 'C') {$Code = 'CYS'} # Cysteine
		elsif ($Code eq 'D') {$Code = 'ASP'} # Aspartic acid
		elsif ($Code eq 'E') {$Code = 'GLU'} # Glutamic acid
		elsif ($Code eq 'F') {$Code = 'PHE'} # Phenylalanine
		elsif ($Code eq 'G') {$Code = 'GLY'} # Glycine
		elsif ($Code eq 'H') {$Code = 'HIS'} # Histidine
		elsif ($Code eq 'I') {$Code = 'ILE'} # Isoleucine
		elsif ($Code eq 'K') {$Code = 'LYS'} # Lysine
		elsif ($Code eq 'L') {$Code = 'LEU'} # Leucine
		elsif ($Code eq 'M') {$Code = 'MET'} # Methionine
		elsif ($Code eq 'N') {$Code = 'ASN'} # Asparagine
		elsif ($Code eq 'P') {$Code = 'PRO'} # Proline
		elsif ($Code eq 'Q') {$Code = 'GLN'} # Glutamine
		elsif ($Code eq 'R') {$Code = 'ARG'} # Arginine
		elsif ($Code eq 'S') {$Code = 'SER'} # Serine
		elsif ($Code eq 'T') {$Code = 'THR'} # Threonine
		elsif ($Code eq 'U') {$Code = 'SEC'} # Selenocysteine
		elsif ($Code eq 'V') {$Code = 'VAL'} # Valine
		elsif ($Code eq 'W') {$Code = 'TRP'} # Tryptophan
		elsif ($Code eq 'X') {$Code = 'UNK'} # Unknown amino acid
		elsif ($Code eq 'Y') {$Code = 'TYR'} # Tyrosine
		elsif ($Code eq 'Z') {$Code = 'GLX'} # Glutamic acid or Glutamine
		else {
			$self->_ReportWarning ("AminoAcidConvert", "UnknownAminoAcid", $Code);
			return undef;
		}
	} # of elsif (length ($Code) == 1)
	else {
		$self->_ReportWarning ("AminoAcidConvert", "InvalidAminoAcid");
		return undef;
	}
	return $Code;
} # of sub AminoAcidConvert


sub ChainLabelsValid { # returns true, if the chain IDs are valid and may be used for accessing the chains
	my $self = shift;
	my %args = @_;

	my $Model = $self->_GetModelOrDefault (\%args);
	return $self->{Models}{$Model}{ChainLabelsValid};
} # of sub ChainLabelsValid


sub FormatLine { # returns a PDB format line from a given atom hash
	my $self = shift;
	my %args = @_;
	
	my $Atom = $args{Atom};
	my ($TypForm, $Line);
	
	if ($Atom !~ m/^HASH\(0x[\da-z]+\)$/) {
		$self->_ReportWarning ("FormatLine", "Misc", "Parameter Atom needs to be an atom hash with PDB entries!");
		return undef;
	}
	
	if    ($Atom->{AtomType} =~ m/^\d/) {
		# if it starts with a number
		$TypForm = "%4s"; # no leading blanks
	}
	elsif ($Atom->{AtomType} =~ m/^CA$/ and $Atom->{Race} eq "HETATM") {
		# if it is a calcium atom
		$TypForm = "%4s"; # no leading blanks
	}
	else {
		$TypForm = " %-3s"; # one leading blank
	}
	
	$Line = sprintf "%-6s%5s $TypForm%1s%3s %1s%4s%s   %8.3f%8.3f%8.3f  %.2f %5s %13s\n",
					$Atom->{Race},
					$Atom->{AtomNumber},
					$Atom->{AtomType},
					$Atom->{AltLoc},
					$Atom->{ResidueLabel},
					$Atom->{ChainLabel},
					$Atom->{ResidueNumber},
					$Atom->{InsResidue},
					$Atom->{x},
					$Atom->{y},
					$Atom->{z},
					$Atom->{Occupancy},
					$Atom->{Temp},
					$Atom->{Rest};
	
	return $Line;
} # of sub FormatLine


#######################################################################################################################
## Accessor Methods for Global Variables
#######################################################################################################################

sub SetChainLabelAsLetter { # to set the global variable ChainLabelAsLetter
	my $self = shift;
	my $SetChainLabelAsLetter = shift;
	
	if (defined $SetChainLabelAsLetter) {
		$self->{ChainLabelAsLetter} = $SetChainLabelAsLetter
	}
	else {
		$self->{ChainLabelAsLetter} = $ChainLabelAsLetterDefault;
	}
} # of sub SetChainLabelAsLetter


sub SetChainSuffix { # to set the global variable ChainSuffix
	my $self = shift;
	my $SetChainSuffix = shift;
	
	if (defined $SetChainSuffix) {
		$self->{ChainSuffix} = $SetChainSuffix;
	}
	else {
		$self->{ChainSuffix} = $ChainSuffixDefault;
	}
} # of sub SetChainSuffix


sub SetModelSuffix { # to set the global variable ModelSuffix
	my $self = shift;
	my $SetModelSuffix = shift;
	
	if (defined $SetModelSuffix) {
		$self->{ModelSuffix} = $SetModelSuffix;
	}
	else {
		$self->{ModelSuffix} = $ModelSuffixDefault;
	}
} # of sub SetModelSuffix


sub SetHeaderRemark { # to set the global variable HeaderRemark
	my $self = shift;
	my $SetHeaderRemark = shift;
	
	if (defined $SetHeaderRemark) {
		$self->{HeaderRemark} = $SetHeaderRemark;
	}
	else {
		$self->{HeaderRemark} = $HeaderRemarkDefault;
	}
} # of sub SetHeaderRemark


sub SetAtomLocations { # to set the global variable AtomLocations
	my $self = shift;
	my $SetAtomLocations = shift;
	
	if (defined $SetAtomLocations) {
		$self->{AtomLocations} = uc ($SetAtomLocations);
	}
	else {
		$self->{AtomLocations} = uc ($AtomLocationsDefault);
	}
} # of sub SetAtomLocations


sub SetVerbose { # to set the global variable Verbose
	my $self = shift;
	my $Verbose = shift;
	
	if (defined $Verbose) {
		$self->{Verbose} = $Verbose;
	}
	else {
		$self->{Verbose} = $VerboseDefault;
	}
} # of sub SetVerbose


#######################################################################################################################
## Error Handling
#######################################################################################################################

sub Warning { # returns true if warnings were issued and false if not
	my $self = shift;
	return $self->{Warning}{Warning}
} # of sub Warning


sub Warning_NoENDMDL {
	my $self = shift;
	
	# if a MODEL without a corresponding ENDMDL has been found
	return $self->{Warning}{NoENDMDL};
} # of sub Warning_NoENDMDL


sub Warning_NoChainLabel {
	my $self = shift;
	
	# if a chain ID is given
	return $self->{Warning}{NoChainLabel};
} # of sub Warning_NoChainLabel


sub Warning_MultipleChainLabel {
	my $self = shift;
	
	# if a chain ID has been found more than once in a model
	return $self->{Warning}{MultipleChainLabel};
} # of sub Warning_MultipleChainLabel


sub Warning_ChainLabelsInvalid {
	my $self = shift;
	
	# if it was tried to access chains with external identifieres and invalid chain labels
	return $self->{Warning}{ChainLabelsInvalid};
} # of sub Warning_ChainLabelsInvalid


sub Warning_UnknownModel {
	my $self = shift;
	
	# if the requested model is not defined
	return $self->{Warning}{UnknownModel};
} # of sub Warning_UnknownModel


sub Warning_UnknownChain {
	my $self = shift;
	
	# if the requested chain is not defined
	return $self->{Warning}{UnknownChain};
} # of sub Warning_UnknownChain


sub Warning_UnknownResidue {
	my $self = shift;
	
	# if the requested residue is not defined
	return $self->{Warning}{UnknownResidue};
} # of sub Warning_UnknownResidue


sub Warning_UnknownAtom {
	my $self = shift;
	
	# if the requested atom is not defined
	return $self->{Warning}{UnknownAtom};
} # of sub Warning_UnknownAtom


sub Warning_UnknownLabel {
	my $self = shift;
	
	# if the requested label is not defined
	return $self->{Warning}{UnknownLabel};
} # of sub Warning_UnknownLabel


sub Warning_UnknownNumber {
	my $self = shift;
	
	# if the requested number is not defined
	return $self->{Warning}{UnknownNumber};
} # of sub Warning_UnknownNumber


sub Warning_UnknownAminoAcid {
	my $self = shift;
	
	# if a 1- or 3-letter-code given to AminoAcidConvert has not been recognized
	return $self->{Warning}{UnknownAminoAcid};
} # of sub Warning_UnknownAminoAcid


sub Warning_InvalidAminoAcid {
	my $self = shift;
	
	# if a code given to AminoAcidConvert has not 1 or 3 letters
	return $self->{Warning}{InvalidAminoAcid};
} # of sub Warning_InvalidAminoAcid


sub Warning_UnkownParameter {
	my $self = shift;
	
	# if a parameter given to one of the routines was not recognized 
	return $self->{Warning}{UnknownParameter};
} # of sub Warning_UnknownParameter


sub GetWarnings { # prints all warnings
	my $self = shift;
	
	if ($self->{Warning}{Warning}) {
		return @{$self->{WarningMsg}}
	}
} # of sub GetWarnings


sub PrintWarnings { # prints all warnings
	my $self = shift;
	
	if ($self->{Warning}{Warning}) {
		print "\nThe following warnings were reported:\n\n";
		
		foreach (@{$self->{WarningMsg}}) {
			print "$_";
		}
		
		return @{$self->{WarningMsg}}
	}
} # of sub PrintWarnings


#######################################################################################################################
#######################################################################################################################
## "Internal" subroutines
#######################################################################################################################
#######################################################################################################################


#######################################################################################################################
## BEGIN Check methods
#######################################################################################################################


sub _CheckChainLabel { # check for correct chain descriptor
	my $self = shift;
	my %args = @_;
	
	my ($Model, $Repair, $ChainIndex, $ChainLabel, %AllChainLabels);
	
	my $WarningNoChainLabel = 0;         # whether a chain ID error has already been reported
	my $WarningMultipleChainLabel = 0;   # whether a wrong sequence error has already been reported
	
	if (defined $args{Repair}) { $Repair = $args{Repair} }
	                      else { $Repair = 0 }
	
	$Model = $self->_GetModelOrDefault (\%args);
	
	$ChainIndex = $self->{Models}{$Model}{ChainIndex};
	
	$self->{Models}{$Model}{ChainLabelsValid} = 1;       # let's assume the best for a start...

	foreach $ChainLabel ( @{$ChainIndex} ) {
		# if even one single chain ID is not defined, the parser does not rely on them
        #if ($ChainLabel eq ' ' or not $ChainLabel) {
		if ($ChainLabel eq ' ' or !defined $ChainLabel) { # G.Postic (2014): "the chain label can be '0' (e.g. 2ZKR0)"
			# if a chain ID warning has not yet been reported
			if ( (not $WarningNoChainLabel) and (not $Repair) ) {
				$self->{Models}{$Model}{ChainLabelsValid} = 0;
				$WarningNoChainLabel = 1; # => warning is only reported once
				$self->_ReportWarning ("_CheckChainLabel", "NoChainLabel");
			}
		}
		# if this chain ID has already been found in another chain and is not a blank
		elsif ( $AllChainLabels{$ChainLabel} ) {
			if ( (not $WarningMultipleChainLabel) and (not $Repair) ) {
				$self->{Models}{$Model}{ChainLabelsValid} = 0;
				$WarningMultipleChainLabel = 1; # => warning is only reported once
				$self->_ReportWarning ("_CheckChainLabel", "MultipleChainLabel", "'$ChainLabel'");
			}
		}
		# if the chain label is defined and wasn't found before, save it as hash key to test for it later
		else {
			$AllChainLabels{$ChainLabel} = 1;
		}
	}

	if ( $Repair and (not $self->ChainLabelsValid) ) { $self->RenumberChains }
} # of sub _CheckChainLabel


#######################################################################################################################
## END Check Methods
#######################################################################################################################


#######################################################################################################################
## BEGIN Retrieve Methods for accessing the content
#######################################################################################################################

# The "Retrieve"-methods locate the requested content in the main PDB hash. Except the RetrieveModel methods they all
# all take two parameters, the $args reference to the arguments-hash of the calling method and a true or false value
# which tells the routine whether to filter the retrieved content or return it as found.

# If the retrieved content is filtered (depending on the second parameter), it is added to $args->{Content}
# and $args->{AtomIndex}. If it is not filtered, the "Retrieve"-method returns the reference to the original content
# and "translates" the identifiers, if necessary. That is, if 20 residues are in the first chain and "Residue" 23 is
# requested (without specifying a chain!), then chain 2 is returned and "Residue" changed to 3


sub _RetrieveAllModels { # combines models and all chains in one array
	my $self = shift;
	my $args = shift;
	my (@Models, $Model);
	
	@Models = $self->IdentifyModels;
	
	foreach $Model (@Models) {
		$args->{Model} = $Model;
		$self->_RetrieveModel ($args);
	}
	
	return;
} # of sub _RetrieveAllModels


sub _RetrieveModel { # adds a whole model to $args->{Content} and $args->{AtomIndex}
	my $self = shift;
	my $args = shift;
	my $DoNotFilter = shift;
	my $Model = $args->{Model};
	my (@Chains, $Chain);
	
	if (not defined $Model) { return undef }

	# if Atoms and Residues have not been indexed in this model yet, do it now
	if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
	
	if ( (not defined $Model) or (not defined $self->{Models}{$Model}) ) {
		$self->_ReportWarning ("_RetrieveModel", "UnknownModel", $Model);
		return undef;
	}
	
	@Chains = $self->IdentifyChains (Model => $Model);
	
	if (defined $self->{Models}{$Model}{MODEL}) {
		push @{$args->{Content}}, $self->{Models}{$Model}{MODEL};
		push @{$args->{AtomIndex}}, {
		                              Race => "MODEL",
		                              ModelNumber => $self->{Models}{ModelIndex}[$Model],
		                              Model => $Model,
		                            };
	}
	
	foreach $Chain (@Chains) {
		$args->{Chain} = $Chain;
		$self->_FilterAllLines ($args, $self->{Models}{$Model}{$Chain});
	}
	
	if (defined $self->{Models}{$Model}{ENDMDL}) {
		push @{$args->{Content}}, $self->{Models}{$Model}{ENDMDL};
		push @{$args->{AtomIndex}}, { Race => "ENDMDL" };
	}
	
	return;
} # of sub _RetrieveModel


sub _RetrieveChain { # adds a chain with an internal number to $args->{Content} and $args->{AtomIndex}
	my $self = shift;
	my $args = shift;
	my $DoNotFilter = shift;
	my ($Model, $Chain);
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	# This should already be caught by the calling method, but just in case...
	if (defined $args->{ChainLabel} and not $self->ChainLabelsValid) {
		$self->_ReportWarning ("_RetrieveChain", "ChainLabelsInvalid");
		return undef;
	}
	
	$Model = $self->_GetModelOrDefault ($args);
	$Chain = $self->_GetArgsInternal ($args, "Chain");
	
	# if Atoms and Residues have not been indexed in this model yet, do it now
	if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
	
	if (not defined $Chain) {
		return;
	}
	
	if ($args->{Chain} !~ m/\d/g) {
		$self->_ReportWarning ("_RetrieveChain", "Misc", "Keyword \"Chain\" is only for numeric identifiers. Use \"ChainLabel\" for the chain letters.");
		return;
	}
	
	if ($DoNotFilter) {
		return $self->{Models}{$Model}{$Chain};
	}
	else {
		$self->_FilterAllLines ($args, $self->{Models}{$Model}{$Chain});
	}
} # of sub _RetrieveChain


sub _RetrieveResidue { # adds a residue to $args->{Content} and $args->{AtomIndex}
	my $self = shift;
	my $args = shift;
	my $DoNotFilter = shift;
	
	my ($Model, @Chains, $Chain, $Residue, $Atom);
	my ($ResidueIndex, $CurResidues, $PrevResidues);
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	# This should already be caught by the calling method, but just in case...
	if (defined $args->{ChainLabel} and not $self->ChainLabelsValid) {
		$self->_ReportWarning ("_RetrieveResidue", "ChainLabelsInvalid");
		return undef;
	}
	
	$Model   = $self->_GetModelOrDefault ($args);
	$Chain   = $self->_GetArgsInternal ($args, "Chain");
	$Residue = $self->_GetArgsInternal ($args, "Residue");
	
	# if Atoms and Residues have not been indexed in this model yet, do it now
	if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
	
	if (not defined $Residue) { # if the residue has not been found
		return undef;
	}
	
	if (defined $Residue and $Residue !~ m/\d/g) {
		$self->_ReportWarning ("_RetrieveResidue", "Misc", "Keyword \"Residue\" is only for numeric identifiers. Use \"ResidueLabel\" for the chain letters.");
		return undef;
	}
	
	# "Residue" is relative to "Model" and "Chain"
	if (defined $Chain) { # the easy case, a chain is defined
		# check whether the given number is larger than the present residue number
		
		if ($Residue > $#{$self->{Models}{$Model}{$Chain}{ResidueIndex}}) {
			$self->_ReportWarning ("_RetrieveResidue", "UnknownResidue", $Residue);
			return undef;
		}
		
		if ($DoNotFilter) {
			return $self->{Models}{$Model}{$Chain};
		}
		else {
			$self->_FilterAllLines ($args, $self->{Models}{$Model}{$Chain});
			if ( (not defined $args->{AtomIndex}) or ($#{$args->{AtomIndex}} < 0) ) { return undef }
			
			return undef;
		}
	}
	else { # if no chain is defined
		@Chains = $self->IdentifyChains (Model => $Model);
		$PrevResidues = 0;
		
		foreach $Chain (@Chains) {
			$CurResidues = scalar @{$self->{Models}{$Model}{$Chain}{ResidueIndex}};
			
			if ($Residue > $PrevResidues + $CurResidues - 1) {
				$PrevResidues = $PrevResidues + $CurResidues;
				next;
			}
			
			$ResidueIndex = $self->{Models}{$Model}{$Chain}{ResidueIndex};
			
			# Convert $Residue, the wanted residue in this chain is $Residue minus
			# the number of residues in the previous chains
			$Residue = $Residue - $PrevResidues;
			
			# save the "converted" residue number for the filter routine
			$args->{Residue} = $Residue;
			
			if ($DoNotFilter) {
				return $self->{Models}{$Model}{$Chain};
			}
			else {
				$self->_FilterAllLines ($args, $self->{Models}{$Model}{$Chain});
				return;
			}
		}
	}
	
	$self->_ReportWarning ("_RetrieveResidue", "UnknownResidue", $Residue);
} # of sub _RetrieveResidue


sub _RetrieveAtom { # adds a requested atom to $args->{Content} and $args->{AtomIndex}
	my $self = shift;
	my $args = shift;
	my $DoNotFilter = shift;
	
	my ($Model, $Chain, @Chains, $Residue, $CurResidue, $Atom);
	my ($ChainHash, $ResidueHash, $CurAtoms, $PrevAtoms);
	
	# Check whether an external chain label was specified and block this at once if the chain labels are invalid.
	# This should already be caught by the calling method, but just in case...
	if (defined $args->{ChainLabel} and not $self->ChainLabelsValid) {
		$self->_ReportWarning ("_RetrieveAtom", "ChainLabelsInvalid");
		return undef;
	}
	
	$Model   = $self->_GetModelOrDefault ($args);
	$Chain   = $self->_GetArgsInternal ($args, "Chain");
	$Residue = $self->_GetArgsInternal ($args, "Residue");
	$Atom    = $self->_GetArgsInternal ($args, "Atom");
	
	# if Atoms and Residues have not been indexed in this model yet, do it now
	if (not $self->{Models}{$Model}{AtomsIndexed}) { $self->_IndexAtoms ($Model) }
	
	if (defined $Chain and $Chain !~ m/\d/g) {
		$self->_ReportWarning ("_RetrieveAtom", "Misc",
		                       "Keyword \"Chain\" is only for numeric identifiers. Use \"ChainLabel\" for the chain letters.");
		return undef;
	}
	
	if (defined $Residue and $Residue !~ m/\d/g) {
		$self->_ReportWarning ("_RetrieveAtom", "Misc",
		                       "Keyword \"Residue\" is only for numeric identifiers. Use \"ResidueLabel\" for the chain letters.");
		return undef;
	}
	
	if (defined $Chain) {
		if (defined $Residue) { # if chain and residue are given
			$ChainHash = $self->_RetrieveResidue ($args, 1);
			
			# check if enough atoms are present in the residue
			if ($Atom > $#{$ChainHash->{ResidueIndex}[$Residue]{Atoms}}) {
				$self->_ReportWarning ("_RetrieveAtom", "UnknownAtom", $Atom);
				return undef;
			}
			
			# now determine the absolute atom number in that chain
			$PrevAtoms = 0;
			if ($Residue > 0) {
				for $CurResidue (0 .. $Residue - 1) {
					$CurAtoms = scalar @{$ChainHash->{ResidueIndex}[$CurResidue]{Atoms}};
					$PrevAtoms = $PrevAtoms + $CurAtoms;
				}
			}
			
			$Atom = $Atom + $PrevAtoms;
			$args->{Atom} = $Atom;
			delete $args->{Residue};
			
			$self->_FilterAllLines ($args, $self->{Models}{$Model}{$Chain});
			return;
		}
		else { # if a chain but no residue is given
			if ($Atom > $#{$self->{Models}{$Model}{$Chain}{AtomIndex}}) {
				$self->_ReportWarning ("_RetrieveAtom", "UnknownAtom", $Atom);
				return undef;
			}
			
			$self->_FilterAllLines ($args, $self->{Models}{$Model}{$Chain});
			return;
		}
	}
	else { # if no chain is defined
		if (defined $Residue) { # no chain but residue given (there cannot be a TER)
			$ChainHash = $self->_RetrieveResidue ($args, 1);
			$Residue = $args->{Residue}; # update $Residue, might be changed by _RetrieveResidue
			
			if (not defined $Residue or not defined $Atom) { return undef }
			
			# check if enough atoms are present in the residue
			if ($Atom > $#{$ChainHash->{ResidueIndex}[$Residue]{Atoms}}) {
				$self->_ReportWarning ("_RetrieveAtom", "UnknownAtom", $Atom);
				return undef;
			}
			
			# now determine the absolute atom number in that chain
			$PrevAtoms = 0;
			
			for $CurResidue (0 .. $Residue-1) {
				$CurAtoms = scalar @{$ChainHash->{ResidueIndex}[$CurResidue]{Atoms}};
				$PrevAtoms = $PrevAtoms + $CurAtoms;
			}
			
			$Atom = $Atom + $PrevAtoms;
			$args->{Atom} = $Atom;
			$args->{Residue} = undef;
			
			$self->_FilterAllLines ($args, $ChainHash);
			return;
		}
		else { # no chain and no residue given
			@Chains = $self->IdentifyChains (Model => $Model);
			
			$PrevAtoms = 0;
			foreach $Chain (@Chains) {
				# "Atom" starts counting at 0, but AtomCount states the absolute number of atoms
				# => for 3 atoms in a chain, Atom can be 0..2, i.e. one less than AtomCount.
				
				# AtomCount also states the absolute number of atoms with respect to whether a TER
				# is included in the chain or not
				
				$CurAtoms = $self->{Models}{$Model}{$Chain}{AtomCount} - 1;
				
				if ($Atom > $CurAtoms + $PrevAtoms) {
					# now the absolute number of atoms is needed, i.e. +1
					$PrevAtoms = $PrevAtoms + $CurAtoms + 1;
					next;
				}
				
				$Atom = $Atom - $PrevAtoms;
				
				$args->{Atom} = $Atom;
				$self->_FilterAllLines ($args, $self->{Models}{$Model}{$Chain});
				return;
			}
		}
	}
	
	$self->_ReportWarning ("_RetrieveAtom", "UnknownResidue", $Residue);
} # of sub _RetrieveAtom


sub _GetModelOrDefault { # returns either the defined model or a default value
	my $self = shift;
	my $args = shift;
	my ($Model, @ModelKeys);
	
	if ( (defined $args->{Model}) and (defined $args->{ModelNumber}) ) {
		$self->_ReportError ("BadParameter", "Model and ModelNumber");
	}
	
	if (defined $args->{Model}) {
		$Model = $args->{Model};
		
		if ($self->_ExistsModel (Model => $Model)) {
			return $Model;
		}
		else {
			$self->_ReportWarning ("_GetModelOrDefault", "UnknownModel", $Model);
			return undef;
		}
	}
	elsif (defined $args->{ModelNumber}) {
		$Model = $self->GetModel (%{$args});
		
		if (not defined $Model) { return undef }
		
		delete $args->{ModelNumber};     # remove the external key
		$args->{Model} = $Model;        # update the args-hash for the following routines
	}
	else {
		# if no model is specified, model 1 is taken by default
		# if verbose mode is enabled and more than one model is present,
		# a warning is issued
		@ModelKeys = $self->_GetNumericKeys ($self->{Models});
		
		if ( $self->{Verbose} and (scalar @ModelKeys > 1) ) {
			$self->_ReportWarning ("_GetModelOrDefault", "ModelDefault");
		}
		
		$args->{Model} = 0;  # add the default value to the argument hash
		return 0;            # and also return it
	}
} # of sub _GetModelOrDefault


sub _GetArgsInternal { # returns the internal identifier number from a given hash
	my $self = shift;
	my $args = shift;
	my $Internal = shift;
	my ($External, $Value);
	
	if    ($Internal eq "Model")   { $External = "ModelNumber"   }
	elsif ($Internal eq "Chain")   { $External = "ChainLabel"    }
	elsif ($Internal eq "Residue") { $External = "ResidueNumber" }
	elsif ($Internal eq "Atom")    { $External = "AtomNumber"    }
	
	if ( (defined $args->{$Internal}) and (defined $args->{$External}) ) {
		$self->_ReportError ("BadParameter", "$Internal and $External");
	}
	
	if (defined $args->{$Internal}) {
		return $args->{$Internal};
	}
	elsif (defined $args->{$External}) {
		if    ($Internal eq "Model")   { $Value = $self->GetModel   (%{$args}) }
		elsif ($Internal eq "Chain")   { $Value = $self->GetChain   (%{$args}) }
		elsif ($Internal eq "Residue") { $Value = $self->GetResidue (%{$args}) }
		elsif ($Internal eq "Atom")    { $Value = $self->GetAtom    (%{$args}) }
		
		if (not defined $Value) { return undef }
		
		delete $args->{$External};     # remove the external key
		$args->{$Internal} = $Value;   # update the args-hash for the following routines
		return $args->{$Internal};
	}
	else {
		return undef;
	}
} # of sub _GetArgsInternal


sub _ExistsModel {
	my $self = shift;
	my %args = @_;
	
	my (@Models, $Model);
	
	@Models = $self->IdentifyModels;
	$Model = $args{Model};
	
	if ($self->_InArray (\@Models, $Model) ) { return 1 }
	                                    else { return 0 }
} # of sub _ExistsModel


sub _ExistsChain {
	my $self = shift;
	my %args = @_;
	my (@Chains, $Chain);
	
	@Chains = $self->IdentifyChains (%args);
	$Chain = $args{Chain};
	
	if ($self->_InArray (\@Chains, $Chain) ) { return 1 }
	                                    else { return 0 }
} # of sub _ExistsChain


#######################################################################################################################
## END Methods for accessing the content
#######################################################################################################################


#######################################################################################################################
## BEGIN Filter Methods
#######################################################################################################################

sub _FilterAllLines { # filters the ATOM lines in the a given array (lines not starting with ATOM et.al. are also returned)
	my $self = shift;
	my $args = shift;
	my $Hash = shift;
	
	my ($Element, $FilterOutElement, $ResidueLabel, $FilterOutResidueLabel);
	my ($AtomType, $FilterOutAtomType, $Race, $FilterOutRace, $Calcium);
	my ($CurLine, $AtomPattern, @SearchArray, $FilterAtomLocations, $Code);
	my ($AtomLocations, $AltLocIndicator, $AltLocField, $AltLocPattern);
	my ($OldResidueNumber, $CurResidueNumber, $Residue, $NewResidue);
	
	# _FilterAllLines is given an array with the respective chain. However, it does not know
	# *which* chain this is (1st, 2nd, ...). Thus relative identifiers for a residue or an atom
	# have already been "converted" by _RetrieveResidue or _RetrieveAtom.
	# This means, if the 1st chain has 10 residues, and residue 12 is requested (without specifying
	# a chain), chain 2 will be forwarded to _FilterAllLines with Residue => 2.
	
	# Atom indeces are also converted to absolute indeces, if they had been relative to a residue before.
	# That is, if residue 0 has 5 atoms and "Atom 0" in "Residue 1" is requested, "Residue 1" has been
	# deleted and "Atom" was converted 5 (since 0..4 are in the first residue).
	
	if (defined $args->{Residue}) {
		if (not defined $Hash->{ResidueIndex}[ $args->{Residue} ]) { return }
		
		@SearchArray = ( @{$Hash->{ResidueIndex}[$args->{Residue}]{Atoms}} );
	}
	elsif (defined $args->{Atom}) {
		push @SearchArray, $args->{Atom};
	}
	else {
		@SearchArray = (0 .. $#{$Hash->{Content}});
	}
	
	####################################################################################################
	# Read the given parameters and adjust some things for the filter process
	####################################################################################################
	
	if (defined $args->{Element}) {
		$Element = $args->{Element};
		$Element = uc $Element; # upper case
		
		if ($Element =~ m/^-/) {
			$FilterOutElement = 1;
			$Element = substr ($Element, 1);
		}
		
		if (defined $ELEMENTFILTERS{$Element}) {
			$Element = $ELEMENTFILTERS{$Element};
		}
		else {
			$self->_ReportError ("UnknownElement", $Element);
		}
	}
		
	if (defined $args->{ResidueLabel}) {
		$ResidueLabel = $args->{ResidueLabel};
		
		if ($ResidueLabel =~ m/^-/) {
			$FilterOutResidueLabel = 1;
			# take everything from position 1 till the end (cut away the "-")
			$ResidueLabel = substr ($ResidueLabel, 1);
		}
		
		if (length ($ResidueLabel) == 1) {
			$ResidueLabel = $self->ResidueLabelConvert ($ResidueLabel)
		}
	}
		
	if (defined $args->{AtomType}) {
		$AtomType = $args->{AtomType};
		
		if ($AtomType =~ m/^-/) {
			$FilterOutAtomType = 1;
			$AtomType = substr ($AtomType, 1);
		}
		
		if ($AtomType eq "Ca") {
			$AtomType = "CA";
			$Calcium = 1;
		}
	}
		
	if (defined $args->{Race}) {
		$Race = $args->{Race};
		
		if ($Race =~ m/^-/) {
			$FilterOutRace = 1;
			$Race = substr ($Race, 1);
		}
	}
	
	
	####################################################################################################
	# Prepare the filtering of alternative atom locations
	####################################################################################################
	
	$AtomLocations   = $args->{AtomLocations};
	$AtomLocations   = uc $AtomLocations;
	$AltLocIndicator = $self->{AltLocIndicator};
	
	# 0 = no indicator found
	# 1 = col 16 with letters
	# 2 = col 16 with digits
	# 3 = col 13 with digits
	# 4 = col 15 with digits
	
	if ( (defined $AtomLocations) and ( ($AtomLocations ne "ALL") and ($AltLocIndicator > 0) ) ) {
		$FilterAtomLocations = 1;
		
		if ( ($AltLocIndicator == 1) or ($AltLocIndicator == 2) ) {
			$AltLocField = "AltLoc";
		}
		
		if ( ($AltLocIndicator == 3) or ($AltLocIndicator == 4) ) {
			$AltLocField = "AtomType";
		}
		
		if ($AtomLocations eq "FIRST") {
			# get only atoms with 'A' or no alternate location
			if ($AltLocIndicator == 1) { $AltLocPattern = "(A| )" }
			
			# get only atoms with '1' or no alternate location
			if ($AltLocIndicator == 2) { $AltLocPattern = "(1| )" }
			
			# get only atoms with '1' or no alternate location
			if ($AltLocIndicator == 3) { $AltLocPattern = "^(1| )[A-Za-z]" }
			
			# get only atoms with '1' or no alternate location
			if ($AltLocIndicator == 4) { $AltLocPattern = "[1-9](1| )" }
		}
		elsif ($AtomLocations eq "NONE")  {
			if ($AltLocIndicator == 1) { $AltLocPattern = "" }
			
			if ($AltLocIndicator == 2) { $AltLocPattern = "" }
			
			# get only atoms with '1' or no alternate location
			if ($AltLocIndicator == 3) { $AltLocPattern = "^[A-Za-z]" }
			
			# get only atoms with '1' or no alternate location
			if ($AltLocIndicator == 4) { $AltLocPattern = "[A-Za-z][1-9]" }
		}
		else  {
			$AltLocPattern = $AtomLocations; # take $AtomLocations as it was entered
		}
	}
	else {
		$FilterAtomLocations = 0;
	}
	
	
	####################################################################################################
	# start the filtering
	####################################################################################################
	
	$AtomPattern = "$self->{AllAtoms}|TER";
	
	# this is compared to the residue number of the current atom, to check whether a new residue began
	$OldResidueNumber = "-1000";
	
	# this denotes the current residue in $Hash->{ResidueIndex} to add it to $args->{ResidueIndex}
	$Residue = -1;
	
	foreach $CurLine ( @SearchArray ) {
		if ($Hash->{AtomIndex}[$CurLine]{Race} =~ m/^($AtomPattern)/) {
			
			if (defined $args->{ResidueIndex}) {
				# if the ResidueIndex has to be returned later on, it needs to be kept track of the index of the current
				# residue in the ResdiueIndex of the main hash
				
				# combine ResidueNumber and InsResidue to catch cases like residue 2 and 2A
				$CurResidueNumber = $Hash->{AtomIndex}[$CurLine]{ResidueNumber} . $Hash->{AtomIndex}[$CurLine]{InsResidue};
				$NewResidue = 0;
				
				if ($Hash->{AtomIndex}[$CurLine]{Race} ne "TER" and $OldResidueNumber ne $CurResidueNumber) {
					$OldResidueNumber = $CurResidueNumber;
					$NewResidue = 1;
					++$Residue;
				}
			}
			
			if ($Element) {
				if ($FilterOutElement) {
					if ($Hash->{AtomIndex}[$CurLine]{AtomType} =~ m/$Element/) { next }
				}
				else {
					if (not $Hash->{AtomIndex}[$CurLine]{AtomType} =~ m/$Element/) { next }
				}
			}
			
			if ($ResidueLabel) {
				if ($FilterOutResidueLabel) {
					if ($Hash->{AtomIndex}[$CurLine]{ResidueLabel} =~ m/$ResidueLabel/) { next }
				}
				else {
					if (not $Hash->{AtomIndex}[$CurLine]{ResidueLabel} =~ m/$ResidueLabel/) { next }
				}
			}
			
			if ($AtomType) {
				if ($FilterOutAtomType) {
					if ($Hash->{AtomIndex}[$CurLine]{AtomType} =~ m/^$AtomType$/) {
						if ($Calcium) { # if only calciums are affected check whether it is a HETATM
							if ($Hash->{AtomIndex}[$CurLine]{Race} eq "HETATM") { next }
						}
						else { # if all atoms are affected and this AtomType matches
							next;
						}
					}
				}
				else {
					if ($Hash->{AtomIndex}[$CurLine]{AtomType} !~ m/^$AtomType$/) { next }
					
					if ($Calcium and $Hash->{AtomIndex}[$CurLine]{Race} eq "ATOM") { next }
				}
			}
			
			if ($Race) {
				if ($FilterOutRace) {
					if ($Hash->{AtomIndex}[$CurLine]{Race} eq $Race) { next }
				}
				else {
					if ($Hash->{AtomIndex}[$CurLine]{Race} ne $Race) { next }
				}
			}
			
			if ($FilterAtomLocations) {
				if ($Hash->{AtomIndex}[$CurLine]{$AltLocField} !~ m/$AltLocPattern/) { next }
			}
			
			# if the loop got until here (without a "next"), the line is OK
			push @{$args->{Content}}, $Hash->{Content}[$CurLine];
			push @{$args->{AtomIndex}}, $Hash->{AtomIndex}[$CurLine];
			
			if (defined $args->{ResidueIndex}) {
				if ($NewResidue) { # if a new residue has begun 
					# Create a copy of the ResidueIndex (pushing the reference would return the main hash's
					# variable to the user and result in corruption of the data in the memory)
					my %ResidueIndexCopy = %{$Hash->{ResidueIndex}[$Residue]};
					
					push @{$args->{ResidueIndex}}, \%ResidueIndexCopy;
					
					if ($args->{OneLetterCode}) {
						$Code = $self->AminoAcidConvert ($Hash->{ResidueIndex}[$Residue]{ResidueLabel});
						$args->{ResidueIndex}[ $#{$args->{ResidueIndex}} ]{ResidueLetter} = $Code;
					}
				}
			}
			
		} # of if ($Hash->{AtomIndex}[$CurLine]{Race} =~ m/^($AtomPattern)/)
	} # of foreach $CurLine ( @SearchArray )
} # of sub _FilterAllLines


sub _FilterLines { # returns an array containing only the lines beginning with a given pattern
	my $self = shift;
	my $args = shift;
	my $Pattern = $args->{Pattern};
	my (@NewContent, $Line);
	
	if (not defined $Pattern) { return } # if no Pattern is given, all lines are returned again
	
	$Line = 0;
	
	while ($Line <= $#{$args->{Content}}) {
		if ($args->{Content}[$Line] =~ m/^$Pattern/) { # if the line does not start with the pattern
			push @NewContent, $args->{Content}[$Line];
		}
		
		++$Line;
	} # of while ($Line <= $#{$args->{Content}}
	
	$args->{Content} = \@NewContent;
} # of sub _FilterLines


sub _DetermineAltLocIndicator { # check which alternate location indicator is used
	my $self = shift;
	my $Content = shift;
	
	my ($CurType, $Line);
	
	foreach $Line (@{$Content}) {
		# if it is no atom line, go to next line, some HETATMs have two-digit numbers,
		# so at least the third check could return wrong results
		if ($Line !~ m/^ATOM/) { next }
		
		# if column 16 contains something else than a blank
		if (substr ($Line, 16, 1) ne " ") {
			if (substr ($Line, 16, 1) =~ m/[A-Za-z]/) {
				$self->{AltLocIndicator} = 1;
				last;
			}
			
			if (substr ($Line, 16, 1) =~ m/[0-9]/) {
				$self->{AltLocIndicator} = 2;
				last;
			}
		}
		
		#     The following "possibilities" of assigning the alternative atom locations
		#     were found in some really crappy files. Determination can lead to caveats
		#     and conflicts with "correct" hydrogen names and was therefore commented out.
		#
		#		$CurType = substr ($Line, 12, 4);
		#		
		#		# if the type field starts with a digit between 1 and 9, it is
		#		# assumed that the first is the alternate location indicator
		#		if ($CurType =~ m/^[1-9]/) {
		#			$self->{AltLocIndicator} = 3;
		#			last;
		#		}
		#		
		#		# if the type field ends with two digits, it is assumed that
		#		# the last is the alternate location indicator
		#		if ($CurType =~ m/[1-9][1-9]$/) {
		#			$self->{AltLocIndicator} = 4;
		#			last;
		#		}
	} # of foreach $Line (@{$Content})
	
	# if no alternate location indicator has been found
	if (not defined $self->{AltLocIndicator}) {
		$self->{AltLocIndicator} = 0;
	}
} # of sub _DetermineAltLocIndicator


#######################################################################################################################
## END Filter Methods
#######################################################################################################################


#######################################################################################################################
## BEGIN Calculation of the angles
#######################################################################################################################


sub _DihedralAngle { # calculates the dihedral angle between two residues (i.e. 4 atoms)
	my $self = shift;
	my $RefAtom1 = shift;
	my $RefAtom2 = shift;
	my $RefAtom3 = shift;
	my $RefAtom4 = shift;
	
	if (not defined $RefAtom1 or not defined $RefAtom2 or
		 not defined $RefAtom3 or not defined $RefAtom4) {
		return undef;
	}
	
	my %Atom1 = %{$RefAtom1};
	my %Atom2 = %{$RefAtom2};
	my %Atom3 = %{$RefAtom3};
	my %Atom4 = %{$RefAtom4};
	
	my $Result = 360.000;
	
	my (%Vector12, %Vector23, %Vector43);
	my (%CrossP, %CrossX, %CrossY, $DotX, $DotY, $DotPX, $DotPY);
	
	my $Pi = 3.141593;
	my $Radian = 57.29578;  # = 180 / Pi
	
	%Vector12 = $self->_Diff (\%Atom1, \%Atom2);
	%Vector23 = $self->_Diff (\%Atom2, \%Atom3);
	%Vector43 = $self->_Diff (\%Atom4, \%Atom3);
	
	%CrossP = $self->_CrossProduct (\%Vector23, \%Vector12);
	%CrossX = $self->_CrossProduct (\%Vector23, \%Vector43);
	%CrossY = $self->_CrossProduct (\%Vector23, \%CrossX);
	
	$DotX = $self->_DotProduct (\%CrossX, \%CrossX);
	$DotY = $self->_DotProduct (\%CrossY, \%CrossY);
	
	if ( ($DotX <= 0.0) or ($DotY <= 0.0) ) {
		return sprintf ("%.3f", $Result);
	}
	
	$DotPX = $self->_DotProduct (\%CrossP, \%CrossX) / sqrt ($DotX);
	$DotPY = $self->_DotProduct (\%CrossP, \%CrossY) / sqrt ($DotY);
	
	if ( ($DotPX == 0.0) or ($DotPY == 0.0) ) {
		return sprintf ("%.3f", $Result);
	}
	
	$Result = atan2 ($DotPY, $DotPX) * $Radian;
	
	return sprintf ("%.3f", $Result);
} # of sub _DihedralAngle


sub _CrossProduct { # calculates the cross product of two vectors
	my $self = shift;
	my $Vector1 = shift;
	my $Vector2 = shift;
	my %CrossProduct;
	
	$CrossProduct{x} =   $Vector1->{y} * $Vector2->{z}
	                   - $Vector2->{y} * $Vector1->{z};
	
	$CrossProduct{y} =   $Vector1->{z} * $Vector2->{x}
	                   - $Vector2->{z} * $Vector1->{x};
	
	$CrossProduct{z} =   $Vector1->{x} * $Vector2->{y}
	                   - $Vector2->{x} * $Vector1->{y};
	
	return %CrossProduct;
} # of sub _CrossProduct


sub _DotProduct { # calculates the dot product of two vectors
	my $self = shift;
	my $Vector1 = shift;
	my $Vector2 = shift;
	my $DotProduct;
	
	$DotProduct = $Vector1->{x} * $Vector2->{x} +
	              $Vector1->{y} * $Vector2->{y} +
	              $Vector1->{z} * $Vector2->{z};
	
	return $DotProduct
} # of sub _DotProduct


sub _Diff {
	my $self = shift;
	my $Vector1 = shift;
	my $Vector2 = shift;
	my %ResultVector;
	
	if (not defined $Vector1 or not defined $Vector2) { return undef }
	
	$ResultVector{x} = $Vector1->{x} - $Vector2->{x};
	$ResultVector{y} = $Vector1->{y} - $Vector2->{y};
	$ResultVector{z} = $Vector1->{z} - $Vector2->{z};
	
	return %ResultVector;
} # of sub _Diff


sub _Distance {
	my $self = shift;
	my $Vector1 = shift;
	my $Vector2 = shift;
	my $Distance;
	
	if (not defined $Vector1 or not defined $Vector2) { return undef }
	
	$Distance = sqrt ( ($Vector1->{x} - $Vector2->{x})**2
	                 + ($Vector1->{y} - $Vector2->{y})**2
	                 + ($Vector1->{z} - $Vector2->{z})**2 );
	
	return $Distance;
} # of sub _Distance

#######################################################################################################################
## END Calculation of the dihedral angles
#######################################################################################################################


sub _IndexChains { # indexes the ModelNumbers and ChainLabels
	my $self = shift;
	my (@ModelKeys, @ModelIndex, $ModelNumber, $Model);
	my (@ChainKeys, @ChainIndex, $ChainLabel, $Chain);
	
	# returns e.g. 0, 1, 2 when 3 Models are present
	@ModelKeys = $self->_GetNumericKeys ($self->{Models});
	
	foreach $Model (@ModelKeys) {
		if (defined $self->{Models}{$Model}{MODEL}) {
			$ModelNumber = substr ($self->{Models}{$Model}{MODEL}, 10, 4);
			$ModelNumber =~ s/\s//g;
		}
		else {
			$ModelNumber = '';
		}
		
		my @ChainIndex;
		
		$ModelIndex[$Model] = $ModelNumber;
		
		@ChainKeys = $self->_GetNumericKeys ($self->{Models}{$Model});
		
		foreach $Chain (@ChainKeys) {
			push @ChainIndex, $self->{Models}{$Model}{$Chain}{ChainLabel};
		}
		
		$self->{Models}{$Model}{ChainIndex} = \@ChainIndex;
	
		# determine the "status" of the chain IDs in the file
		$self->_CheckChainLabel (Model => $Model, Repair => 0);
	}
	
	$self->{Models}{ModelIndex} = \@ModelIndex;
} # of sub _IndexChains


sub _IndexResidues { # indexes only the residues after atoms have been deleted from the hash
	my $self  = shift;
	my $Model = shift;
	
	my (@Models, @Chains, $Chain);
	my ($Race, $OldInsResidue, $ResidueLabel, $ChainLabel, $AtomIndex);
	my ($CurResidueNumber, $OldResidueNumber, $ResidueIndex, $InsResidue);
	my ($Content, $LineCount, $AtomCount);
	
	if (not defined $Model or not defined $self->{Models}{$Model} ) {
		return undef;
	}
	
	@Chains = $self->IdentifyChains (Model => $Model);
	
	foreach $Chain (@Chains) {
		# reinitialize for each chain
		my @ResidueIndex;
		
		$Content = $self->{Models}{$Model}{$Chain}{Content};
		$AtomIndex = $self->{Models}{$Model}{$Chain}{AtomIndex};
		
		$ResidueIndex = -1;
		$OldResidueNumber = -1;
		$OldInsResidue = "";
		
		$AtomCount = 0;
		
		for $LineCount (0 .. $#{$AtomIndex}) {
			$Race = $AtomIndex->[$LineCount]{Race};
			
			if ($Race ne "TER") { ++$AtomCount }
			
			$CurResidueNumber = $AtomIndex->[$LineCount]{ResidueNumber};
			$ResidueLabel = $AtomIndex->[$LineCount]{ResidueLabel};
			$InsResidue = $AtomIndex->[$LineCount]{InsResidue};
			
			# if the line belongs to a new residue
			# string comparisons used, since in some crappy files the numbers are missing
			# (especially in TER lines)
			if ( ($CurResidueNumber ne $OldResidueNumber)
					or
					( ($CurResidueNumber eq $OldResidueNumber) and ($InsResidue ne $OldInsResidue) ) )
			{
				++$ResidueIndex;
				
				# add the atom to the residue array, a TER does NOT belong to the residue, thus the
				# number of atoms in a residue can be determined by $#{$Residue{Atoms}}+1.
				if ($Race ne "TER") {
					push @{$ResidueIndex[$ResidueIndex]->{Atoms}}, $LineCount;
				
					# add its number and label
					$ResidueIndex[$ResidueIndex]->{ResidueNumber} = $CurResidueNumber;
					$ResidueIndex[$ResidueIndex]->{ResidueLabel} = $ResidueLabel;
					$ResidueIndex[$ResidueIndex]->{InsResidue} = $InsResidue;
					$OldResidueNumber = $CurResidueNumber;
					$OldInsResidue = $InsResidue;
				}
			}
			else { # if the line belongs to the same residue
				if ($Race ne "TER") { push @{$ResidueIndex[$ResidueIndex]->{Atoms}}, $LineCount }
			}
		}
		
		# A TER does not count as atom, thus the absolute number of atoms in a chain cannot be simplily
		# determined bz $#{$Chain}, but depends on whether a TER is included or not. To avoid this check,
		# the absolute number of atoms is saved in {AtomCount}. A TER does not belong to a residues, thus
		# in the case of a residue, the $# can be used to determine the AtomCount.
		# ANISOU, SIGATM and SIGUIJ are regarded as "full" atoms, otherwise an array of some lines had to
		# be returned, if the atom is requested.
		$self->{Models}{$Model}{$Chain}{ResidueIndex} = \@ResidueIndex;
		$self->{Models}{$Model}{$Chain}{ResidueCount} = $ResidueIndex + 1;
		$self->{Models}{$Model}{$Chain}{AtomCount} = $AtomCount;
	}

	return 1;
} # of _IndexResidues


sub _IndexAtoms { # indexes the ResidueNumbers and AtomLabels
	my $self = shift;
	my $Model = shift;
	
	my (@Models, @Chains, $Chain);
	my ($Race, $AtomNumber, $AtomType, $AltLoc, $OldInsResidue, $ResidueLabel, $ChainLabel);
	my ($CurResidueNumber, $OldResidueNumber, $ResidueIndex, $InsResidue, $DisableAtomTypes);
	my ($xCoord, $yCoord, $zCoord, $Occupancy, $Temp, $Rest);
	my ($Content, $LineCount, $AtomCount, $TER);
	
	if (not defined $Model or not defined $self->{Models}{$Model} ) {
		return undef;
	}
	
	@Chains = $self->IdentifyChains (Model => $Model);
	
	foreach $Chain (@Chains) {
		# reinitialize for each chain
		my (@AtomIndex, @ResidueIndex);
		
		$Content = $self->{Models}{$Model}{$Chain}{Content};
		
		$ResidueIndex     = -1;
		$OldResidueNumber = -1;
		$OldInsResidue    = "";
		$DisableAtomTypes =  0;  # true, if multiple atoms with the same atom type are found in a residue
		
		$AtomCount = 0;
		
		for $LineCount (0 .. $#{$Content}) {
			if ($Content->[$LineCount] !~ m/^($self->{AllAtoms}|TER)/) { next }
			
			my %CurIndex;
			
			# with to identify to which chain an atom belongs in case no
			# chainID is present (used e.g. by GetAngles)
			$CurIndex{Chain} = $Chain;
			
			# It has been tried whether the routine could be speeded up by replacing the
			# numerous substr commands with one unpack command (many google threads write this).
			# However, in each benchmark the unpack command was slightly slower and is of course
			# much more confusing.
			
#				 ($Race, $AtomNumber, $AtomType, $AltLoc, $ResidueLabel, $ChainLabel,
#				 $CurResidueNumber, $InsResidue, $xCoord, $yCoord, $zCoord, $Occupancy, $Temp, $Rest) =
#				 unpack ("A6 A5 x1 A4 A1 A3 x1 A1 A4 A1 x3 A8 A8 A8 A5 A5 A13", $Content->[$LineCount]);
#				 dp ($Race, $AtomNumber, $AtomType, $AltLoc, $ResidueLabel, $ChainLabel,
#				 $CurResidueNumber, $InsResidue, $xCoord, $yCoord, $zCoord, $Occupancy, $Temp, $Rest);
			
			$Race             = substr ($Content->[$LineCount],  0, 6);
			$AtomNumber       = substr ($Content->[$LineCount],  6, 5);
			$AtomType         = substr ($Content->[$LineCount], 12, 4);
			$AltLoc           = substr ($Content->[$LineCount], 16, 1);
			$ResidueLabel     = substr ($Content->[$LineCount], 17, 3);
			$ChainLabel       = substr ($Content->[$LineCount], 21, 1);
			$CurResidueNumber = substr ($Content->[$LineCount], 22, 4);
			$InsResidue       = substr ($Content->[$LineCount], 26, 1);
			$xCoord           = substr ($Content->[$LineCount], 30, 8);
			$yCoord           = substr ($Content->[$LineCount], 38, 8);
			$zCoord           = substr ($Content->[$LineCount], 46, 8);
			$Occupancy        = substr ($Content->[$LineCount], 54, 5);
			$Temp             = substr ($Content->[$LineCount], 61, 5);
			$Rest             = substr ($Content->[$LineCount], 67,13);
			
			$Race             =~ s/\s//g;
			$AtomNumber       =~ s/\s//g;
			$AtomType         =~ s/\s//g;
			# $AltLoc           =~ s/\s//g;
			$ResidueLabel     =~ s/\s//g;
			# $ChainLabel       =~ s/\s//g;
			$CurResidueNumber =~ s/\s//g;
			# $InsResidue       =~ s/\s//g;
			$xCoord           =~ s/\s//g;
			$yCoord           =~ s/\s//g;
			$zCoord           =~ s/\s//g;
			$Occupancy        =~ s/\s//g;
			$Temp             =~ s/\s//g;
			# $Rest             =~ s/\s//g;
			
			%CurIndex = (
				Race => $Race,
				AtomNumber => $AtomNumber,
				AtomType => $AtomType,
				AltLoc => $AltLoc,
				ResidueLabel => $ResidueLabel,
				ChainLabel => $ChainLabel,
				ResidueNumber => $CurResidueNumber,
				InsResidue => $InsResidue,
				x => $xCoord,
				y => $yCoord,
				z => $zCoord,
				Occupancy => $Occupancy,
				Temp => $Temp,
				Rest => $Rest
			);
			
			if ($Race ne "TER") { ++$AtomCount }
			
			# if the line belongs to a new residue
			# string comparisons used, since in some crappy files the numbers are missing
			# (especially in TER lines)
			if ( ($CurResidueNumber ne $OldResidueNumber)
					or
					( ($CurResidueNumber eq $OldResidueNumber) and ($InsResidue ne $OldInsResidue) ) )
			{
				++$ResidueIndex;
				$DisableAtomTypes = 0; # reset for the new residue
				
				# add the atom to the residue array, a TER does NOT belong to the residue, thus the
				# number of atoms in a residue can be determined by $#{$Residue{Atoms}}+1.
				if ($Race ne "TER") {
					# create/extend the array of all atoms in that residue
					push @{$ResidueIndex[$ResidueIndex]->{Atoms}}, $LineCount;
					
					# create a hash using the atom types as keys
					if (not $ResidueIndex[$ResidueIndex]->{AtomTypes}{$AtomType}
					    and not $DisableAtomTypes) {
						$ResidueIndex[$ResidueIndex]->{AtomTypes}{$AtomType} = \%CurIndex;
					}
					else { # if the same atom type was added before, access by atom types is disabled for this residue
						delete $ResidueIndex[$ResidueIndex]->{AtomTypes};
						$DisableAtomTypes = 1;
					}
				
					# add its number and label
					$ResidueIndex[$ResidueIndex]->{ResidueNumber} = $CurResidueNumber;
					$ResidueIndex[$ResidueIndex]->{ResidueLabel} = $ResidueLabel;
					$ResidueIndex[$ResidueIndex]->{InsResidue} = $InsResidue;
					$OldResidueNumber = $CurResidueNumber;
					$OldInsResidue = $InsResidue;
				}
			}
			else { # if the line belongs to the same residue
				# the TER is not part of the residue
				if ($Race ne "TER") {
					# creat/extend  the array of all atoms in that residue
					push @{$ResidueIndex[$ResidueIndex]->{Atoms}}, $LineCount;
					
					# create a hash using the atom types as keys
					if (not $ResidueIndex[$ResidueIndex]->{AtomTypes}{$AtomType}
					    and not $DisableAtomTypes) {
						$ResidueIndex[$ResidueIndex]->{AtomTypes}{$AtomType} = \%CurIndex;
					}
					else { # if the same atom type was added before, access by atom types is disabled for this residue
						delete $ResidueIndex[$ResidueIndex]->{AtomTypes};
						$DisableAtomTypes = 1;
					}
				}
			}
			
			push @AtomIndex, \%CurIndex;
		}
		
		# A TER does not count as atom, thus the absolute number of atoms in a chain cannot be simply
		# determined by $#{$Chain}, but depends on whether a TER is included or not. To avoid this check,
		# the absolute number of atoms is saved in {AtomCount}. A TER does not belong to a residues, thus
		# in the case of a residue, the $# can be used to determine the AtomCount.
		# ANISOU, SIGATM and SIGUIJ are regarded as "full" atoms, otherwise an array of some lines had to
		# be returned, if the atom is requested.
		$self->{Models}{$Model}{$Chain}{ResidueIndex} = \@ResidueIndex;
		$self->{Models}{$Model}{$Chain}{AtomIndex} = \@AtomIndex;
		$self->{Models}{$Model}{$Chain}{ResidueCount} = $ResidueIndex + 1;
		$self->{Models}{$Model}{$Chain}{AtomCount} = $AtomCount;
	}
	
	$self->{Models}{$Model}{AtomsIndexed} = 1;
	
	# determine the "status" of the chain IDs in the file
	$self->_CheckChainLabel (Model => $Model, Repair => 0);
	
	return 1;
} # of sub _IndexAtoms


sub _GetSingleElement { # returns the element from an AtomType in a given AtomHash (a parsed atom line)
	my $self = shift;
	my $AtomHash = shift;
	my ($AtomType, $Key, $Element);
	
	if ($AtomHash->{Race} eq "TER") { return undef }
	
	$AtomType = $AtomHash->{AtomType};
	
	if (not defined $AtomType) {
		$self->_ReportWarning ("_GetSingleElement", "Misc", "AtomType is not defined in one atom and the mass could not be determined.");
		return undef;
	}
	
	if ($AtomType eq "CA") {
		if ($AtomHash->{Race} eq "ATOM") { $Element = "C"  }
	                               else { $Element = "CA" }
	}
	else {
		foreach $Key ( keys %ELEMENTFILTERS ) {
			if ($AtomType =~ m/^$ELEMENTFILTERS{$Key}$/) { return "$Key" }
		}
		
		# if no key matched the atom type
		$self->_ReportWarning ("_GetSingleElement", "Misc", "The Element $AtomType has not been recognized.");
		return undef;
	}
} # of sub _GetSingleElement


sub _CheckParameters { # checks for unknown keywords (to track typos)
	my $self = shift;
	my %args = @_;
	my (@Keywords, @Keys, $Key);
	
	@Keywords =
		qw/ Model ModelNumber Chain ChainLabel
		    Residue ResidueLabel ResidueNumber
		    Atom AtomNumber AtomType
		    ModelStart ChainStart ResidueStart AtomStart
		    CharmmFormat PDB2CHARMM CHARMM2PDB
		    SetChainLabel KeepInsertions Race IgnoreTER
		    Header MinHeader Footer MinFooter AtomLocation Element
		    FileName Content ChainSuffix ModelSuffix OneLetterCode
		    AtomLocations AtomIndex ResidueIndex GetReference Pattern
		/;
		
	@Keys = keys %args;
	
	foreach $Key (@Keys) {
		if (not $self->_InArray (\@Keywords, $Key) ) {
			$self->_ReportWarning ("_CheckParameters", "Misc", "$Key");
		}
	}
} # of sub _CheckParameters


sub _NumberConvert { # converts 1 into A, 2 into B, ...
	my $self = shift;
	my $Var = shift;
	
	if ($Var =~ m/[0-9]+/) {
		$Var = chr($Var+64); # A has the ASCII code 65
	}
	
	return $Var;
} # of sub _NumberConvert


sub _PDB2CHARMM { # converts PDB format to CHARMM
	my $self = shift;
	my %args = @_;
	
	if (not defined $self->{Models}) { $self->Parse }
	if (not defined $args{Content} or not defined $args{AtomIndex}) { return undef }
	
	$self->_ReplaceAtomLabels ($args{AtomIndex}, $args{Content}, "CHARMM");
	
	return 1;
} # of sub _PDB2CHARMM


sub _CHARMM2PDB { # converts CHARMM to PDB
	my $self = shift;
	my %args = @_;
	
	if (not defined $self->{Models}) { $self->Parse }
	if (not defined $args{Content} or not defined $args{AtomIndex}) { return undef }
	
	$self->_ReplaceAtomLabels ($args{AtomIndex}, $args{Content}, "PDB");
	
	return 1;
} # of sub _CHARMM2PDB


sub _ReplaceAtomLabels { # converts the atom labels from or to CHARMM format
	my $self      = shift;
	my $AtomIndex = shift;
	my $Content   = shift;
	my $Direction = shift;  # this is either PDB to convert to PDB or CHARMM to convert to CHARMM
	
	my ($LineCount, $TERCount, $FirstResidue, $LastResidue);
	
	$TERCount  = 0;
	$LineCount = 0;
	
	until ($LineCount > $#{$AtomIndex}) {
		# if converting to CHARMM, remove all lines not starting with the following tags
		if ($Direction eq "CHARMM" and $AtomIndex->[$LineCount]{Race} !~ m/HEADER|COMPND|ATOM|HETATM/) {
			splice (@{$AtomIndex}, $LineCount, 1);
			splice (@{$Content},   $LineCount, 1);
			next;
		}
		
		# count the TERs to trigger a warning for multiple chain files
		if ($AtomIndex->[$LineCount]{Race} eq "TER") {
			++$TERCount;
			++$LineCount;
			next;
		}
		
		if ($AtomIndex->[$LineCount]{Race} =~ m/ATOM|HETATM/) {
			
			# those replacements were taken from http://www.bmrb.wisc.edu/ref_info/atom_nom.tbl
			
			if ($AtomIndex->[$LineCount]{ResidueLabel} eq "ALA") {
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HB ", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HB ", " HB2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "3HB ", " HB3");
			}
			
			if ($AtomIndex->[$LineCount]{ResidueLabel} eq "ARG") {
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HB ", " HB2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HB ", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HG ", " HG2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HG ", " HG1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HD ", " HD2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HD ", " HD1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HH1", "HH11");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HH1", "HH12");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HH2", "HH21");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HH2", "HH22");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "HB3 ", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "HG3 ", " HG1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "HD3 ", " HD1");
			}
			
			if ($AtomIndex->[$LineCount]{ResidueLabel} eq "ASP") {
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HB ", " HB2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HB ", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "HB3 ", " HB1");
			}
			
			if ($AtomIndex->[$LineCount]{ResidueLabel} eq "ASN") {
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HB ", " HB2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HB ", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HD2", "HD21");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HD2", "HD22");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "HB3 ", "HB1" );
			}
			
			if ($AtomIndex->[$LineCount]{ResidueLabel} eq "CYS") {
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HB ", " HB2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HB ", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "HB3 ", " HB1");
			}
			
			if ($AtomIndex->[$LineCount]{ResidueLabel} eq "GLU") {
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HB ", " HB2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HB ", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HG ", " HG2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HG ", " HG1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " HB3", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " HG3", " HG1");
			}
			
			if ($AtomIndex->[$LineCount]{ResidueLabel} eq "GLN") {
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HB ", " HB2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HB ", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HG ", " HG2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HG ", " HG1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HE2", "HE21");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HE2", "HE22");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " HB3", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " HG3", " HG1");
			}
			
			if ($AtomIndex->[$LineCount]{ResidueLabel} eq "GLY") {
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HA ", " HA2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HA ", " HA1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " HA3", " HA1");
			}
			
			if ($AtomIndex->[$LineCount]{ResidueLabel} eq "HIS") {
				$Content->[$LineCount] =~ s/HIS/HSD/;
				$AtomIndex->[$LineCount]{AtomType} = "HSD";
				
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HB ", " HB2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HB ", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "HB3 ", " HB1");
			}
			
			if ($AtomIndex->[$LineCount]{ResidueLabel} eq "ILE") {
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HG1", "HG12");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HG1", "HG11");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HG2", "HG21");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HG2", "HG22");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HG2", "HG22");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "3HG2", "HG23");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HD1", "HD11");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HD1", "HD12");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "3HD1", "HD13");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " CD1", " CD ");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "HG13", "HG11");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "HD11", " HD1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "HD12", " HD2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "HD13", " HD3");
			}
			
			if ($AtomIndex->[$LineCount]{ResidueLabel} eq "LEU") {
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HB ", " HB2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HB ", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HD1", "HD11");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HD1", "HD12");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "3HD1", "HD13");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HD2", "HD21");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HD2", "HD22");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "3HD2", "HD23");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " HB3", " HB1");
			}
			
			if ($AtomIndex->[$LineCount]{ResidueLabel} eq "LYS") {
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HB ", " HB2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HB ", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HG ", " HG2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HG ", " HG1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HD ", " HD2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HD ", " HD1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HE ", " HE2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HE ", " HE1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HZ ", " HZ1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HZ ", " HZ2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "3HZ ", " HZ3");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " HB3", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " HG3", " HG1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " HD3", " HD1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " HE3", " HE1");
			}
			
			if ($AtomIndex->[$LineCount]{ResidueLabel} eq "MET") {
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HB ", " HB2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HB ", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HG ", " HG2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HG ", " HG1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HE ", " HE1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HE ", " HE2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "3HE ", " HE3");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " HB3", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " HG3", " HG1");
			}
			
			if ($AtomIndex->[$LineCount]{ResidueLabel} eq "PHE") {
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HB ", " HB2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HB ", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " HB3", " HB1");
			}
			
			if ($AtomIndex->[$LineCount]{ResidueLabel} eq "PRO") {
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HB ", " HB2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HB ", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HG ", " HG2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HG ", " HG1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HD ", " HD2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HD ", " HD1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " HB3", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " HG3", " HG1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " HD3", " HD1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " H1 ", " HT1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " H2 ", " HT2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " H3 ", " HT3");
			}
			
			if ($AtomIndex->[$LineCount]{ResidueLabel} eq "SER") {
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HB ", " HB2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HB ", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " HG ", " HG1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " HB3", " HB1");
			}
			
			if ($AtomIndex->[$LineCount]{ResidueLabel} eq "THR") {
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HG2", "HG21");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HG2", "HG22");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "3HG2", "HG23");
			}
			
			if ($AtomIndex->[$LineCount]{ResidueLabel} eq "TRP") {
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HB ", " HB2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HB ", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "HB3 ", " HB1");
			}
			
			if ($AtomIndex->[$LineCount]{ResidueLabel} eq "TYR") {
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HB ", " HB2");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HB ", " HB1");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " HB3", " HB1");
			}
			
			if ($AtomIndex->[$LineCount]{ResidueLabel} eq "VAL") {
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HG1", "HG11");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HG1", "HG12");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "3HG1", "HG13");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1HG2", "HG21");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2HG2", "HG22");
				$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "3HG2", "HG23");
			}
			
			$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " H  ", " HN ");
		}
		
		++$LineCount;
	}
	
	# if the last line is not TER, add it
	if ($Direction eq "CHARMM" and $AtomIndex->[ $#{$AtomIndex} ]{Race} ne "TER") {
		push @{$AtomIndex}, {Race => "TER"};
		push @{$Content}, "TER                                                                            \n";
	}
	
	$LineCount = 0;                                                      # start at the first line
	until ($AtomIndex->[$LineCount]{Race} eq "ATOM") { ++$LineCount }    # forward till the first ATOM line
	$FirstResidue = $AtomIndex->[$LineCount]{ResidueNumber};             # read the ResidueNumber
	
	# and replace the terminal hydrogen labels
	until ($AtomIndex->[$LineCount]{ResidueNumber} ne $FirstResidue) {
		$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "1H  ", " HT1");
		$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "2H  ", " HT2");
		$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, "3H  ", " HT3");
		
		++$LineCount;
	}
	
	$LineCount = $#{$AtomIndex};                                         # start the the end
	until ($AtomIndex->[$LineCount]{Race} eq "ATOM") { --$LineCount }    # rewind till the last ATOM line
	$LastResidue = $AtomIndex->[$LineCount]{ResidueNumber};              # read the ResidueNumber
	
	# now go backwards until the residue number changes (2nd to last residue)
	# and replace the terminal oxygen labels
	until ($AtomIndex->[$LineCount]{ResidueNumber} ne $LastResidue) {
		$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " OXT", " OT1");
		$self->_ReplaceLabel ($AtomIndex->[$LineCount], \$Content->[$LineCount], $Direction, " O  ", " OT2");
		
		--$LineCount;
	}

	$self->RenumberAtoms    (Content => $Content, AtomIndex => $AtomIndex);
	$self->RenumberResidues (Content => $Content, AtomIndex => $AtomIndex);
	
	if ($Direction eq "CHARMM" and $TERCount > 1) {
		$self->_ReportWarning ("PDB2CHARMM", "Misc", "More than one chain was converted into CHARMM format! " .
		                       "CHARMM can only process one chain at a time, so use ->WriteChains instead");
	}
} # of sub _ReplaceAtomLabels


sub _ReplaceLabel { # replaces an atom label in a given atom index line and a content line
	my $self        = shift;
	my $AtomIndex   = shift;
	my $Content     = shift;
	my $Direction   = shift;
	my $PDBLabel    = shift;
	my $CHARMMLabel = shift;
	
	my ($PDBType, $CHARMMType);
	
	# copy the atom labels and remove the blanks for the comparison
	$PDBType = $PDBLabel;
	$PDBType =~ s/\s//g;
	$CHARMMType = $CHARMMLabel;
	$CHARMMType =~ s/\s//g;
	
	if ($Direction eq "PDB") {
		if ($AtomIndex->{AtomType} eq $CHARMMType) {
			$AtomIndex->{AtomType} = $PDBType;
			substr (${$Content}, 12, 4) = $PDBLabel;
		}
	}
	else {
		if ($AtomIndex->{AtomType} eq $PDBType) {
			$AtomIndex->{AtomType} = $CHARMMLabel;
			substr (${$Content}, 12, 4) = $CHARMMLabel;
		}
	}
} # of sub _ReplaceLabel


sub _GetBaseName { # returns the base name of a given file, without path and extension
	my $self = shift;
	my ($FileName, $BaseName, $DotPosition, $Extension);
	
	$FileName = shift;
	if (not $FileName) { return undef }
	
	$Extension = $self->_GetExtension ($FileName); # get the extension
	
	if (defined $Extension) { # if an extension has been found
		$DotPosition = rindex ($FileName, "."); # determine the position of the last dot
		$BaseName = substr ($FileName, 0, $DotPosition); # take everything till the last dot
		return $BaseName;
	}
	else { # if no extension has been found, it is assumed that no one was there :-)
		return $FileName;
	}
} # of sub _GetBaseName


sub _GetExtension { # returns the extension (without the dot!) of the file name, just to be not forced to use .pdb
	my $self = shift;
	my $FileName = shift;
	my ($DotPosition, $Extension);
	
	# separate path, file name and extension
	$DotPosition = rindex ($FileName, "."); # determine the position of the last dot
	
	if ($DotPosition != -1) { # if a dot has been found
		$Extension = substr ($FileName, $DotPosition+1); # take everything after the last dot til the end
		return $Extension;
	}
	else { # if no dot has been found
		return undef;
	}
} # of sub _GetExtension


sub _HeaderRemark { # adds a remark to the header that the file may have been altered by the parser
	my $self = shift;
	my %args = @_;
	my ($Line, $Content, $AtomIndex, @RemarkIndex, $Position, $Race, $Rest);
	
	my @weekdays = qw(Sunday Monday Tuesday Wednesday Thursday Friday Saturday);
	my @HeaderRemark = ();
	my ($Date, $Time);
	my ($sec, $min, $hour, $mday, $mon, $year,  $wday) = localtime(time);
	
	$year = $year + 1900;
	$mon = $mon + 1;
	
	foreach ($mday, $mon, $hour, $min){ # add a leading zero if needed
		if (length($_) == 1) { $_ = "0" . $_ }
	}
	
	# Mind that the lines should be shorter than 80 characters!
	$HeaderRemark[0] = "REMARK This file was created by ParsePDB.pm\n";
	$HeaderRemark[1] = "REMARK Date: $weekdays[$wday], $mday/$mon/$year - $hour:$min\n";
	$HeaderRemark[2] = "REMARK Some data may have been altered and header and footer may not be valid!\n";
	
	for $Line (0 .. $#HeaderRemark) {
		$Race = substr ($HeaderRemark[$Line], 0, 6);
		$Race =~ s/\s//g;
		$Rest = substr ($HeaderRemark[$Line], 6);
		
		$RemarkIndex[$Line] = { Race => $Race, Rest => $Rest };
	}
		
	if (not defined $args{Content} and not defined $args{AtomIndex}) {
		return @HeaderRemark;
	}
	else {
		$Content = $args{Content};
		$AtomIndex = $args{AtomIndex};
		$Position = 0;
		
		for $Line (0 .. $#{$Content}) {
			if ($Content->[$Line] =~ m/^HEADER/) { $Position = $Line+1 }
			if ($Content->[$Line] =~ m/^COMPND/) { $Position = $Line+1 }
			
			# this determines the line BEFORE which the HeaderRemark will be added
			if ($Content->[$Line] =~ m/^$HeaderRemarkLine/) { $Position = $Line+1 };
			if ($Content->[$Line] =~ m/^(MODEL|ATOM|HETATM)/) { last }
		}
		
		splice (@{$Content}, $Position, 0, @HeaderRemark);
		splice (@{$AtomIndex}, $Position, 0, @RemarkIndex);
	}
} # sub _HeaderRemark


sub _InArray { # returns true, if the string is already in the array
	my $self = shift;
	my $RefArray = shift;
	my $String = shift;
	my $i;
	
	for $i (0 .. $#{$RefArray}) {
		if ($String eq $RefArray->[$i]) { return 1 } # if the string is found, return true
	}
	
	return 0; # else return false
} # of sub _InArray


sub _GetNumericKeys { # returns only the keys from a given hash that are numeric and omits words
	my $self = shift;
	my $Hash = shift;
	my $KeyCount = 0;
	my @Keys = keys (%{$Hash});
	
	until ($KeyCount > $#Keys) {
		# if the key does not consist entirely of digits, delete it
		if ($Keys[$KeyCount] !~ m/\d/g) { splice (@Keys, $KeyCount, 1) }
		                           else { ++$KeyCount                  }
	}
	
	return sort _numerically @Keys;
} # of sub _GetNumericKeys


sub _numerically { # used with sort, to sort a list of numbers numerically
	$a <=> $b;
} # of sub _numerically


sub _ReportWarning { # updates the warning hash and prints the warning message if verbose mode is active
	my $self    = shift;
	my $Method  = shift;
	my $Warning = shift;
	my $Var     = shift;
	
	my $WarningMsg;
	
	if ($self->{ReportNoWarnings}) { return 1 }
	
	if (not defined $Var) { $Var = " " } # to avoid a "Use of uninitialized value in concatenation"-error
	
	if ($Warning eq "Misc") { $WarningMsg = "Warning by $Method: " . $Var . "\n" }
	                   else { $WarningMsg = "Warning by $Method: " . $WARNINGS{$Warning} . $Var . "\n" }
	
	if ($self->{Verbose}) { print $WarningMsg }
	
	$self->{Warning}{Warning} = 1;  # any warning was reported
	$self->{Warning}{$Warning} = 1; # set the value of the respective warning true
	push @{$self->{WarningMsg}}, $WarningMsg; # save the warning message
	
	return 1;
} # of sub _ReportWarning


sub _ReportError { # updates the error hash and prints the error message if verbose mode is active
	my $self = shift;
	my $Error = shift;
	my $Var = shift;
	my $ErrorMsg;
	
	if (not defined $Var) { $Var = " " } # to avoid a "Use of uninitialized value in concatenation"-error
	
	$ErrorMsg = $ERRORS{$Error} . $Var . "\n";
	$self->{ErrorMsg} = $ErrorMsg; # save the error message
	
	if    (  ($Error eq "NoFile") or
	         ($Error eq "FileNotFound") or
	         ($Error eq "ParamNumeric") or 
	         ($Error eq "BadParameter") ) {
	         throw Exception::Config ("$ErrorMsg");
	}
	elsif (  ($Error eq "CorruptFile") or
	         ($Error eq "UnknownElement") ) {
	         throw Exception::PDB ("$ErrorMsg");
	}
	elsif ($Error eq "IOError") {
	         throw Exception::IO ("$ErrorMsg");
	}
	else {
		throw Exception::PDB ("$Error");
	}
} # of sub _ReportError


sub dp { # DebugPrint
	my @args = @_;
	my ($Var, $File);
	
	if ($args[$#args] eq "DumpIntoFile") { # if it was called by dpf
		open (OUTPUT, ">debug.txt");        # dump variable into file
		pop @args;                          # remove last element
		$File = 1;
	}
	else {
		open (OUTPUT, ">&STDOUT");    # dump variable to STDOUT
		$File = 0;
	}
	
	while (@args) {
		$Var = shift @args;
		
		print OUTPUT "\n\n";
		
		if (ref $Var) {
			if		($Var =~ m/ARRAY/)  {
				chomp @{$Var};
				print OUTPUT Dumper @{$Var};
			}
			elsif ($Var =~ m/HASH/)   { print OUTPUT Dumper %{$Var} }
			elsif ($Var =~ m/SCALAR/) { print OUTPUT Dumper ${$Var} }
			else  { print "Error determining reference type!\n" }
		}
		else {
			print OUTPUT Dumper $Var;
		}
		
		print OUTPUT "\n";
	
	} # of while (@args)
	
	close OUTPUT;
	
	if ($File) { die "\nDied happily a controlled death in debug print. Output saved to debug.txt\n\n" }
			else { die "\nDied happily a controlled death in debug print.\n\n" }
} # of sub dp


sub dpf { # DebugPrintFile, prints to debug.txt
	my @args = @_;
	
	push @args, "DumpIntoFile";
	dp (@args);
} # of sub dpf


1;

@Exception::Config::ISA = qw( Error::Simple ); # define @Exception::Config
@Exception::IO::ISA     = qw( Error::Simple ); # define @Exception::IO
@Exception::PDB::ISA    = qw( Error::Simple ); # define @Exception::PDB


