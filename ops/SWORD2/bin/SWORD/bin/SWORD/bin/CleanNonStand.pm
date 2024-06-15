sub CleanNonStand
{

	my $pathpdbs  = shift;
	my $PDB  = shift;
	my $chain = '';
	
	$chain = substr($PDB,4,1);
	$PDB 	= substr($PDB,0,4);
    $PDB = $PDB.$chain;

	my $temp = '';
	$temp = `mkdir -p $pathpdbs/PDBs_Stand/`; # directory that will contain "standardized" PDB files

	my %corres_aa_nonstand=(
	 "  A" => "ALA","  B" => "UNK","  C" => "CYS","  D" => "ASP","  E" =>  "GLU",
	 "  F" => "PHE","  G" => "GLY","  H" => "HIS","  I" => "ILE","  J" =>  "UNK",
	 "  K" => "LYS","  L" => "LEU","  M" => "MET","  N" => "ASN","  O" =>  "LYS",
	 "  P" => "PRO","  Q" => "GLN","  R" => "ARG","  S" => "SER","  T" =>  "THR",
	 "  U" => "CYS","  V" => "VAL","  W" => "TRP","  X" => "UNK","  Y" =>  "TYR",

	 "OCS" => "SER"," DA" => "ALA"," DC" => "CYS","CSD" => "CYS","1AC" => "UNK",
	 "HYP" => "PRO","LVN" => "UNK","UNK" => "UNK","FOX" => "UNK","KPI" => "UNK",
	 "5IU" => "UNK","KCX" => "LYS"," DG" => "GLY"," DT" => "TYR","  Z" => "UNK",
	 "8OG" => "UNK","SNN" => "UNK","ACY" => "UNK","TRW" => "UNK","BRU" => "UNK",
	 "CSE" => "UNK","B3Q" => "UNK","B3L" => "UNK","BIL" => "UNK","PRN" => "UNK",
	 "CSS" => "UNK","DDZ" => "UNK","LP6" => "UNK","FHU" => "UNK","ASX" => "UNK",
	 "5CM" => "UNK","A1P" => "UNK","TTD" => "UNK","5BU" => "UNK","DHA" => "UNK",
	 "6MA" => "UNK","5PY" => "UNK","MIR" => "UNK","6OG" => "UNK","GL3" => "UNK",
	 "MHS" => "UNK","AGM" => "UNK","SMC" => "UNK","MEN" => "UNK","2HF" => "UNK",
	 "K1R" => "CYS","MEQ" => "CYS","PED" => "CYS","O12" => "CYS","CAV" => "UNK",
	 "APY" => "UNK",
	 
	 "ARG" => "ARG","ASN" => "ASN","ASP" => "ASP","ALA" => "ALA","CYS" => "CYS",
	 "GLU" => "GLU","GLN" => "GLN","GLY" => "GLY","HIS" => "HIS","ILE" => "ILE",
	 "LEU" => "LEU","LYS" => "LYS","MET" => "MET","PHE" => "PHE","PRO" => "PRO",
	 "SER" => "SER","THR" => "THR","TRP" => "TRP","TYR" => "TYR","VAL" => "VAL",
	 "GLY" => "GLY", "PAQ" => "TYR", "B2V" => "VAL", "SIC" => "ASP", "B2I" => "ILE",
	 "B2F" => "PHE", "AGT" => "CYS", "KYN" => "ALA", "R1A" => "CYS", "0CS" => "ALA", 
	 "NYS" => "CYS", "SCY" => "CYS", "KYQ" => "LYS", "VAF" => "VAL", "VAD" => "ALA", 
	 "AHB" => "ARG", "SNC" => "CYS", "OSE" => "SER", "AHO" => "ALA", "MAA" => "ALA", 
	 "AHP" => "ALA", "NLO" => "LEU", "XYG" => "ASP", "MAI" => "ARG", "DTH" => "THR", 
	 "CMH" => "CYS", "PCS" => "PHE", "CME" => "CYS", "1TY" => "TYR", "DTR" => "TRP", 
	 "CMT" => "CYS", "B3Y" => "TYR", "VLL" => "VAL", "TTS" => "TYR", "LSO" => "LYS", 
	 "TTQ" => "TRP", "CHG" => "ALA", "NYG" => "ASN", "M2L" => "LYS", "FTR" => "TRP", 
	 "BLY" => "LYS", "NFA" => "PHE", "FTY" => "TYR", "P1L" => "CYS", "DSN" => "SER", 
	 "DSG" => "ASN", "IGL" => "GLY", "HLU" => "LEU", "HSE" => "SER", "LME" => "GLU", 
	 "GLH" => "GLN", "5CS" => "CYS", "B3D" => "ASP", "GLZ" => "GLY", "AKL" => "ASP", 
	 "SMF" => "PHE", "SME" => "MET", "CTH" => "THR", "DLY" => "LYS", "CS3" => "CYS", 
	 "MDO" => "ALA", "B3A" => "ALA", "FCL" => "PHE", "R7A" => "CYS", "MVA" => "VAL", 
	 "ALN" => "ALA", "ALO" => "THR", "ALM" => "ALA", "ALG" => "ARG", "ALC" => "ALA", 
	 "ORQ" => "ARG", "ORN" => "ALA", "SBG" => "SER", "SBD" => "SER", "ALY" => "LYS", 
	 "ALT" => "ALA", "GYS" => "SER", "ALS" => "ALA", "SBL" => "SER", "IAM" => "ALA", 
	 "TOX" => "TRP", "DPR" => "PRO", "DPQ" => "TYR", "DPP" => "ALA", "CAS" => "CYS", 
	 "NC1" => "SER", "B2A" => "ALA", "CAB" => "ALA", "DPN" => "PHE", "BAL" => "ALA", 
	 "DPL" => "PRO", "CAF" => "CYS", "HIP" => "HIS", "HIQ" => "HIS", "KOR" => "CYS", 
	 "4IN" => "TRP", "DYS" => "CYS", "HIA" => "HIS", "HIC" => "HIS", "AB7" => "GLU", 
	 "C5C" => "CYS", "ABA" => "ALA", "143" => "CYS", "PG1" => "SER", "NZH" => "HIS", 
	 "AZK" => "LYS", "PG9" => "GLY", "3MD" => "ASP", "23F" => "PHE", "BHD" => "ASP", 
	 "1AB" => "PRO", "QLG" => "GLN", "NCB" => "ALA", "A0A" => "ASP", "23S" => "TRP", 
	 "GHP" => "GLY", "GHG" => "GLN", "SLZ" => "LYS", "MIS" => "SER", "TY2" => "TYR", 
	 "TY3" => "TYR", "SAH" => "CYS", "SAC" => "SER", "FGP" => "SER", "4PH" => "PHE", 
	 "FGL" => "GLY", "PRQ" => "PHE", "PRS" => "PRO", "CH7" => "LYS", "SAR" => "GLY", 
	 "4F3" => "GLY", "LAL" => "ALA", "BFD" => "ASP", "TPO" => "THR", "H5M" => "PRO", 
	 "AYG" => "ALA", "TPL" => "TRP", "AYA" => "ALA", "4FB" => "PRO", "CHP" => "GLY", "CHS" => "PHE", "TPQ" => "TYR", "IIL" => "ILE", "4FW" => "TRP", "LPD" => "PRO", "NPH" => "CYS", "PYX" => "CYS", "TYN" => "TYR", "TYI" => "TYR", "PYR" => "SER", "PYT" => "ALA", "TYB" => "TYR", "MBQ" => "TYR", "DBY" => "TYR", "APK" => "LYS", "API" => "LYS", "TYX" => "CYS", "TYY" => "TYR", "TYZ" => "ARG", "TYT" => "TYR", "PYA" => "ALA", "PYC" => "PRO", "TYQ" => "TYR", "TYS" => "TYR", "KST" => "LYS", "CZ2" => "CYS", "OCY" => "CYS", "P2Y" => "PRO", "HMR" => "ARG", "NLE" => "LEU", "BMT" => "THR", "CS4" => "CYS", "HMF" => "ALA", "IT1" => "LYS", "CS1" => "CYS", "NLN" => "LEU", "SHP" => "GLY", "CSO" => "CYS", "CSI" => "GLY", "DLE" => "LEU", "CSA" => "CYS", "AFA" => "ASN", "CSB" => "CYS", "DLS" => "LYS", "SHC" => "CYS", "CSY" => "THR", "CSX" => "CYS", "CSZ" => "CYS", "CSU" => "CYS", "B3E" => "GLU", "CSW" => "CYS", "CSP" => "CYS", "CSR" => "CYS", "CZZ" => "CYS", "N10" => "SER", "MGG" => "ARG", "BTR" => "TRP", "MGN" => "GLN", "C1X" => "LYS", "TBM" => "THR", "TBG" => "GLY", "CFY" => "PHE", "CML" => "CYS", "175" => "ALA", "2TY" => "TYR", "CRW" => "ALA", "SET" => "SER", "C12" => "GLY", "SEP" => "SER", "LNT" => "LEU", "1TQ" => "TRP", "DGN" => "GLN", "DGL" => "GLU", "ASB" => "ASP", "ASA" => "ASP", "SEC" => "ALA", "SEB" => "SER", "ASK" => "ASP", "GT9" => "CYS", "ASI" => "ASP", "SEL" => "SER", "ASL" => "ASP", "CLH" => "LYS", "FRF" => "PHE", "CLD" => "ALA", "CLE" => "LEU", "LA2" => "LYS", "BBC" => "CYS", "IML" => "ILE", "CLB" => "ALA", "OXX" => "ASP", "CLV" => "ALA", "VME" => "VAL", "DTY" => "TYR", "CY4" => "CYS", "OMT" => "MET", "2MR" => "ARG", "BPE" => "CYS", "2MT" => "PRO", "AAR" => "ARG", "MC1" => "SER", "CYR" => "CYS", "CYQ" => "CYS", "GYC" => "CYS", "BIF" => "PHE", "POM" => "PRO", "CYG" => "CYS", "CYF" => "CYS", "CYD" => "CYS", "CYJ" => "LYS", "CY3" => "CYS", "CY1" => "CYS", "CLG" => "LYS", "CY7" => "CYS", "DHI" => "HIS", "MLE" => "LEU", "CWR" => "ALA", "32S" => "TRP", "32T" => "TRP", "MLZ" => "LYS", "MLY" => "LYS", "MCL" => "LYS", "AA4" => "ALA", "NRQ" => "MET", "MCS" => "CYS", "PAT" => "TRP", "TH5" => "THR", "4HT" => "TRP", "HSO" => "HIS", "HSL" => "SER", "DCY" => "CYS", "GPL" => "LYS", "BNN" => "ALA", "THC" => "THR", "IAS" => "ASP", "PHI" => "PHE", "PHM" => "PHE", "PHL" => "PHE", "PHA" => "PHE", "PHD" => "ASP", "SYS" => "CYS", "BG1" => "SER", "NIY" => "TYR", "OAS" => "SER", "CRQ" => "GLN", "CRU" => "GLU", "DMH" => "ASN", "PRR" => "ALA", "CRX" => "ALA", "DMK" => "ASP", "AEI" => "THR", "CRK" => "MET", "CRO" => "GLY", "MFC" => "GLY", "B3X" => "ASN", "LDH" => "LYS", "TCQ" => "TYR", "BUC" => "CYS", "BUG" => "LEU", "MHO" => "MET", "MHL" => "LEU", "FOE" => "CYS", "FOG" => "PHE", "EFC" => "CYS", "DDE" => "HIS", "OHI" => "HIS", "CR2" => "GLY", "CR0" => "THR", "CR7" => "LYS", "OPR" => "ARG", "CR5" => "GLY", "CR8" => "HIS", "PSA" => "PHE", "LCX" => "LYS", "NVA" => "VAL", "LCK" => "LYS", "ILG" => "GLU", "6CL" => "LYS", "DVA" => "VAL", "BLE" => "LEU", "BCX" => "CYS", "ILX" => "ILE", "BCS" => "CYS", "DSE" => "SER", "OLT" => "THR", "LED" => "LEU", "C6C" => "CYS", "IEY" => "HIS", "CH6" => "MET", "6CW" => "TRP", "PSH" => "HIS", "B3K" => "LYS", "S1H" => "SER", "FZN" => "LYS", "B3S" => "SER", "OBS" => "LYS", "DYG" => "ASP", "PLE" => "LEU", "NMC" => "GLY", "NMM" => "ARG", "DIL" => "ILE", "FHL" => "LYS", "AR2" => "ARG", "AR4" => "GLU", "RCY" => "CYS", "CY0" => "CYS", "IYR" => "TYR", "DIV" => "VAL", "AIB" => "ALA", "SOC" => "CYS", "XXY" => "THR", "PPH" => "LEU", "PPN" => "PHE", "2ML" => "LEU", "XX1" => "LYS", "MTY" => "TYR", "YOF" => "TYR", "HPQ" => "PHE", "M0H" => "CYS", "C99" => "THR", "HPH" => "PHE", "4DP" => "TRP", "HPE" => "PHE", "SDP" => "SER", "RC7" => "HIS", "ARM" => "ARG", "ARO" => "ARG", "PR3" => "CYS", "BOR" => "ARG", "SVY" => "SER", "CCS" => "CYS", "NAM" => "ALA", "TIH" => "ALA", "NAL" => "ALA", "HRG" => "ARG", "LYP" => "LYS", "SXE" => "SER", "LYR" => "LYS", "ML3" => "LYS", "PIA" => "ALA", "LYX" => "LYS", "LYZ" => "LYS", "LYM" => "LYS", "LYN" => "LYS", "OHS" => "ASP", "SVA" => "SER", "TYO" => "TYR", "SVZ" => "SER", "SVX" => "SER", "SVV" => "SER", "CQR" => "GLY", "IPG" => "GLY", "3TY" => "TYR", "BIU" => "ILE", "TMD" => "THR", "CXM" => "MET", "SOY" => "SER", "FLT" => "TYR", "DHL" => "SER", "TYC" => "TYR", "FLE" => "LEU", "VDL" => "VAL", "PTM" => "TYR", "SCH" => "CYS", "PTH" => "TYR", "LLP" => "LYS", "AME" => "MET", "MNV" => "VAL", "DBZ" => "ALA", "MNL" => "LEU", "C3Y" => "CYS", "PEC" => "CYS", "PTR" => "TYR", "SCS" => "CYS", "TRO" => "TRP", "TRN" => "TRP", "R1F" => "CYS", "TRX" => "TRP", "PBF" => "PHE", "R1B" => "CYS", "MLL" => "LEU", "YCM" => "CYS", "PBB" => "CYS", "TRQ" => "TRP", "SC2" => "CYS", "HTI" => "CYS", "M3L" => "LYS", "HTR" => "TRP", "MPQ" => "GLY", "ACB" => "ASP", "SUN" => "SER", "SUI" => "ASP", "ACL" => "ARG", "NEP" => "HIS", "NBQ" => "TYR", "TEE" => "CYS", "HY3" => "PRO", "CGU" => "GLU", "DNE" => "LEU", "DNG" => "LEU", "UMA" => "ALA", "PM3" => "PHE", "DNM" => "LEU", "DNL" => "LYS", "GMA" => "GLU", "1LU" => "LEU", "MEG" => "GLU", "MEA" => "PHE", "LEF" => "LEU", "MEU" => "GLY", "LET" => "LYS", "X9Q" => "ALA", "3AH" => "HIS", "1PA" => "PHE", "MSO" => "MET", "CIR" => "ARG", "PF5" => "PHE", "MSE" => "MET", "MSA" => "GLY", "DAB" => "ALA", "OTY" => "TYR", "DAL" => "ALA", "SGR" => "SER", "DAH" => "PHE", "DAR" => "ARG", "DAS" => "ASP", "SGB" => "SER", "TRF" => "TRP", "GVL" => "SER", "LLY" => "LYS", "ESC" => "MET", "ESB" => "TYR", "DA2" => "ARG", "TNB" => "CYS", "PFF" => "PHE", "KGC" => "LYS", "HHK" => "ALA");

	my $line = '';
	my @lines = ();
	my $i = 0;
	my $j = 0;
	my $flag1 = 0;
	my $flag2 = 0;
	my $count = 0;
	my $counter1 = 0;
	my $counter2 = 0;

	my $residu = '';
	my $residu_stand = '';
	
	open (IN, '<', "$pathpdbs/$PDB") || die "Error while reading $pathpdbs/$PDB:\n\t$!\n";
	open (OUT, '>', "$pathpdbs/PDBs_Stand/$PDB") || die "Error while writing $pathpdbs/PDBs_Stand/$PDB:\n\t$!\n";

	while ($line = <IN>)
	{
		my @buff = split(' ',$line);
		if($buff[0] ne 'ANISOU')
		{
			push(@lines,$line);
		}
	}

	for ($i = 0; $i < $#lines; $i++)
	{
		my @buff = split(' ',$lines[$i]);
		my @buff2 = ();
		my @buff3 = ();

		if(($buff[0] ne 'HETATM') and ($buff[0] ne 'ATOM')){		
			$flag1 = 0;	
			$flag2 = 0;	
			$count = 0;					
		}
		
		if($buff[0] eq 'ATOM'){			
			$flag1 = 1;	
			$counter1 = (substr($lines[$i], 22, 4));

			$residu = substr($lines[$i],17,3);
			$residu_stand = $corres_aa_nonstand{$residu};
			
			@buff2 = split('',$lines[$i]);
			@buff3 = split('',$residu_stand);
			
			$buff2[17] = $buff3[0];
			$buff2[18] = $buff3[1];
			$buff2[19] = $buff3[2];
			$lines[$i] = join('',@buff2);

			if(($buff[0] eq 'ATOM') and ($flag2 == 1) and (($counter1-1) == ($counter2))){
				$flag2 = 0;	
				for ($j = 1 ; $j <= $count; $j++)
				{
					
					# replce hetatm and res3
					$residu = substr($lines[$i-$j],17,3);
					$residu_stand = $corres_aa_nonstand{$residu};
					
					@buff2 = split('',$lines[$i-$j]);
					@buff3 = split('',$residu_stand);
					
					$buff2[17] = $buff3[0];
					$buff2[18] = $buff3[1];
					$buff2[19] = $buff3[2];
					$buff2[0] = 'A';
					$buff2[1] = 'T';
					$buff2[2] = 'O';
					$buff2[3] = 'M';
					$buff2[4] = ' ';
					$buff2[5] = ' ';
					$lines[$i-$j] = join('',@buff2);
					
				}
			$count = 0;
			}
		}
		
		if(($buff[0] eq 'HETATM') and ($flag1 == 1) and ((substr($lines[$i], 22, 4)) == ($counter1+1))){	
			$flag1 = 0;									
			$flag2 = 1;	
			$counter2 = (substr($lines[$i], 22, 4));
		}
		
		if(($buff[0] eq 'HETATM') and ($flag2 == 1)){
			$count++;	
			$counter2 = (substr($lines[$i], 22, 4));								
		}
		
		
	}
		
	foreach (@lines)
	{
		print OUT $_;
	}
	close OUT;
	close IN;
}


1;
