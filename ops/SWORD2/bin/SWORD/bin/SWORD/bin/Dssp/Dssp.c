/* Output from p2c, the Pascal-to-C translator */
/* From input file "dssp.p" */

/* This is file DSSP.PAS and includes:

   1. DSSP header / license agreement (as PASCAL comment)
   2. DSSP source code
   3. sample  input file (as PASCAL comment),
      becomes logical file TAPEIN under VMS.
   4. sample output file (as PASCAL comment),
      becomes logical TAPEOUT under DEC/VMS).

Compile this file DSSP.PAS as is, using your favorite PASCAL compiler.

   DEC/VMS command file to run DSSP:

   $   assign myprotein.brk tapein    ! input coordinates
   $   assign myprotein.dssp tapeout  ! output
   $   run dssp                                                       */

/* ------------------------------------------------------------------ */
/*

DSSP version October 1988.
This copy for $name at $place
who have agreed to the following software license agreement:


An academic license for the DSSP program
((c) W. Kabsch, C. Sander and MPI-MF, 1983, 1985, 1988)
is granted to in exchange for the following commitments:

I hereby certify that

        (1) I am an academic user at an academic research institution. In
            using the software, I will respect the interests of the authors
            and their institutions.

        (2) I will not use the software in commercial activities without
            a written commercial license agreement; commercial activities
            include, in particular, work under contract from a commercial
            company.

        (3) I will not redistribute the software to others outside of my
            immediate research group. I will suggest to other interested
            research groups to contact the authors directly.

        (4) I will not alter or suppress the run-time copyright message.

        (5) I will acknowledge the program authors on any publication of
            scientific results based in part on use of the program and
            cite the article in which the program was described.

        (6) I will report evidence of program bugs to the authors.

        (7) I will send the source code of any bug corrections and program
            extensions, major or minor, to the original authors, for free
            academic use. If I have made major extensions which are incor-
            porated by the authors, I reserve the right to be appropriately
            included in any future commercial license agreement.

        (8) I will not extract part of the software, e.g. modules or sub-
            routines, for use in other contexts without permission by the
            authors.

        (9) I will not use the program in the context of classified research.
*/
/* PREVIOUS RELEASE: VERSION OCTOBER 1985                             */
/* PREVIOUS RELEASE: VERSION JUNE 1983                                */
/* LANGUAGE: STANDARD PASCAL WITH 128 CHARACTER ASCII SET             */
/* AUTHORS AND COPYRIGHT (1983,1985,1988,1993):
   Wolfgang Kabsch and Chris Sander, Max Planck Institut
   fuer Medizinische Forschung, Jahnstr. 29, 6900 Heidelberg, Germany
   Telephone: +49-6221-486 276  Telex: 461505 mpimf d
   Bitnet:    KABSCH@EMBL
   Current address for Chris Sander:
   Biocomputing, EMBL, 6900 Heidelberg, Germany
   Telephone: +49-6221-387 361 Telex: 461613 embl d
   Telefax:   +49-6221-387 306
   Bitnet:    SANDER@EMBL
   Do report errors if you find any.
   Reference: Kabsch,W. and Sander,C. (1983) Biopolymers 22, 2577-2637*/
/*--------------------------------------------------------------------*/
/* DEFINES SECONDARY STRUCTURE AND SOLVENT EXPOSURE OF PROTEINS FROM
   ATOMIC COORDINATES AS GIVEN BY THE BROOKHAVEN PROTEIN DATA BANK.   */
/*--------------------------------------------------------------------*/
/* This program including sample input and output files for dataset 1PPT
   is available from the authors in exchange for an academic or
   commercial license agreement. The program is no longer available
   from the Brookhaven Protein Data Bank */
/*--------------------------------------------------------------------*/
/* CORRECTION AND MODIFICATION LOG SINCE JUNE 1983 */
/* (1) MODIFICATIONS THAT AFFECT OUTPUT ON FILE TAPEOUT FOR AT LEAST ONE
       OF THE 62 PROTEIN DATA SETS IN THE 1983 BIOPOLYMERS PAPER:
   - SIDECHAIN ATOMS MORE THAN MAXDIST ANGSTROM DISTANT FROM ATOM CA ARE
     DECLARED ILLEGAL AND IGNORED. OUTPUT CHANGE: ACCESSIBILITY VALUES
     FOR ASN 76 OF 1SBT (ILLEGAL ATOM OD1) AND PRO 49 OF 156B (ILLEGAL
     ATOM UNK).
   - ANY RESIDUE WITH INCOMPLETE BACKBONE IS IGNORED. OUTPUT CHANGE:
     CHAIN BREAK BETWEEN RESIDUE SER 11 AND ILE 16 IN 2GCH
     (DUE TO INCOMPLETE COORDINATES FOR SER 11) IS NOW CORRECT.
   (2) MODIFICATIONS THAT DO NOT AFFECT OUTPUT ON FILE TAPEOUT FOR ANY
       OF THE 62 PROTEIN DATA SETS IN THE 1983 BIOPOLYMERS PAPER:
   - SPELLING OF FLAGCHIRALITY AND TESTSSBOND CORRECTED.
   - WARNING MESSAGE FOR RESIDUES WITH NON-STANDARD NUMBER OF
     SIDECHAIN ATOMS HAS BEEN ADDED. FOR EXAMPLE, THIS ALERTS THE USER
     TO BAD DATA FOR RESIDUES 8,12,21,24 AND 44 OF DATA SET 2RXN.
   - WARNING MESSAGE FOR RESIDUES IGNORED DUE TO NON-STANDARD RESIDUE
     NAME SUCH AS 'ACE' AND 'FOR' HAS BEEN ADDED.
   - WARNING MESSAGE FOR ALTERNATE ATOM LOCATION IDENTIFIER HAS BEEN
     ADDED. FOR EXAMPLE, THE USER IS NOW WARNED THAT ATOM CH2 IN ANY
     TRP OF DATA SET 1APP IS IGNORED DUE TO BAD ALTERNATE LOCATION
     IDENTIFIER '2'.
   WE THANK STEVEN SHERIFF, FRANCOIS COLONNA AND JOHN MOULT FOR
   REPORTING PROBLEMS AND STEVEN SHERIF, JANET THORNTON AND
   WILLIE TAYLOR FOR RESULTS OF TEST RUNS ON VAX COMPUTERS.

   Changes after 1985:

   - program speeded up by a factor of two or three by avoiding use
     of square root.
   - hydrogen atoms in data set ignored on input (H of NH of backbone
     is built as before)
   - 18-AUG-1988: CADIST=9.0, replacing CADIST=8.0. Has affected output
     for 63/300 proteins in a minor way. Thanks to Jean Richelle (Bruxelles)
     for pointing out this bug.

     Output changes due to change in parameter CADIST (8 to 9 Angstrom) :
     additional backbone-backbone Hbonds found with slight
     adjustments in secondary structure summary. In about 300 protein
     data sets from the Fall 1988 PDB release, 63 additional
     Hbonds were found, i.e. 0.2 Hbonds per protein (29 of type
     i,i+5;  16 of type i,i+4; 6 of type i,i+3; 10 in antiparallel beta
     bridges and 2 in a parallel beta bridge). These additional
     Hbonds affected the secondary structure summary of 26 of these
     protein data sets in a minor way, typically joining a 4-turn to
     an alpha-helix, changing a geometrical turn to a hydrogen-
     bonded turn or adding an extra residue pair to a beta ladder.
     The changes are (using _ for blank):

     [protein id, old secstruc summary > corrected summary]

     1FC2       _E > EE   and  _T > ET
     1GP1       GGG > HHH
     1HBS       S > T
     1HDS       S > T and  GGGGGG > TTHHHH
     1HFM       __ > TT
     1HKG       SSS > TTT
     1IG2       S_ > TT
     1LDX        GGG > HTT
     1MEV       __ > TT  and  _BS > TBS  and  SSS > TTS
     1PFC       SSS > TTS
     1PP2       _E > EE  and  _S > ES
     1RN3       E_SS_E > EEEEEE  and _E > EE  (>3-res beta bulge)
     1RNS       same as 1RN3
     2ATC       HH > TT
     2CAB       B_ > EE
     2CPP       SS > TT  and  GGGGGG > HHHHTT
     2LYZ       T > H
     2MDH       SSS > TTT
     3CPA       TTT > HHH
     4CAT       TTT > HHH
     4SBV       S > T
     5API       _ > B
     5CPA       TTT > HHH
     7LYZ       S > H
     8CAT       _ > B  and  _ > B
     8LYZ       T > H

     Note that this bugfix results in a small variation in the total
     number of Hbonds, compared to the variation which would
     result, say, from changing the (somewhat arbitrary) cutoff of
     -0.5 kcal/mol for the Hbond electrostatic potential energy. We
     cannot here solve the fundamental difficulty of arbitrary
     cutoffs involved in extracting binary classifications (an Hbond
     exists, yes/no) from real numbers (coordinates, energies).
     However, for most purposes the secondary structure summary agrees
     will with anyone's intuitive definition, especially for well-refined and
     high resolution structures. For a more clearcut assignment of protein
     substructure, we recommend using the detailed H-bond and other assignments
     in the columns following the summary column, i.e. columns 19-38 (xxx):

     ....;....1....;....2....;....3....;....4....;....5....;....6....;....7..
                       xxxxxxxxxxxxxxxxxxxx
                       .-- 3-turns/helix
                       |.-- 4-turns/helix
                       ||.-- 5-turns/helix
                       |||.-- geometrical bend
                       ||||.-- chirality
                       |||||.-- beta bridge label
                       ||||||.-- beta bridge label
                       |||||||   .-- beta bridge partner resnum
                       |||||||   |   .-- beta bridge partner resnum
                       |||||||   |   |.-- beta sheet label
                       |||||||   |   ||   .-- solvent accessibility
                       |||||||   |   ||   |
        35   47   I  E     +     0   0    2
        36   48   R  E >  S- K   0  39C  97
        37   49   Q  T 3  S+     0   0   86    (example from 1EST)
        38   50   N  T 3  S+     0   0   34
        39   51   W  E <   -KL  36  98C   6
                                                                           */
/*--------------------------------------------------------------------*/
/* GENERAL PROGRAM INSTALLATION GUIDE. */
/* (1) THE PROGRAM REQUIRES THE FULL STANDARD ASCII 128 CHARACTER SET,
       IN PARTICULAR LOWER CASE LETTERS 'abcdefg....'.
   (2) STANDARD PASCAL MAY NOT RECOGNIZE REAL NUMBERS SUCH AS .1, +.1,
       -.1 ON INPUT. CHANGE TO 0.1,+0.1,-0.1.
   (3) THE NON-STANDARD PROCEDURE 'DATE' RETURNS THE CURRENT DAY, MONTH,
       AND YEAR. IF THE PROCEDURE IS NOT CALLED (LINE COMMENTED OUT)
       THE PSYCHEDELIC DATE DEC 24, 2001 IS RETURNED. YOU MAY  REPLACE
       'DATE' BY THE CORRESPONDING PROCEDURE FROM YOUR PASCAL
       IMPLEMENTATION. THE EXAMPLE GIVEN WORKS IN DEC VAX VMS 5.0.
   (4) DUE TO INCOMPATIBLE ASCII CODES, SQUARE BRACKETS '[' AND ']'
       MAY APPEAR AS '!','?' ETC. USE YOUR EDITOR TO CONVERT THESE.   */
/* INSTALLATION GUIDE FOR VAX/VMS USERS. */
/* (1) THE /OPTIMIZE OPTION OF THE PASCAL COMPILER PRODUCED
       INCORRECT CODE ON THE VAX 8600 AT EMBL RUNNING UNDER VMS V4.2.
       LATER VERSIONS OF VMS (E.G. VMS 5.0) PRODUCED CORRECT CODE.
       IF IN DOUBT, COMPILE USING PASCAL /NOOPTIMIZE.
   (2) COPY BROOKHAVEN DATA BANK COORDINATE INPUT TO A FILE NAMED
       TAPEIN.DAT . OUTPUT WILL BE IN A FILE NAMED TAPEOUT.DAT        */
/* IMPLEMENTATION ON OTHER COMPUTERS */
/* (1) NORD-500. EXECUTION TIME COMPARABLE TO VAX 780.
   (2) SUN-3.    EXECUTION TIME COMPARABLE TO VAX 780.
                 Compile using: pc -L
                 in ORDER to map upper case letters in keywords
                 and identifiers to lower case.
   (3) ATARI 520 ST. RUNS FACTOR 60 SLOWER THAN NORD-500 DUE TO
       SOFTWARE-EMULATED FLOATING POINT OPERATIONS ON MC68000.        */
/*--------------------------------------------------------------------*/
/* INPUT/OUTPUT FILES. */
/* INPUT:   DEFAULT  INPUT UNIT, E.G. YOUR TERMINAL
   OUTPUT:  DEFAULT OUTPUT UNIT, E.G. YOUR TERMINAL,
            USED FOR RUN-TIME MESSAGES. WARNINGS AND ERRORS LOOK
            LIKE THIS: !!! TEXT !!!
   TAPEIN:  FILE WITH PROTEIN DATA BANK COORDINATES, E.G. PDB3PTI.COO
   TAPEOUT: DSSP OUTPUT OF LINE LENGTH 128, E.G. PAPER PRINTER        */
/*--------------------------------------------------------------------*/
/* DESCRIPTION OF OUTPUT ON FILE TAPEOUT:
   LINE LENGTH OF OUTPUT IS 128 CHARCTERS.
   FOR DEFINITONS, SEE ABOVE BIOPOLYMERS ARTICLE.
   IN ADDITION NOTE:
   HISTOGRAMS - E.G. 2 UNDER COLUMN '8' IN LINE 'RESIDUES PER ALPHA
            HELIX' MEANS: THERE ARE 2 ALPHA HELICES OF LENGTH  8
            RESIDUES IN THIS DATA SET.
   #  RESIDUE AA STRUCTURE BP1 BP2 ACC ..ETC..FOR EACH RESIDUE I:
   #  RESIDUE - TWO COLUMNS OF RESIDUE NUMBERS. FIRST COLUMN IS DSSP'S
            SEQUENTIAL RESIDUE NUMBER, STARTING AT THE FIRST
            RESIDUE ACTUALLY IN THE DATA SET AND INCLUDING CHAIN BREAKS;
            THIS NUMBER IS USED TO REFER TO RESIDUES THROUGHOUT. SECOND
            COLUMN GIVES CRYSTALLOGRAPHERS' 'RESIDUE SEQUENCE
            NUMBER','INSERTION CODE' AND 'CHAIN IDENTIFIER' (SEE PROTEIN
            DATA BANK FILE RECORD FORMAT MANUAL), GIVEN FOR REFERENCE
            ONLY AND NOT USED FURTHER..
   AA -     ONE LETTER AMINO ACID CODE, LOWER CASE FOR SS-BRIDGE CYS.
   STRUCTURE - SEE BIOPOLYMERS
   BP1 BP2  - RESIDUE NUMBER OF FIRST AND SECOND BRIDGE PARTNER
            FOLLOWED BY ONE LETTER SHEET LABEL
   ACC -    NUMBER OF WATER MOLECULES IN CONTACT WITH THIS RESIDUE *10.
            OR RESIDUE WATER EXPOSED SURFACE IN ANGSTROM**2.
   N-H-->O ETC. -  HYDROGEN BONDS. E.G. -3,-1.4 MEANS: IF THIS RESIDUE
            IS RESIDUE I THEN N-H OF I IS H-BONDED TO C=O OF I-3
            WITH AN ELECTROSTATIC H-BOND ENERGY OF -1.4 KCAL/MOL.
   TCO -    COSINE OF ANGLE BETWEEN C=O OF RESIDUE I AND C=O OF
            RESIDUE I-1. FOR ALPHA-HELICES, TCO IS NEAR +1, FOR
            BETA-SHEETS TCO IS NEAR -1. NOT USED FOR STRUCTURE
            DEFINITION.
   KAPPA -  VIRTUAL BOND ANGLE (BEND ANGLE) DEFINED BY THE THREE
            C-ALPHA ATOMS OF RESIDUES I-2,I,I+2. USED TO DEFINE
            BEND (STRUCTURE CODE 'S').
   ALPHA -  VIRTUAL TORSION ANGLE (DIHEDRAL ANGLE) DEFINED BY THE FOUR
            C-ALPHA ATOMS OF RESIDUES I-1,I,I+1,I+2. USED TO DEFINE
            CHIRALITY (STRUCTURE CODE '+' OR '-').
   PHI PSI - IUPAC PEPTIDE BACKBONE TORSION ANGLES
   X-CA Y-CA Z-CA -  ECHO OF C-ALPHA ATOM COORDINATES              */
/*--------------------------------------------------------------------*/
/* WORDS OF CAUTION */
/* THE VALUES FOR SOLVENT EXPOSURE MAY NOT MEAN WHAT YOU THINK!
    (A) EFFECTS LEADING TO LARGER THAN EXPECTED VALUES:
     SOLVENT EXPOSURE CALCULATION IGNORES UNUSUAL RESIDUES, LIKE ACE,
     OR RESIDUES WITH INCOMPLETE BACKBONE, LIKE ALA 1 OF DATA SET 1CPA.
     IT ALSO IGNORES HETATOMS, LIKE A HEME OR METAL LIGANDS.
     ALSO, SIDE CHAINS MAY BE INCOMPLETE (AN ERROR MESSAGE IS WRITTEN).
    (B) EFFECTS LEADING TO SMALLER THAN EXPECTED VALUES:
     IF YOU APPLY THIS PROGRAM TO PROTEIN DATA BANK DATA SETS
     CONTAINING OLIGOMERS, SOLVENT EXPOSURE IS FOR THE ENTIRE ASSEMBLY,
     NOT FOR THE MONOMER. ALSO, ATOM OXT OF C-TERMINAL RESIDUES IS
     TREATED LIKE A SIDE CHAIN ATOM IF IT IS LISTED AS PART OF THE LAST
     RESIDUE. ALSO, PEPTIDE SUBSTRATES, WHEN LISTED AS ATOMS RATHER THAN
     HETATOMS, ARE TREATED AS PART OF THE PROTEIN, E.G. RESIDUES 499 S
     AND 500 S IN 1CPA.                                               */
/* UNKNOWN OR UNUSUAL RESIDUES ARE NAMED X ON OUTPUT AND THEY ARE
   NOT CHECKED FOR STANDARD NUMBER OF SIDECHAIN ATOMS.                */
/* ALL EXPLICIT WATER MOLECULES, LIKE OTHER HETATOMS, ARE IGNORED.    */
/* END OF INTRODUCTORY COMMENTS */
/**********************************************************************/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
/* #include "License.h" * JL */
/* #include "lock.h" * JL */
#include "p2c.h"
#include "Date.h"
#include "Dssp.h"
#include "AccSurf.h"
#include "Vector.h"

/*--------------------------------------------------------------------*/
/* PROGRAM FATAL ERROR EXIT LABEL */
/*******************  MATHEMATICAL CONSTANTS  **************************
  --------------------------------------------------------------------*/
#define DSSP_VERSION_OLD  "**** SECONDARY STRUCTURE DEFINITION \
BY THE PROGRAM DSSP, VERSION JULY 1995 ****"
#define DSSP_VERSION      "==== Secondary Structure Definition \
by the program DSSP, Version July 1995 ===="
#define PIHALF          1.570796
#define PI              3.141593
#define TWOPI           6.283185
#define RADIAN          57.29578
/***/
/***************  ARRAY DIMENSIONING CONSTANTS  ***********************
 MAXBRIDGE- MAXIMUM NUMBER OF BRIDGES IN ARRAY BRIDGETABLE
 NFACE,   - NUMBER OF FACES OF POLYHEDRON. THE COORDINATES OF THE CENTRE
 ORDER      OF EACH TRIANGULAR FACE ARE STORED IN ARRAY P, THE AREA
             IS STORED IN ARRAY WP IN PROCEDURE FLAGACCESS. NFACE MUST
             BE OF THE FORM NFACE=20*(4**ORDER), ORDER=0,1,2,...
             THE ACCURACY OF THE SOLVENT ACCESSIBLE SURFACE OF EACH
             AMINOACID RESIDUE IS ONE ANGSTROM**2 FOR ORDER=2,NFACE=320.
 MAXHIST  - NUMBER OF SLOTS IN ARRAYS HELIXHIST AND BETAHIST USED FOR
             LENGTH STATISTICS OF SECONDARY STRUCTURE.
 MAXSS    - MAXIMUM NUMBER OF SSBOND RECORDS ON INPUT FILE. THE
             DISULFIDE BOND ARE SAVED IN ARRAY SSBONDS.
  --------------------------------------------------------------------*/
#define MAXBRIDGE       300
#define ORDER           2
#define MAXHIST         30
#define MAXSS           100
/***/
/*********************  PHYSICAL CONSTANTS   **************************
 RN       - RADIUS OF PEPTIDE NITROGEN ATOM
 RCA      - RADIUS OF PEPTIDE ALPHA-CARBON ATOM
 RC       - RADIUS OF PEPTIDE C'-CARBON ATOM
 RO       - RADIUS OF PEPTIDE OXYGEN ATOM
 RSIDEATOM- RADIUS OF SIDECHAIN ATOM
 RWATER   - RADIUS OF WATER MOLECULE
 SSDIST   - MAXIMUM ALLOWED DISTANCE OF DISULFIDE BRIDGE
 BREAKDIST- MAXIMUM ALLOWED PEPTIDE BOND LENGTH. IF DISTANCE IS
             GREATER A POLYPEPTIDE CHAIN INTERRUPTION IS ASSUMED.
 RESRAD   - MAXIMUM RADIUS OF A SPHERE AROUND C-ALPHA CONTAINING
             ALL ATOMS OF A RESIDUE
 CADIST   - MINIMUM DISTANCE BETWEEN ALPHA-CARBON ATOMS SUCH THAT NO
             BACKBONE HYDROGEN BONDS CAN BE FORMED
 DIST     - SMALLEST ALLOWED DISTANCE BETWEEN ANY ATOMS
 MAXDIST  - LARGEST ALLOWED DISTANCE BETWEEN SIDECHAIN ATOM AND C-ALPHA
             WITHIN A RESIDUE
 Q        - COUPLING CONSTANT FOR ELECTROSTATIC ENERGY
                    Q=-332*0.42*0.2*1000.0
 HBLOW    - LOWEST ALLOWED  ENERGY OF A HYDROGEN BOND IN CAL/MOL
 HBHIGH   - HIGHEST ALLOWED ENERGY OF A HYDROGEN BOND IN CAL/MOL
  --------------------------------------------------------------------*/

#define RN              1.65
#define RCA             1.87
#define RC              1.76
#define RO              1.4
#define RSIDEATOM       1.8
#define RWATER          1.4
#define SSDIST          3.0
#define BREAKDIST       2.5
#define RESRAD          10.0
#define CADIST          9.0
#define DIST            0.5
#define MAXDIST         10.0
#define Q               (-27888.0)

#define HBLOW           (-9900)
#define HBHIGH          (-500)


/***/
/***************** GLOBAL DATA TYPE DEFINITIONS ***********************/


typedef enum {
  parallel, antiparallel, nobridge
} bridgetyp;

typedef long bridgeset[MAXBRIDGE / 32 + 2];


typedef struct Bridge {
  char sheetname, laddername;
  bridgetyp btyp;
  bridgeset linkset;
  long ib, ie, jb, je, from, towards;
} Bridge;


Backbone chain[NMAX + 1];
long lchain;
Vector sidechain[MAXATOM];

static long nss, nssintra, nssinter, nbridge;
static Char6 ssbonds[MAXSS][2];

static FILE *tapein, *tapeout;
static Bridge bridgetable[MAXBRIDGE];

static int gClassicFormatFlag;

static char* gDSSP_Version=DSSP_VERSION;

/***/
/*-------------------------------------------------------------------*/
/* PROCEDURE DATE(VAR DATESTRING:PACKED ARRAY[1..11] OF CHAR);EXTERN; */
/* activate DATE by removing comment brackets if necessary */
/***/

static double Atan2(double y, double x)
{
  double z;

  if (x != 0.0)
    z = atan(y / x);
  else if (y > 0.0)
    z = PIHALF;
  else if (y < 0.0)
    z = -PIHALF;
  else
    z = TWOPI;
  if (x >= 0.0)
    return z;
  if (y > 0.0)
    z += PI;
  else
    z -= PI;
  return z;
}  /* Atan2 */


/***/

static double Dihedralangle(double *v1, double *v2, double *v3, double *v4)
{
  /*CALCULATES TORSION ANGLE OF A SET OF 4 ATOMS V1-V2-V3-V4.
    DIHEDRALANGLE IS THE ANGLE BETWEEN THE PROJECTION OF
    V1-V2 AND THE PROJECTION OF V4-V3 ONTO A PLANE NORMAL TO
    BOND V2-V3.*/
  /***/
  double Result, u, v;
  Vector v12, v43, x, y, z, p;

  Diff(v1, v2, v12);
  Diff(v4, v3, v43);
  Diff(v2, v3, z);
  Cross(z, v12, p);
  Cross(z, v43, x);
  Cross(z, x, y);
  u = Dot(x, x);
  v = Dot(y, y);
  Result = 360.0;
  if (u <= 0.0 || v <= 0.0)
    return Result;
  u = Dot(p, x) / sqrt(u);
  v = Dot(p, y) / sqrt(v);
  if (u != 0.0 || v != 0.0)
    return (Atan2(v, u) * RADIAN);
  return Result;
}  /* Dihedralangle alob*/


/***/

static double Cosangle(double *v1, double *v2, double *v3, double *v4)
{
  Vector u, v;
  double x;

  Diff(v1, v2, u);
  Diff(v3, v4, v);
  x = Dot(u, u) * Dot(v, v);
  if (x > 0.0)
    return (Dot(u, v) / sqrt(x));
  else
    return 0.0;
}  /* Cosangle */


/***/



/*--------------------------------------------------------------------*/

static int Nochainbreak(long i, long j)
{
  long k;
  int test;

  test = (i >= 1 && j <= NMAX && i <= j);
  k = i;
  while (test && k <= j) {
    if (chain[k].aa == '!')
      test = false;
    else
      k++;
  }
  return test;
}  /* Nochainbreak */


/***/
/*--------------------------------------------------------------------*/

static void Writeresidue(Backbone res)
{
  long i;

  for (i = 0; i <= 3; i++)
    fprintf(stderr,"%c",res.threelettercode[i]);
  for (i = 0; i <= 5; i++)
    fprintf(stderr,"%c",res.aaident[i]);
}  /* Writeresidue */



#define MAXSIDEATOMS    20


typedef enum {
  headercard, compndcard, sourcecard, authorcard, ssbondcard, atomcard,
  tercard, endcard, othercard
} cardtype;
/***/

typedef struct cardcontents {
  cardtype art;
  union {
    char z[128];
    Char6 r[2];
    struct {
      Char4 atomname, aaname;
      char altloc, residuename;
      Char6 reseqnum;
      Vector coordinates;
    } U5;
    char ch;
  } UU;
} cardcontents;   /* CARDCONTENTS TYPE DEFINITION */

/***/


/* static variables for Inputcoordinates: */
struct LOC_Inputcoordinates {
  long *lchain, latom, hatoms;
  int nmissing, camissing, cmissing, omissing, modelfound;
      /*MS Flag to indicate that a MODEL line has been found. */
  int corelimit;
  Vector sidecoordinates[MAXSIDEATOMS];
  double dco;
  Char4 sideatomnames[MAXSIDEATOMS];
  Backbone reszero, resinfo;
} ;

/***/

static char Onelettercode(char *aaa, struct LOC_Inputcoordinates *LINK)
{
  char aasymbol[50];
  char aminoacid[150];
  char string[5][30];
  long i, l, k;
  char a;

  memcpy(aasymbol, "ARNDCEQGHILKMFPSTWYVBZXXXXXXXXXXXXXXXX--CCCCIPPPW-", 50);
  memcpy(string[0], "ALAARGASNASPCYSGLUGLNGLYHISILE", 30);
  memcpy(string[1], "LEULYSMETPHEPROSERTHRTRPTYRVAL", 30);
  memcpy(string[2], "ASXGLXACDALBALIABUAROBASBETHSE", 30);
  memcpy(string[3], "HYPHYLORNPCASARTAUTHYUNKACEFOR", 30);
  memcpy(string[4], "CYHCSHCSSCYXILUPRZPR0CPRTRYHOH", 30);
  l = 0;
  for (k = 0; k <= 4; k++) {
    for (i = 0; i <= 29; i++) {
      l++;
      aminoacid[l - 1] = string[k][i];
    }
  }
  a = '-';
  i = 1;
  k = 1;
  while (k < 51 && a == '-') {
    if (aminoacid[i - 1] == aaa[0]) {
      if (aminoacid[i] == aaa[1]) {
	if (aminoacid[i + 1] == aaa[2])
	  a = aasymbol[k - 1];
      }
    }
    i += 3;
    k++;
  }
  if (a == '-')   /*MS Let's assume that anything is an aminoacid. */
    a = 'X';
  return a;
}  /* Onelettercode */

/* static variables for Checksideatoms: */
struct LOC_Checksideatoms {
  struct LOC_Inputcoordinates *LINK;
} ;

/***/

static void Checkdist(Backbone *resinfo, struct LOC_Checksideatoms *LINK)
{
  long i, j, FORLIM;

  i = 1;
  while (i <= resinfo->nsideatoms) {
    if (Distance(resinfo->ca, LINK->LINK->sidecoordinates[i - 1]) <= MAXDIST) {
      i++;
      continue;
    }
    fprintf(stderr," !!! Residue ");
    Writeresidue(*resinfo);
    fprintf(stderr," has illegal sidechain atom named ");
    for (j = 0; j <= 3; j++)
      fprintf(stderr,"%c",LINK->LINK->sideatomnames[i - 1][j]);
    fprintf(stderr,".\n");
    fprintf(stderr,"     This atom will be ignored !!!\n\n");
    FORLIM = resinfo->nsideatoms;
    for (j = i + 1; j <= FORLIM; j++) {
      memcpy(LINK->LINK->sideatomnames[j - 2],
	     LINK->LINK->sideatomnames[j - 1], sizeof(Char4));
      memcpy(LINK->LINK->sidecoordinates[j - 2],
	     LINK->LINK->sidecoordinates[j - 1], sizeof(Vector));
    }
    resinfo->nsideatoms--;
  }
}  /* Checkdist */

/***/

static void Checksideatoms(Backbone *resinfo, struct LOC_Inputcoordinates *LINK)
{
  struct LOC_Checksideatoms V;
  long i, j;
  char c;

  /***/

  V.LINK = LINK;
  Checkdist(resinfo, &V);
  i = -1;
  c = resinfo->aa;
  if (c == 'G')
    i = 0;
  if (c == 'A')
    i = 1;
  if (c == 'S' || c == 'C')
    i = 2;
  if (c == 'P' || c == 'T' || c == 'V')
    i = 3;
  if (c == 'B' || c == 'M' || c == 'L' || c == 'I' || c == 'D' || c == 'N')
    i = 4;
  if (c == 'Z' || c == 'K' || c == 'Q' || c == 'E')
    i = 5;
  if (c == 'H')
    i = 6;
  if (c == 'F' || c == 'R')
    i = 7;
  if (c == 'Y')
    i = 8;
  if (c == 'W')
    i = 10;
  if (resinfo->nsideatoms < i) {
    fprintf(stderr," !!! Residue ");
    Writeresidue(*resinfo);
    fprintf(stderr," has%3ld instead of expected ", resinfo->nsideatoms);
    fprintf(stderr,"%3ld sidechain atoms.\n", i);
    fprintf(stderr,
      "     Calculated solvent accessibility refers to incomplete sidechain !!!\n\n");
  }
  if (i == -1 || resinfo->nsideatoms <= i)
    return;
  fprintf(stderr," !!! Residue ");
  Writeresidue(*resinfo);
  fprintf(stderr," has%3ld instead of expected ", resinfo->nsideatoms);
  fprintf(stderr,"%3ld sidechain atoms.\n", i);
  fprintf(stderr,"     last sidechain atom name is ");
  for (j = 0; j <= 3; j++)
    fprintf(stderr,"%c",LINK->sideatomnames[resinfo->nsideatoms - 1][j]);
  fprintf(stderr,"\n     calculated solvent accessibility includes extra atoms !!!\n\n");
}  /* Checksideatoms */

/***/

static void Putresidue(struct LOC_Inputcoordinates *LINK)
{
  /* insert residue into protein chain */
  long i;
  int complete;
  long FORLIM;

  complete = !(LINK->nmissing || LINK->camissing || LINK->cmissing ||
	       LINK->omissing);
  if (!complete &&
      strncmp(LINK->reszero.aaident, LINK->resinfo.aaident, sizeof(Char6)) &&
      LINK->resinfo.aa != 'X') {
    fprintf(stderr," !!! Backbone incomplete for residue ");
    Writeresidue(LINK->resinfo);
    fprintf(stderr,"\n     residue will be ignored !!!\n\n");
  }
  /*MS if X has no Backbone it is no aminoacid!*/
  LINK->corelimit = (LINK->latom + LINK->resinfo.nsideatoms > MAXATOM ||
		     *LINK->lchain > NMAX - 2);
  if (complete && !LINK->corelimit) {
    Checksideatoms(&LINK->resinfo, LINK);
    memcpy(LINK->resinfo.h, LINK->resinfo.n, sizeof(Vector));
    if (Nochainbreak(*LINK->lchain, *LINK->lchain)) {
      if (Distance(chain[*LINK->lchain].c, LINK->resinfo.n) > BREAKDIST)
	  /* keep ! at LCHAIN */
	  {  /* CS Oct 1987 */
	fprintf(stderr," !!! Excessive C to N distance ");
	fprintf(stderr,"% .5E>% .5E\n",
	       Distance(chain[*LINK->lchain].c, LINK->resinfo.n), BREAKDIST);
	fprintf(stderr,"     before residue ");
	Writeresidue(LINK->resinfo);
	fprintf(stderr,". chain break residue inserted !!!\n\n");
	(*LINK->lchain)++;
      }
    }
    if (Nochainbreak(*LINK->lchain, *LINK->lchain) && LINK->resinfo.aa != 'P') {
      LINK->dco = Distance(chain[*LINK->lchain].c, chain[*LINK->lchain].o);
      for (i = 0; i <= 2; i++)
	LINK->resinfo.h[i] = LINK->resinfo.n[i] +
	    (chain[*LINK->lchain].c[i] - chain[*LINK->lchain].o[i]) / LINK->dco;
    }
    (*LINK->lchain)++;
    chain[*LINK->lchain] = LINK->resinfo;
    FORLIM = LINK->resinfo.nsideatoms;
    for (i = 0; i < FORLIM; i++) {
      memcpy(sidechain[LINK->latom + i], LINK->sidecoordinates[i],
	     sizeof(Vector));
      sidechain[LINK->latom + i][3] = RSIDEATOM;
    }
    LINK->latom += LINK->resinfo.nsideatoms;
  }
  if (Nochainbreak(*LINK->lchain, *LINK->lchain) && !complete)
    (*LINK->lchain)++;
  LINK->resinfo = LINK->reszero;
  LINK->nmissing = true;
  LINK->camissing = true;
  LINK->cmissing = true;
  LINK->omissing = true;
}  /* Putresidue */

/***/

static void Getresidue(char *atomname, double *coordinates,
		      struct LOC_Inputcoordinates *LINK)
{
  int hydrogenatom;

  hydrogenatom = ((atomname[0] == '9' || atomname[0] == '8' ||
		   atomname[0] == '7' || atomname[0] == '6' ||
		   atomname[0] == '5' || atomname[0] == '4' ||
		   atomname[0] == '3' || atomname[0] == '2' ||
		   atomname[0] == '1' || atomname[0] == '0' ||
		   atomname[0] == ' ') &&
		  (atomname[1] == 'D' || atomname[1] == 'H'));
  if (hydrogenatom) {
    LINK->hatoms++;
    return;
  }
  if (!strncmp(atomname, " N  ", sizeof(Char4))) {
    LINK->nmissing = false;
    memcpy(LINK->resinfo.n, coordinates, sizeof(Vector));
    LINK->resinfo.n[3] = RN;
    return;
  }
  if (!strncmp(atomname, " CA ", sizeof(Char4))) {
    LINK->camissing = false;
    memcpy(LINK->resinfo.ca, coordinates, sizeof(Vector));
    LINK->resinfo.ca[3] = RCA;
    return;
  }
  if (!strncmp(atomname, " C  ", sizeof(Char4))) {
    LINK->cmissing = false;
    memcpy(LINK->resinfo.c, coordinates, sizeof(Vector));
    LINK->resinfo.c[3] = RC;
    return;
  }
  if (!strncmp(atomname, " O  ", sizeof(Char4))) {
    LINK->omissing = false;
    memcpy(LINK->resinfo.o, coordinates, sizeof(Vector));
    LINK->resinfo.o[3] = RO;
    return;
  }
  if (LINK->resinfo.nsideatoms >= MAXSIDEATOMS)
    return;
  LINK->resinfo.nsideatoms++;
  memcpy(LINK->sidecoordinates[LINK->resinfo.nsideatoms - 1], coordinates,
	 sizeof(Vector));
  memcpy(LINK->sideatomnames[LINK->resinfo.nsideatoms - 1], atomname,
	 sizeof(Char4));
}  /* Getresidue */

/***/

static void Readcard(cardcontents *cardinfo, struct LOC_Inputcoordinates *LINK)
{
  char c;
  long k, l, m;
  Char6 key;

  cardinfo->art = othercard;
  do {
    if (!P_eof(tapein)) {
      *key = getc(tapein);
      if (key[0] == '\n')
	key[0] = ' ';
    }
  } while (!(isupper(key[0]) | P_eof(tapein)));
  if (P_eof(tapein)) {
    cardinfo->art = endcard;
    return;
  }
  for (l = 1; l <= 5; l++) {
    if (!P_eoln(tapein)) {
      key[l] = getc(tapein);
      if (key[l] == '\n')
	key[l] = ' ';
    }
  }
  if (!strncmp(key, "HEADER", sizeof(Char6)))
    cardinfo->art = headercard;
  if (!strncmp(key, "COMPND", sizeof(Char6)))
    cardinfo->art = compndcard;
  if (!strncmp(key, "SOURCE", sizeof(Char6)))
    cardinfo->art = sourcecard;
  if (!strncmp(key, "AUTHOR", sizeof(Char6)))
    cardinfo->art = authorcard;
  if (!strncmp(key, "SSBOND", sizeof(Char6)))
    cardinfo->art = ssbondcard;
  if (!strncmp(key, "ATOM  ", sizeof(Char6)))
    cardinfo->art = atomcard;
  if (!strncmp(key, "HETATM", sizeof(Char6)))
	/*MS Let's look also at HETATM's (e.g.: 1amt) */
	  cardinfo->art = atomcard;
  if (!strncmp(key, "TER   ", sizeof(Char6)))
    cardinfo->art = tercard;
  if (!strncmp(key, "END   ", sizeof(Char6)))
    cardinfo->art = endcard;
  if (!strncmp(key, "MODEL ", sizeof(Char6))) {  /*MS deal with NMR models */
    if (LINK->modelfound)  /* don't read more than 1 model */
      cardinfo->art = endcard;
    else
      LINK->modelfound = true;
  }
  if (!strncmp(key, "ENDMDL", sizeof(Char6)))
	/*MS Terminate: read only one MODEL! */
	  cardinfo->art = endcard;
  switch (cardinfo->art) {

  case headercard:
  case compndcard:
  case sourcecard:
  case authorcard:
    for (l = 0; l <= 5; l++)
      cardinfo->UU.z[l] = key[l];
    for (l = 6; l <= 126; l++)
      cardinfo->UU.z[l] = ' ';
    cardinfo->UU.z[127] = '.';
    if (cardinfo->art == headercard)
      m = 66;
    else
      m = 70;
    for (l = 6; l < m; l++) {
      if (!P_eoln(tapein)) {
	cardinfo->UU.z[l] = getc(tapein);
	if (cardinfo->UU.z[l] == '\n')
	  cardinfo->UU.z[l] = ' ';
      }
    }
    break;

  case ssbondcard:
    for (l = 7; l <= 8; l++) {
      c = getc(tapein);
      if (c == '\n')
	c = ' ';
    }
    for (k = 0; k <= 1; k++) {
      for (l = 1; l <= 7; l++) {
	c = getc(tapein);
	if (c == '\n')
	  c = ' ';
      }
      cardinfo->UU.r[k][5] = getc(tapein);
      c = getc(tapein);
      if (cardinfo->UU.r[k][5] == '\n')
	cardinfo->UU.r[k][5] = ' ';
      if (c == '\n')
	c = ' ';
      /* minor modification suggested by Steven Sheriff */
      for (l = 0; l <= 3; l++) {
	cardinfo->UU.r[k][l] = getc(tapein);
	if (cardinfo->UU.r[k][l] == '\n')
	  cardinfo->UU.r[k][l] = ' ';
      }
      if (P_eoln(tapein))
	cardinfo->UU.r[k][4] = ' ';
      else {
	cardinfo->UU.r[k][4] = getc(tapein);
	if (cardinfo->UU.r[k][4] == '\n')
	  cardinfo->UU.r[k][4] = ' ';
      }
    }
    /* end minor modification suggested by Steven Sheriff */
    break;

  case atomcard:
    for (l = 7; l <= 12; l++) {
      c = getc(tapein);
      if (c == '\n')
	c = ' ';
    }
    for (l = 0; l <= 3; l++) {
      cardinfo->UU.U5.atomname[l] = getc(tapein);
      if (cardinfo->UU.U5.atomname[l] == '\n')
	cardinfo->UU.U5.atomname[l] = ' ';
    }
    cardinfo->UU.U5.altloc = getc(tapein);
    if (cardinfo->UU.U5.altloc == '\n')
      cardinfo->UU.U5.altloc = ' ';
    for (l = 0; l <= 2; l++) {
      cardinfo->UU.U5.aaname[l] = getc(tapein);
      if (cardinfo->UU.U5.aaname[l] == '\n')
	cardinfo->UU.U5.aaname[l] = ' ';
    }
    cardinfo->UU.U5.aaname[3] = ' ';
    cardinfo->UU.U5.residuename = Onelettercode(cardinfo->UU.U5.aaname, LINK);
    c = getc(tapein);
    cardinfo->UU.U5.reseqnum[5] = getc(tapein);
    if (c == '\n')
      c = ' ';
    if (cardinfo->UU.U5.reseqnum[5] == '\n')
      cardinfo->UU.U5.reseqnum[5] = ' ';
    for (l = 0; l <= 4; l++) {
      cardinfo->UU.U5.reseqnum[l] = getc(tapein);
      if (cardinfo->UU.U5.reseqnum[l] == '\n')
	cardinfo->UU.U5.reseqnum[l] = ' ';
    }
    for (l = 0; l <= 2; l++)
      fscanf(tapein, "%lf", &cardinfo->UU.U5.coordinates[l]);
    break;

  case tercard:
  case endcard:
  case othercard:
    /* blank case */
    break;
  }
  fscanf(tapein, "%*[^\n]");
  getc(tapein);
}  /* Readcard */


/*--------------------------------------------------------------------*/
/* SEE BROOKHAVEN PROTEIN DATA BANK ATOMIC COORDINATE ENTRY FORMAT
                                    OF DEC. 1981.
   -------------------------------------------------------------------*/

static void Inputcoordinates(long *lchain_)
{
  struct LOC_Inputcoordinates V;
  char datestring[11];
  long i, j;
  int finish;
  structure s;
  cardtype ctype;
  cardcontents cardinfo;
  long cardhist[(long)othercard - (long)headercard + 1];

  /***/

  V.lchain = lchain_;
  nss = 0;
  V.latom = 0;
  V.hatoms = 0;
  V.modelfound = false;   /*MS init */
  for (j = 0; j <= 5; j++)
    V.reszero.aaident[j] = ' ';
  V.reszero.aa = '!';
  V.reszero.access = 0;
  memmove(V.reszero.threelettercode, "    ", sizeof(Char4));

/*
  struct LOC_Inputcoordinates {
  long *lchain, latom, hatoms;
  int nmissing, camissing, cmissing, omissing, modelfound;
  int corelimit;
  Vector sidecoordinates[MAXSIDEATOMS];
  double dco;
  Char4 sideatomnames[MAXSIDEATOMS];
  Backbone reszero, resinfo;
} ;
  */
  for (s = symbol; (long)s <= (long)beta2; s = (structure)((long)s + 1))
    V.reszero.ss[(long)s - (long)symbol] = ' ';
  V.reszero.sheetlabel = ' ';
  V.reszero.partner[0] = 0;
  V.reszero.partner[(long)beta2 - (long)beta1] = 0;
  V.reszero.alpha = 360.0;
  V.reszero.kappa = 360.0;
  for (j = 0; j <= 1; j++) {
    V.reszero.acceptor[j].residue = 0;
    V.reszero.acceptor[j].energy = 0;
    V.reszero.donor[j].residue = 0;
    V.reszero.donor[j].energy = 0;
  }
  V.reszero.atompointer = 0;
  V.reszero.nsideatoms = 0;
  for (j = 0; j <= 2; j++) {
    V.reszero.h[j] = 0.0;
    V.reszero.h[3] = 0.0;
    V.reszero.n[j] = 0.0;
    V.reszero.n[3] = 0.0;
    V.reszero.ca[j] = 0.0;
    V.reszero.ca[3] = 0.0;
    V.reszero.c[j] = 0.0;
    V.reszero.c[3] = 0.0;
    V.reszero.o[j] = 0.0;
    V.reszero.o[3] = 0.0;
  }
  for (i = 0; i <= NMAX; i++)
    chain[i] = V.reszero;
  Date(datestring);    /* DATE(DAY-MONTH-YEAR); */
  /* comment out this line if necessary */
  fprintf(tapeout, gDSSP_Version);
  fprintf(tapeout, " DATE=%.11s", datestring);
  for (i = 106; i <= 127; i++)
    putc(' ', tapeout);
  fprintf(tapeout, ".\n");
  fprintf(tapeout, "REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS ");
  fprintf(tapeout, "22 (1983) 2577-2637");
  for (i = 66; i <= 127; i++)
    putc(' ', tapeout);
  fprintf(tapeout, ".\n");
  for (ctype = headercard;
       (long)ctype <= (long)othercard;
       ctype = (cardtype)((long)ctype + 1))
    cardhist[(long)ctype - (long)headercard] = 0;
  V.corelimit = false;
  finish = false;
  V.resinfo = V.reszero;
  V.nmissing = true;
  V.camissing = true;
  V.cmissing = true;
  V.omissing = true;

  do {
    Readcard(&cardinfo, &V);
    cardhist[(long)cardinfo.art - (long)headercard]++;
    switch (cardinfo.art) {

    case headercard:
    case compndcard:
    case sourcecard:
    case authorcard:
      if (cardhist[(long)cardinfo.art - (long)headercard] == 1) {
	for (i = 0; i <= 127; i++)
	  putc(cardinfo.UU.z[i], tapeout);
	putc('\n', tapeout);
      }
      break;

    case ssbondcard:
      nss++;
      for (i = 0; i <= 1; i++)
	memcpy(ssbonds[nss - 1][i], cardinfo.UU.r[i], sizeof(Char6));
      break;

    case atomcard:
      /*      IF (NOT (nmissing OR camissing OR cmissing OR omissing)
                AND (altloc IN [' ', 'A'])) THEN */
      if (cardinfo.UU.U5.residuename != '-' &&
	  (cardinfo.UU.U5.altloc == 'A' || cardinfo.UU.U5.altloc == ' ')) {
	if (strncmp(V.resinfo.aaident, cardinfo.UU.U5.reseqnum, sizeof(Char6))) {
	  Putresidue(&V);
	  V.resinfo.atompointer = V.latom;
	  memcpy(V.resinfo.aaident, cardinfo.UU.U5.reseqnum, sizeof(Char6));
	  V.resinfo.aa = cardinfo.UU.U5.residuename;
	  memcpy(V.resinfo.threelettercode, cardinfo.UU.U5.aaname,
		 sizeof(Char4));
	}
	Getresidue(cardinfo.UU.U5.atomname, cardinfo.UU.U5.coordinates, &V);
      }
      if (cardinfo.UU.U5.residuename == '-') {
	fprintf(stderr," !!! Residue ");
	for (i = 0; i <= 3; i++)
	  fprintf(stderr,"%c",cardinfo.UU.U5.aaname[i]);
	for (i = 0; i <= 5; i++)
	  fprintf(stderr,"%c",cardinfo.UU.U5.reseqnum[i]);
	fprintf(stderr," has nonstandard name.\n");
	fprintf(stderr,"     residue will be ");
	fprintf(stderr,"ignored !!!\n");
      }
      if (cardinfo.UU.U5.altloc != 'A' && cardinfo.UU.U5.altloc != ' ') {
	/* ????????? */
	fprintf(stderr," !!! In residue");
	for (i = 0; i <= 3; i++)
	  fprintf(stderr," %c", cardinfo.UU.U5.aaname[i]);
	for (i = 0; i <= 5; i++)
	  fprintf(stderr,"%c",cardinfo.UU.U5.reseqnum[i]);
	fprintf(stderr," alternate location indicator ");
	fprintf(stderr,"is %c and\n", cardinfo.UU.U5.altloc);
	fprintf(stderr,"     not blank or A. Atom ");
	fprintf(stderr,"named ");
	for (i = 0; i <= 3; i++)
	  fprintf(stderr,"%c",cardinfo.UU.U5.atomname[i]);
	fprintf(stderr," will be ignored !!!\n\n");
      }
      break;

    case tercard:
      Putresidue(&V);
      break;

    case endcard:
      finish = true;
      Putresidue(&V);
      break;

    case othercard:
      /* blank case */
      break;
    }
  } while (!(V.corelimit || finish));
  if (V.corelimit) {
    fprintf(stderr," !!! Number of atoms or residues exceeds ");
    fprintf(stderr,"storage capacity !!!\n");
  }
  if (!Nochainbreak(*V.lchain, *V.lchain))
    (*V.lchain)--;
  if (V.hatoms > 0) {
    fprintf(stderr," !!! %12ld hydrogen or deuterium atoms were ignored\n", V.hatoms);
    fprintf(stderr,"     in the calculation of side chain solvent \n");
    fprintf(stderr,"     accessibility !!!\n");
  }
  if (cardhist[0] < 1)
    fprintf(stderr," !!! HEADER-card missing !!!\n");
  if (cardhist[(long)compndcard - (long)headercard] < 1)
    fprintf(stderr," !!! COMPOUND-card missing !!!\n");
  if (cardhist[(long)sourcecard - (long)headercard] < 1)
    fprintf(stderr," !!! SOURCE-card missing !!!\n");
  if (cardhist[(long)authorcard - (long)headercard] < 1)
    fprintf(stderr," !!! AUTHOR-card missing !!!\n");
  if (*V.lchain < 1) {
    fprintf(stderr," !!! No residue with complete backbone !!!\n");
    exit(1);
  }
  if (V.latom == 0)
    fprintf(stderr," !!! All sidechain coordinates missing !!!\n");
}  /* Inputcoordinates */

#undef MAXSIDEATOMS


/***/
/*--------------------------------------------------------------------*/

static int Testbond(long i, long j)
{
  /* TESTBOND IS TRUE IF I IS DONOR[=NH] TO J, OTHERWISE FALSE */
  Backbone *WITH;

  WITH = &chain[i];
  return (WITH->acceptor[0].residue == j && WITH->acceptor[0].energy < HBHIGH ||
	  WITH->acceptor[1].residue == j && WITH->acceptor[1].energy < HBHIGH);
}  /* Testbond */


/***/

static int Testssbond(long i, long j)
{
  int ssbond;
  long k;

  ssbond = false;
  k = 1;
  if (!(Nochainbreak(i, i) & Nochainbreak(j, j)))
    return ssbond;
  while (!(ssbond || k > nss)) {
    ssbond = ((!strncmp(chain[i].aaident, ssbonds[k - 1][0], sizeof(Char6)) &&
	       !strncmp(chain[j].aaident, ssbonds[k - 1][1], sizeof(Char6))) ||
	      (!strncmp(chain[i].aaident, ssbonds[k - 1][1], sizeof(Char6)) &&
	       !strncmp(chain[j].aaident, ssbonds[k - 1][0], sizeof(Char6))));
    k++;
  }
  return ssbond;
}  /* Testssbond */


/***/
/*--------------------------------------------------------------------*/

static void Flagssbonds()
{
  int ssbond;
  char cc;
  long i, j, ii, jj;
  double d;
  long FORLIM;
  Backbone *WITH;
  long FORLIM1;

  /***/

  nssintra = 0;
  nssinter = 0;
  cc = 'a' - 1;
  FORLIM = lchain - 2;
  for (i = 1; i <= FORLIM; i++) {
    if (chain[i].aa == 'C' && chain[i].nsideatoms > 1) {
      ii = chain[i].atompointer + 2;
      j = i + 1;
      do {
	j++;
	ssbond = false;
	if (chain[j].nsideatoms > 1 && chain[j].aa == 'C')
	  jj = chain[j].atompointer + 2;
	else
	  jj = 0;
	if (jj > 0)
	  ssbond = (Distance(sidechain[ii - 1], sidechain[jj - 1]) < SSDIST);
      } while (!(ssbond || j == lchain));
      if (ssbond & (!Testssbond(i, j))) {
	fprintf(stderr," !!! Additional ssbond found between ");
	fprintf(stderr,"residues ");
	Writeresidue(chain[i]);
	fprintf(stderr," and ");
	Writeresidue(chain[j]);
	fprintf(stderr," !!!\n\n");
      }
    }
  }
  if (nss > 0) {
    FORLIM = lchain - 2;
    for (i = 1; i <= FORLIM; i++) {
      WITH = &chain[i];
      if (WITH->aa == 'C') {
	FORLIM1 = lchain;
	for (j = i + 2; j <= FORLIM1; j++) {
	  if (chain[j].aa == 'C') {
	    if (Testssbond(i, j)) {
	      if (cc == 'z') {
		fprintf(stderr," !!! SS-bridge label restart at a !!!\n");
		cc = 'a' - 1;
	      }
	      cc++;
	      WITH->aa = cc;
	      chain[j].aa = cc;
	      if (Nochainbreak(i, j))
		nssintra++;
	      else
		nssinter++;
	      if (WITH->nsideatoms > 1) {
		if (chain[j].nsideatoms > 1) {
		  jj = chain[j].atompointer + 2;
		  ii = WITH->atompointer + 2;
		  d = Distance(sidechain[ii - 1], sidechain[jj - 1]);
		  if (d > SSDIST) {
		    fprintf(stderr," !!! SS-bond distance is%5.1f between residues", d);
		    Writeresidue(chain[i]);
		    fprintf(stderr," and ");
		    Writeresidue(chain[j]);
		    fprintf(stderr," !!!\n\n");
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  if (nss != nssintra + nssinter)
    fprintf(stderr," !!! Error in SSBOND data records !!!\n");
}  /* Flagssbonds */


/***/
/*--------------------------------------------------------------------*/

static void Flagchirality()
{
  long i;
  double ckap, skap;
  long FORLIM;
  Backbone *WITH;

  FORLIM = lchain - 2;
  for (i = 2; i <= FORLIM; i++) {
    WITH = &chain[i];
    if (Nochainbreak(i - 1, i + 2)) {
      WITH->alpha = Dihedralangle(chain[i - 1].ca, WITH->ca, chain[i + 1].ca,
				  chain[i + 2].ca);
      if (WITH->alpha < 0.0)
	WITH->ss[(long)chirality - (long)symbol] = '-';
      else
	WITH->ss[(long)chirality - (long)symbol] = '+';
    }
  }
  FORLIM = lchain - 2;
  /***/
  for (i = 3; i <= FORLIM; i++) {
    WITH = &chain[i];
    if (Nochainbreak(i - 2, i + 2)) {
      ckap = Cosangle(chain[i].ca, chain[i - 2].ca, chain[i + 2].ca,
		      chain[i].ca);
      skap = sqrt(1.0 - ckap * ckap);
      WITH->kappa = RADIAN * Atan2(skap, ckap);
    }
  }
}  /* Flagchirality */


/***/

static long Bondenergy(long i, long j)
{
  /*RESIDUE I IS DONOR[=NH],J IS ACCEPTOR[=CO] OF THE PROTON IN THE
     HYDROGEN BOND. THE BONDENERGY IS IN CAL/MOL */
  double dho, dhc, dnc, dno;
  long hbe;
  Backbone *WITH;

  hbe = 0;
  WITH = &chain[i];
  if (WITH->aa == 'P')
    return hbe;
  dho = Distance(WITH->h, chain[j].o);
  dhc = Distance(WITH->h, chain[j].c);
  dnc = Distance(WITH->n, chain[j].c);
  dno = Distance(WITH->n, chain[j].o);
  if (dho < DIST || dhc < DIST || dnc < DIST || dno < DIST)
    hbe = HBLOW;
  else
    hbe = (long)floor(Q / dho - Q / dhc + Q / dnc - Q / dno + 0.5);
  if (hbe > HBLOW)
    return hbe;
  fprintf(stderr," !!! Contact between residues ");
  Writeresidue(chain[i]);
  fprintf(stderr," and ");
  Writeresidue(chain[j]);
  fprintf(stderr,"  too close !!!\n");
  hbe = HBLOW;
  return hbe;
}  /* Bondenergy */

/***/

static void Updatebonds(HydrogenBond *b, HydrogenBond hb)
{
  if (hb.energy < b[0].energy) {
    b[1] = b[0];
    b[0] = hb;
  } else if (hb.energy < b[1].energy)
    b[1] = hb;
}  /* Updatebonds */

/***/

static void Setbonds(long i, long j)
{
  /*I IS NH, J IS CO*/
  HydrogenBond hb;

  hb.energy = Bondenergy(i, j);
  hb.residue = j;
  /* CO(J) IS ACCEPTOR OF NH(I) */
  Updatebonds(chain[i].acceptor, hb);
  hb.residue = i;
  Updatebonds(chain[j].donor, hb);
}  /* Setbond */


/***/
/*--------------------------------------------------------------------*/

static void FlagHydrogenBonds()
{
  long i, j, FORLIM;
  Backbone *WITH;
  long FORLIM1;

  /***/

  FORLIM = lchain;
  for (i = 1; i <= FORLIM; i++) {
    if (Nochainbreak(i, i)) {
      WITH = &chain[i];
      FORLIM1 = lchain;
      for (j = i + 1; j <= FORLIM1; j++) {
	if (Nochainbreak(j, j)) {
	  if (Distance(WITH->ca, chain[j].ca) < CADIST) {
	    Setbonds(i, j);
	    if (j != i + 1)
	      Setbonds(j, i);
	  }
	}
      }
    }
  }
}  /* FlagHydrogenBonds */


/***/

static void Ladder(long i, long j, bridgetyp b)
{
  long k;
  int found;
  Bridge *WITH;

  found = false;
  k = 1;
  if (b == nobridge || i >= j)
    return;
  do {
    WITH = &bridgetable[k - 1];
    if (WITH->ib == 0) {
      WITH->ib = i;
      WITH->ie = i;
      WITH->jb = j;
      WITH->je = j;
      WITH->from = 0;
      WITH->towards = 0;
      WITH->btyp = b;
      nbridge++;
      found = true;
    } else {
      found = (WITH->btyp == b && i == WITH->ie + 1) & Nochainbreak(WITH->ie,
		i) & (((j == WITH->je + 1 && b == parallel) &
		       Nochainbreak(WITH->je, j)) | ((j == WITH->jb - 1 &&
			  b == antiparallel) & Nochainbreak(j, WITH->jb)));
      if (found) {
	WITH->ie++;
	if (b == parallel)
	  WITH->je++;
	else
	  WITH->jb--;
      } else {
	k++;
	if (k > MAXBRIDGE) {
	  fprintf(stderr," !!! Bridgetable overflow !!!\n");
	  exit(1);
	}
      }
    }
  } while (!found);   /* Ladder */
}

/***/

static void Testbridge(long i)
{
  long j1, j2, j;
  bridgetyp b;

  /***/

  j1 = 0;
  j2 = 0;
  j = i + 3;
  if (!Nochainbreak(i - 1, i + 1))
    return;
  while (j2 == 0 && j < lchain) {
    if (Nochainbreak(j - 1, j + 1)) {
      if ((Testbond(i + 1, j) & Testbond(j, i - 1)) |
	  (Testbond(j + 1, i) & Testbond(i, j - 1)))
	b = parallel;
      else if ((Testbond(i + 1, j - 1) & Testbond(j + 1, i - 1)) |
	       (Testbond(j, i) & Testbond(i, j)))
	b = antiparallel;
      else
	b = nobridge;
      if (b != nobridge) {
	if (j1 == 0) {
	  j1 = j;
	  Ladder(i, j, b);
	} else if (j != j1) {
	  j2 = j;
	  Ladder(i, j, b);
	}
      }
    }
    j++;
  }
}  /* Testbridge */

/***/

static void Extendladder()
{
  long i, j, ib1, jb1, je1;
  int bulge;
  long FORLIM;
  Bridge *WITH;

  FORLIM = nbridge;
  for (i = 1; i <= FORLIM; i++) {
    WITH = &bridgetable[i - 1];
    j = i + 1;
    while (j <= nbridge && WITH->towards == 0) {
      ib1 = bridgetable[j - 1].ib;
      jb1 = bridgetable[j - 1].jb;
      je1 = bridgetable[j - 1].je;
      bulge = (Nochainbreak(WITH->ie, ib1) && ib1 - WITH->ie < 6 &&
	       bridgetable[j - 1].btyp == WITH->btyp &&
	       bridgetable[j - 1].from == 0);
      if (bulge) {
	switch (WITH->btyp) {

	case parallel:
	  bulge = (jb1 - WITH->je < 6 && ib1 - WITH->ie < 3 ||
		   jb1 - WITH->je < 3) & Nochainbreak(WITH->je, jb1);
	  break;

	case antiparallel:
	  bulge = (WITH->jb - je1 < 6 && ib1 - WITH->ie < 3 ||
		   WITH->jb - je1 < 3) & Nochainbreak(je1, WITH->jb);
	  break;
	}
      }
      if (bulge) {
	WITH->towards = j;
	bridgetable[j - 1].from = i;
      }
      j++;
    }
  }
  FORLIM = nbridge;
  for (i = 1; i <= FORLIM; i++) {
    WITH = &bridgetable[i - 1];
    if (WITH->from == 0) {
      P_expset(WITH->linkset, 0);
      j = i;
      do {
	P_addset(WITH->linkset, j);
	j = bridgetable[j - 1].towards;
      } while (j != 0);
      j = WITH->towards;
      while (j != 0) {
	P_setcpy(bridgetable[j - 1].linkset, WITH->linkset);
	j = bridgetable[j - 1].towards;
      }
    }
  }
}  /* Extendladder */

/* static variables for Sheet: */
struct LOC_Sheet {
  bridgeset ladderset, sheetset;
} ;

/***/

static int Link(long l1, long l2)
{
  /* LINK IS TRUE IF THERE IS A COMMON RESIDUE IN LADDERS L1 AND L2 */
  long ib1, ie1, jb1, je1, ib2, ie2, jb2, je2;

  ib1 = bridgetable[l1 - 1].ib;
  ie1 = bridgetable[l1 - 1].ie;
  jb1 = bridgetable[l1 - 1].jb;
  je1 = bridgetable[l1 - 1].je;
  ib2 = bridgetable[l2 - 1].ib;
  ie2 = bridgetable[l2 - 1].ie;
  jb2 = bridgetable[l2 - 1].jb;
  je2 = bridgetable[l2 - 1].je;
  return (ie1 >= ib2 && ib1 <= ie2 || ie1 >= jb2 && ib1 <= je2 ||
	  je1 >= ib2 && jb1 <= ie2 || je1 >= jb2 && jb1 <= je2);
}  /* Link */

/***/

static void Findsheet(struct LOC_Sheet *LINK)
{
  long l1, l2;
  int finish;
  long FORLIM, FORLIM1;

  /***/

  P_expset(LINK->sheetset, 0);
  l1 = 0;
  if (*LINK->ladderset != 0) {
    do {
      l1++;
    } while (!P_inset(l1, LINK->ladderset));
  }
  if (l1 > 0)
    P_setcpy(LINK->sheetset, bridgetable[l1 - 1].linkset);
  if (l1 <= 0)
    return;
  do {
    finish = true;
    FORLIM = nbridge;
    for (l1 = 1; l1 <= FORLIM; l1++) {
      if (P_inset(l1, LINK->sheetset)) {
	FORLIM1 = nbridge;
	for (l2 = 1; l2 <= FORLIM1; l2++) {
	  if (P_inset(l2, LINK->ladderset)) {
	    if (Link(l1, l2)) {
	      P_setunion(LINK->sheetset, LINK->sheetset,
			 bridgetable[l2 - 1].linkset);
	      P_setdiff(LINK->ladderset, LINK->ladderset,
			bridgetable[l2 - 1].linkset);
	      finish = false;
	    }
	  }
	}
      }
    }
  } while (!finish);   /* Findsheet */
}

/***/

static void Sheet()
{
  struct LOC_Sheet V;
  long asci, i, j;
  char ccs;
  long FORLIM;
  Bridge *WITH;

  /***/

  P_expset(V.ladderset, 0);
  FORLIM = nbridge;
  for (i = 1; i <= FORLIM; i++)
    P_addset(V.ladderset, i);
  ccs = 'A' - 1;
  asci = 64;
  while (*V.ladderset != 0) {
    ccs++;
    if (ccs == 'Z' + 1)  /* continue with lowercase (don't use []\_')*/
      ccs = 'a';
    if (ccs > 'z') {
      fprintf(stderr," !!! Sheet label restart at A !!!\n");
      ccs = 'A';
    }
    Findsheet(&V);
    FORLIM = nbridge;
    for (i = 1; i <= FORLIM; i++) {
      WITH = &bridgetable[i - 1];
      if (P_inset(i, V.sheetset) && WITH->from == 0) {
	if (asci == 90) {
	  fprintf(stderr," !!! Strand label restart at A !!!\n");
	  asci = 64;
	}
	asci++;
	if (WITH->btyp == parallel)
	  WITH->laddername = (char)(asci + 32);
	else
	  WITH->laddername = (char)asci;
	WITH->sheetname = ccs;
	P_setcpy(WITH->linkset, V.sheetset);
	j = WITH->towards;
	while (j != 0) {
	  bridgetable[j - 1].laddername = WITH->laddername;
	  bridgetable[j - 1].sheetname = WITH->sheetname;
	  P_setcpy(bridgetable[j - 1].linkset, V.sheetset);
	  j = bridgetable[j - 1].towards;
	}
      }
    }
  }
}  /* Sheet */

/***/

static void Markstrands()
{
  long i, j, l, ib0, ie0, jb0, je0;
  structure beta, betai, betaj;
  long iset[(long)beta2 - (long)beta1 + 1][9],
       jset[(long)beta2 - (long)beta1 + 1][9];
  char cc;
  long FORLIM, FORLIM1;
  Bridge *WITH;
  Backbone *WITH1;
  long SET1[3];
  long SET2[3];

  FORLIM = nbridge;
  for (i = 1; i <= FORLIM; i++) {
    if (bridgetable[i - 1].from == 0) {
      j = i;
      for (beta = beta1;
	   (long)beta <= (long)beta2;
	   beta = (structure)((long)beta + 1)) {
	P_expset(iset[(long)beta - (long)beta1], 0);
	P_expset(jset[(long)beta - (long)beta1], 0);
      }
      ib0 = lchain;
      ie0 = 0;
      jb0 = lchain;
      je0 = 0;
      do {
	WITH = &bridgetable[j - 1];
	FORLIM1 = WITH->ie;
	for (l = WITH->ib; l <= FORLIM1; l++) {
	  WITH1 = &chain[l];
	  for (beta = beta1;
	       (long)beta <= (long)beta2;
	       beta = (structure)((long)beta + 1))
	    P_addset(iset[(long)beta - (long)beta1],
		     WITH1->ss[(long)beta - (long)symbol]);
	}
	FORLIM1 = WITH->je;
	for (l = WITH->jb; l <= FORLIM1; l++) {
	  WITH1 = &chain[l];
	  for (beta = beta1;
	       (long)beta <= (long)beta2;
	       beta = (structure)((long)beta + 1))
	    P_addset(jset[(long)beta - (long)beta1],
		     WITH1->ss[(long)beta - (long)symbol]);
	}
	if (WITH->ib < ib0)
	  ib0 = WITH->ib;
	if (WITH->ie > ie0)
	  ie0 = WITH->ie;
	if (WITH->jb < jb0)
	  jb0 = WITH->jb;
	if (WITH->je > je0)
	  je0 = WITH->je;
	j = WITH->towards;
      } while (j != 0);
      j = i;
      if (P_setequal(iset[0], P_addset(P_expset(SET1, 0), ' ')))
	betai = beta1;
      else
	betai = beta2;
      if (P_setequal(jset[0], P_addset(P_expset(SET1, 0), ' ')))
	betaj = beta1;
      else
	betaj = beta2;
      if ((!P_setequal(iset[(long)betai - (long)beta1],
		       P_addset(P_expset(SET1, 0), ' '))) |
	  (!P_setequal(jset[(long)betaj - (long)beta1],
		       P_addset(P_expset(SET2, 0), ' '))))
	fprintf(stderr," !!! Strand column overwritten !!!\n");
      do {
	WITH = &bridgetable[j - 1];
	FORLIM1 = WITH->ie;
	for (l = WITH->ib; l <= FORLIM1; l++) {
	  WITH1 = &chain[l];
	  WITH1->ss[(long)betai - (long)symbol] = WITH->laddername;
	  if (WITH->btyp == parallel)
	    WITH1->partner[(long)betai - (long)beta1] = WITH->jb + l - WITH->ib;
	  else
	    WITH1->partner[(long)betai - (long)beta1] = WITH->je - l + WITH->ib;
	}
	FORLIM1 = WITH->je;
	for (l = WITH->jb; l <= FORLIM1; l++) {
	  WITH1 = &chain[l];
	  WITH1->ss[(long)betaj - (long)symbol] = WITH->laddername;
	  if (WITH->btyp == parallel)
	    WITH1->partner[(long)betaj - (long)beta1] = WITH->ib + l - WITH->jb;
	  else
	    WITH1->partner[(long)betaj - (long)beta1] = WITH->ie - l + WITH->jb;
	}
	j = WITH->towards;
      } while (j != 0);
      if (ib0 == ie0)
	cc = 'B';
      else
	cc = 'E';
      for (j = ib0; j <= ie0; j++) {
	WITH1 = &chain[j];
	if (WITH1->ss[0] != 'E')
	  WITH1->ss[0] = cc;
      }
      for (j = jb0; j <= je0; j++) {
	WITH1 = &chain[j];
	if (WITH1->ss[0] != 'E')
	  WITH1->ss[0] = cc;
      }
    }
  }
  FORLIM = nbridge;
  for (j = 0; j < FORLIM; j++) {
    WITH = &bridgetable[j];
    FORLIM1 = WITH->ie;
    for (l = WITH->ib; l <= FORLIM1; l++)
      chain[l].sheetlabel = WITH->sheetname;
    FORLIM1 = WITH->je;
    for (l = WITH->jb; l <= FORLIM1; l++)
      chain[l].sheetlabel = WITH->sheetname;
  }
}  /* Markstrands */


/***/
/*--------------------------------------------------------------------*/

static void Flagbridge()
{
  long i, FORLIM;
  Bridge *WITH;

  /***/

  for (i = 0; i < MAXBRIDGE; i++) {
    WITH = &bridgetable[i];
    WITH->ib = 0;
    WITH->ie = 0;
    WITH->jb = 0;
    WITH->je = 0;
    WITH->btyp = nobridge;
  }
  nbridge = 0;
  FORLIM = lchain;
  for (i = 2; i < FORLIM; i++)
    Testbridge(i);
  if (nbridge <= 0)
    return;
  Extendladder();
  Sheet();
  Markstrands();
}  /* Flagbridge */


/***/

static void Flagsymbol()
{
  /* FLAGS ALPHA HELICES AND TURNS IN SYMBOL COLUMN */
  long i, j, k;
  char cc;
  long nhset[9];
  structure turn;
  int empty;
  long FORLIM;
  Backbone *WITH;

  P_addset(P_expset(nhset, 0), '>');
  P_addset(nhset, 'X');
  FORLIM = lchain - 4;
  for (i = 2; i <= FORLIM; i++) {
    if (P_inset(chain[i - 1].ss[(long)turn4 - (long)symbol], nhset) &
	P_inset(chain[i].ss[(long)turn4 - (long)symbol], nhset)) {
      for (j = i; j <= i + 3; j++)
	chain[j].ss[0] = 'H';
    }
  }
  FORLIM = lchain - 3;
  for (i = 2; i <= FORLIM; i++) {
    if (P_inset(chain[i - 1].ss[(long)turn3 - (long)symbol], nhset) &
	P_inset(chain[i].ss[(long)turn3 - (long)symbol], nhset)) {
      empty = true;
      for (j = i; j <= i + 2; j++) {
	WITH = &chain[j];
	if (WITH->ss[0] != 'G' && WITH->ss[0] != ' ')
	  empty = false;
      }
      if (empty) {
	for (j = i; j <= i + 2; j++)
	  chain[j].ss[0] = 'G';
      }
    }
  }
  FORLIM = lchain - 5;
  for (i = 2; i <= FORLIM; i++) {
    if (P_inset(chain[i - 1].ss[(long)turn5 - (long)symbol], nhset) &
	P_inset(chain[i].ss[(long)turn5 - (long)symbol], nhset)) {
      empty = true;
      for (j = i; j <= i + 4; j++) {
	WITH = &chain[j];
	if (WITH->ss[0] != 'I' && WITH->ss[0] != ' ')
	  empty = false;
      }
      if (empty) {
	for (j = i; j <= i + 4; j++)
	  chain[j].ss[0] = 'I';
      }
    }
  }
  FORLIM = lchain;
  for (i = 2; i < FORLIM; i++) {
    WITH = &chain[i];
    if (WITH->ss[0] == ' ') {
      cc = ' ';
      j = 1;
      for (turn = turn3;
	   (long)turn <= (long)turn5;
	   turn = (structure)((long)turn + 1)) {
	j++;
	for (k = 1; k <= j; k++) {
	  if (i > k) {
	    if (P_inset(chain[i - k].ss[(long)turn - (long)symbol], nhset))
	      cc = 'T';
	  }
	}
      }
      if (cc == ' ')
	cc = WITH->ss[(long)bend - (long)symbol];
      WITH->ss[0] = cc;
    }
  }
}  /* Flagsymbol */


/***/
/*--------------------------------------------------------------------*/

static void Flagturn()
{
  long i, j, k;
  structure turn;
  char cc;
  long FORLIM1;
  Backbone *WITH;

  /***/

  k = 2;
  cc = '2';
  for (turn = turn3; (long)turn <= (long)turn5; turn = (structure)((long)turn + 1)) {
    k++;
    cc++;
    FORLIM1 = lchain - k;
    for (i = 1; i <= FORLIM1; i++) {
      if (Nochainbreak(i, i + k)) {
	if (Testbond(i + k, i)) {
	  chain[i + k].ss[(long)turn - (long)symbol] = '<';
	  for (j = 1; j < k; j++) {
	    WITH = &chain[i + j];
	    if (WITH->ss[(long)turn - (long)symbol] == ' ')
	      WITH->ss[(long)turn - (long)symbol] = cc;
	  }
	  WITH = &chain[i];
	  if (WITH->ss[(long)turn - (long)symbol] == '<')
	    WITH->ss[(long)turn - (long)symbol] = 'X';
	  else
	    WITH->ss[(long)turn - (long)symbol] = '>';
	}
      }
    }
  }
  FORLIM1 = lchain;
  for (i = 1; i <= FORLIM1; i++) {
    WITH = &chain[i];
    if (WITH->kappa != 360.0 && WITH->kappa > 70.0)
      WITH->ss[(long)bend - (long)symbol] = 'S';
  }
  Flagsymbol();
}  /* Flagturn */


/***/
/*--------------------------------------------------------------------*/


/***/

static void Statistics()
{
  long i, j, k, nchain, nres, nhbond, lhelix;
  bridgetyp b;
  char cc;
  double Surface;
  long nhbturn[11];
  bridgeset ladderset;
  long hbridge[(long)antiparallel - (long)parallel + 1];
  long helixhist[MAXHIST], sheethist[MAXHIST];
  long betahist[(long)antiparallel - (long)parallel + 1][MAXHIST];
  long FORLIM, FORLIM1;
  Backbone *WITH;
  Bridge *WITH1;
  long SET1[257];

  lhelix = 0;
  nhbond = 0;
  nchain = 0;
  nres = 0;
  for (i = 0; i < MAXHIST; i++) {
    for (b = parallel; (long)b <= (long)antiparallel; b = (bridgetyp)((long)b + 1))
      betahist[(long)b - (long)parallel][i] = 0;
    helixhist[i] = 0;
    sheethist[i] = 0;
  }
  Surface = 0.0;
  for (k = 0; k <= 10; k++)
    nhbturn[k] = 0;
  for (b = parallel; (long)b <= (long)antiparallel; b = (bridgetyp)((long)b + 1))
    hbridge[(long)b - (long)parallel] = 0;
  FORLIM = lchain;
  for (i = 0; i <= FORLIM; i++) {
    WITH = &chain[i];
    if (Nochainbreak(i, i)) {
      nres++;
      Surface += WITH->access;
      for (j = 0; j <= 1; j++) {
	if (WITH->donor[j].energy < HBHIGH) {
	  nhbond++;
	  k = WITH->donor[j].residue - i;
	  if (labs(k) < 6)
	    nhbturn[k + 5]++;
	}
      }
    } else
      nchain++;
    if (WITH->ss[0] == 'H')
      lhelix++;
    else if (lhelix > 0) {
      if (lhelix > MAXHIST)
	lhelix = MAXHIST;
      helixhist[lhelix - 1]++;
      lhelix = 0;
    }
  }
  if (nbridge > 0) {
    FORLIM = nbridge;
    for (i = 1; i <= FORLIM; i++) {
      WITH1 = &bridgetable[i - 1];
      hbridge[(long)WITH1->btyp - (long)parallel] += WITH1->ie - WITH1->ib + 2;
      if (WITH1->from == 0) {
	j = i;
	k = 0;
	do {
	  k += bridgetable[j - 1].ie - bridgetable[j - 1].ib + 1;
	  j = bridgetable[j - 1].towards;
	} while (j != 0);
	if (k > MAXHIST)
	  k = MAXHIST;
	betahist[(long)WITH1->btyp - (long)parallel][k - 1]++;
      }
    }
  }
  if (nbridge > 0) {
    P_expset(ladderset, 0);
    FORLIM = nbridge;
    for (i = 1; i <= FORLIM; i++)
      P_addset(ladderset, i);
    FORLIM = nbridge;
    for (i = 1; i <= FORLIM; i++) {
      WITH1 = &bridgetable[i - 1];
      if ((WITH1->from == 0) & P_inset(i, ladderset)) {
	if (!P_setequal(P_addset(P_expset(SET1, 0), i), WITH1->linkset) ||
	    WITH1->ie > WITH1->ib) {
	  k = 0;
	  FORLIM1 = nbridge;
	  for (j = 1; j <= FORLIM1; j++) {
	    if ((bridgetable[j - 1].from == 0) & P_inset(j, WITH1->linkset))
	      k++;
	  }
	  sheethist[k - 1]++;
	}
	P_setdiff(ladderset, ladderset, WITH1->linkset);
      }
    }
  }
  fprintf(tapeout,
    "%5ld%3ld%3ld%3ld%3ld TOTAL NUMBER OF RESIDUES, NUMBER OF CHAINS, NUMBER OF SS-BRIDGES(TOTAL,INTRACHAIN,INTERCHAIN)                .\n",
    nres, nchain, nssinter + nssintra, nssintra, nssinter);
  fprintf(tapeout,
    "%8.1f   ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)                                                                         .\n",
    Surface);
  fprintf(tapeout,
    "%5ld%5.1f   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(J)  , SAME NUMBER PER 100 RESIDUES                              .\n",
    nhbond, 100.0 * nhbond / nres);
  i = hbridge[0];
  j = hbridge[(long)antiparallel - (long)parallel];
  fprintf(tapeout,
    "%5ld%5.1f   TOTAL NUMBER OF HYDROGEN BONDS IN     PARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .\n",
    i, 100.0 * i / nres);
  fprintf(tapeout,
    "%5ld%5.1f   TOTAL NUMBER OF HYDROGEN BONDS IN ANTIPARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .\n",
    j, 100.0 * j / nres);
  for (i = -5; i <= 5; i++) {
    if (i < 0)
      cc = '-';
    else
      cc = '+';
    k = labs(i);
    fprintf(tapeout,
      "%5ld%5.1f   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I%c%ld), SAME NUMBER PER 100 RESIDUES                              .\n",
      nhbturn[i + 5], 100.0 * nhbturn[i + 5] / nres, cc, k);
  }
  for (i = 1; i <= MAXHIST; i++)
    fprintf(tapeout, "%3ld", i);
  fprintf(tapeout, "     *** HISTOGRAMS OF ***           .\n");
  for (i = 0; i < MAXHIST; i++)
    fprintf(tapeout, "%3ld", helixhist[i]);
  fprintf(tapeout, "    RESIDUES PER ALPHA HELIX         .\n");
  for (i = 0; i < MAXHIST; i++)
    fprintf(tapeout, "%3ld", betahist[0][i]);
  fprintf(tapeout, "    PARALLEL BRIDGES PER LADDER      .\n");
  for (i = 0; i < MAXHIST; i++)
    fprintf(tapeout, "%3ld", betahist[(long)antiparallel - (long)parallel][i]);
  fprintf(tapeout, "    ANTIPARALLEL BRIDGES PER LADDER  .\n");
  for (i = 0; i < MAXHIST; i++)
    fprintf(tapeout, "%3ld", sheethist[i]);
  fprintf(tapeout, "    LADDERS PER SHEET                .\n");
}  /* Statistics */

/***/
static void WriteClassicHB(Backbone*resinfo, long i, HydrogenBond hb)
{
  double e;
  int limit=-999;
  if (hb.residue != 0)
    hb.residue -= i;
  e = hb.energy / 1000.0;
  if (hb.residue > -1000 && hb.residue < 10000)
    fprintf(tapeout, "%4ld,%4.1f", hb.residue, e);
  else { /*MS cases where the formated output would be corrupted */
    fprintf(tapeout, "****,%4.1f", e);

    fprintf(stderr," !!! Residue ");
    Writeresidue(*resinfo);
    if(hb.residue>0)
	limit=9999;
    fprintf(stderr," has an hbond partner offset (%ld) exceeding %d that cannot\n", 
	limit,hb.residue);
    fprintf(stderr,"     be printed in classic DSSP format. **** will be printed\n"); 
    fprintf(stderr,"     instead. To avoid this problem, do not use the -c option."); 
    fprintf(stderr," !!!\n\n");
  }
}  /* WriteHB */


static void WriteHB(Backbone*resinfo, long i, HydrogenBond hb)
{
  double e;

  if (hb.residue != 0)
    hb.residue -= i;
  e = hb.energy / 1000.0;
  if (hb.residue > -100000 && hb.residue < 1000000)
    fprintf(tapeout, "%6ld,%4.1f", hb.residue, e);
  else { /*MS cases where the formated output would be corrupted */
    fprintf(tapeout, "******,%4.1f", e);

    fprintf(stderr," !!! Residue ");
    Writeresidue(*resinfo);
    fprintf(stderr," has a to large hbond distance (%ld) to be printed.\n", 
	hb.residue);
    fprintf(stderr,"     To not screw the format, ****** will be printed instead."); 
    fprintf(stderr," !!!\n\n");
  }

  /*   write out **** instead of a too long number... */
}  /* Writehb */


static void Writehb(Backbone*resinfo, long i, HydrogenBond hb)
{
    if(gClassicFormatFlag)
	WriteClassicHB(resinfo,i,hb);
    else
	WriteHB(resinfo,i,hb);
}  /* Writehb */

/***/
/*--------------------------------------------------------------------*/

static void Printout()
{
  long i, j;
  structure s;
  double phi, psi, tco;
  long FORLIM;
  Backbone *WITH;
  char pdb_chain_break;
  /***/

  Statistics();
  if(gClassicFormatFlag) {
      fprintf(tapeout,
        "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC   N-H-->O  O-->H-N  N-H-->O  O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA \n");
  } else {
      fprintf(tapeout,
        "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA \n");
  }
  FORLIM = lchain;

  /* A pdb_chain_break is set if a discontinuity in pdb chain name is detected.
     The pdb_chain_break is flagged by the character '*' in the column that
     follows the one-letter AA code. As a result
       coordinate chain break plus pdb chain break: '!*'
       coordinate chain break only                  '! '
       pdb chain break only (very rare !!)          'C*' 
	   where C is the last chain name 
	   MS and CS 8 July 95 */
  for (i = 1; i <= FORLIM; i++) {
    WITH = &chain[i];
    pdb_chain_break=' '; /* default: no pdb chain break */
    if(i>1 && i<=FORLIM-1 && chain[i-1].aa!='!') { /* check limits */
	int next=i;
	if(WITH->aa=='!') {
	    /* if chain break, check pdb chain name of next residue because
	    chain breaks do not have a chain name! */
	    next=i+1;
	}
	if(!gClassicFormatFlag && chain[i-1].aaident[5]!=chain[next].aaident[5]) {
	    pdb_chain_break='*'; /* pdb chain break detected */
	    if(WITH->aa!='!') {
		fprintf(stderr,"!!! A PDB chain break was found at residue ");
		Writeresidue(*WITH);
		fprintf(stderr,"\n    although the coordinates appear continuous !!!\n");

	    }
	}
    }
    // MAIN OUTPUT!
    fprintf(tapeout, "%5ld ", i);
    for (j = 0; j <= 5; j++)
      putc(WITH->aaident[j], tapeout);
    fprintf(tapeout, " %c%c %c ", WITH->aa, pdb_chain_break, WITH->ss[0]);
    for (s = turn3; (long)s <= (long)beta2; s = (structure)((long)s + 1))
      putc(WITH->ss[(long)s - (long)symbol], tapeout);
    for (s = beta1; (long)s <= (long)beta2; s = (structure)((long)s + 1))
      fprintf(tapeout, "%4ld", WITH->partner[(long)s - (long)beta1]);
    fprintf(tapeout, "%c%4ld ", WITH->sheetlabel, WITH->access);
    for (j = 0; j <= 1; j++) {
      Writehb(WITH, i, WITH->acceptor[j]);
      Writehb(WITH, i, WITH->donor[j]);
    }
    phi = 360.0;
    psi = 360.0;
    tco = 0.0;
    if (Nochainbreak(i - 1, i)) {
      phi = Dihedralangle(chain[i - 1].c, WITH->n, WITH->ca, WITH->c);
      tco = Cosangle(WITH->c, WITH->o, chain[i - 1].c, chain[i - 1].o);
    }
    if (Nochainbreak(i, i + 1))
      psi = Dihedralangle(WITH->n, WITH->ca, WITH->c, chain[i + 1].n);
    fprintf(tapeout, "%8.3f%6.1f%6.1f%6.1f%6.1f%7.1f%7.1f%7.1f\n",
	    tco, WITH->kappa, WITH->alpha, phi, psi, WITH->ca[0], WITH->ca[1],
	    WITH->ca[2]);
  }
}  /* Printout */

static void Copyright()
{
    printf("                           DSSP\n");
    printf("            by Wolfgang Kabsch and Chris Sander\n");
    printf("\n");
    printf("Defines secondary structure and solvent exposure of proteins from\n");
    printf("atomic coordinates as given in Protein Data Bank format.\n");
    printf("\n");
    printf("________________________________________________________________________\n");
    printf("The executable version of DSSP may be downloaded via Internet from\n");
    printf("ftp://ftp.embl-heidelberg.de/pub/software/unix/dssp/ or\n");
    printf("http://www.sander.embl-heidelberg.de/dssp/ \n");
    printf("The executable may be used freely for academic purposes. \n");
    printf("Do not redistribute.\n");
    printf("\n");
    printf("The source code of the DSSP program including sample input and output\n");
    printf("files for dataset 1PPT are available from the authors in exchange for\n");
    printf("a signed academic or commercial license agreement.\n");
    printf("The current version is DSSP-July-95. Original algorithm in 1982 with minor \n");
    printf("changes in 1988. Fast accessibility calculation in 1994 by Michael Scharf. \n");
    printf("Refer to W.Kabsch and C.Sander, Biopolymers 22 (1983) 2577-2637.\n");
    printf("Do report errors if you find any.\n");
    printf("\n");
    printf("Copyright by Wolfgang Kabsch and Chris Sander, 1983, 1985, 1988, 1994, 1995. \n");
    printf("Max Planck Institut fuer Medizinische Forschung\n");
    printf("and EMBL, D-69012 Heidelberg, Germany. Fax: +49-6221-387 306\n");
    printf("\n");
    printf("Email: Sander@embl-heidelberg.de or Sander@embl-ebi.ac.uk\n");
    printf("       Kabsch@embl-heidelberg.de\n");
    printf("\n");
    printf("________________________________________________________________________\n");
    printf("Related databases and software available from the Protein\n");
    printf("Design Group at EMBL-Heidelberg from \n");
    printf("ftp://ftp.embl-heidelberg.de in /pub/databases/protein_extras/\n");
    printf("and /pub/software/unix or from\n");
    printf("http://www.sander.embl-heidelberg.de\n");
    printf("\n");
    printf("pdb_select   Representative set of protein structures.\n");
    printf("             See U. Hobohm, C. Sander, M. Scharf and R. Schneider.\n");
    printf("             Protein Science 1, 409-417 (1992) and \n");
    printf("             U.Hobohm and C.Sander, Protein Science 3, 522-524 (1994).\n");
    printf("DSSP         Dictionary of secondary structures of proteins.\n");
    printf("HSSP         Database of sequence-homology derived protein families.\n");
    printf("             See C. Sander and R.Schneider, Proteins 9, 56-68 (1991)\n");
    printf("FSSP         Database of protein structure families with\n");
    printf("             common folding motifs.\n");
    printf("             See L.Holm, C. Ouzounis, C. Sander, G.Tuparev, G. Vriend\n");
    printf("             See Protein Science 1, 1691-1698 (1992)\n");
    printf("             L.Holm, C.Sander, Nucleic Acids Research 21, 3105-3109 (1994).\n");
    printf("WHAT CHECK   Validation software for protein structure models\n");
    printf("             derived from the WHAT IF modelling package, by G. Vriend\n");
    printf("             and Rob Hooft (vriend or hooft@embl-heidelberg.de).\n");
    printf("\n");
    printf("In the XSSP databases, there is one dataset for each unique or\n");
    printf("             representative PDB protein, e.g., 1PPT.HSSP etc.\n");
    printf("\n");
    printf("************************************************************************\n");
    printf("ftp://ftp.embl-heidelberg.de/pub/software/unix/dssp/ or\n");
    printf("http://www.sander.embl-heidelberg.de/dssp/ \n");
    printf("************************************************************************\n");
}

static void PrintLicense()
{
    printf("DSSP version July 1995.\n");
/*    printf("This executable is licensed to:\n"); * JL */
/*    printf("%s\n",GetLicenseString()); * JL */
/*    printf("\n"); * JL */
/*    printf("The executable may be used freely for academic purposes.\n"); * JL */
/*    printf("\n"); * JL */
/*    printf("Restrictions:Commercial users must apply for a license.\n"); * JL */
/*    printf("             Not to be used for classified research.\n"); * JL */
/*    printf("             Do not redistribute.\n"); * JL */
/*    printf("Obtain additional copies only from http://www.sander.embl-heidelberg.de/dssp/\n"); * JL */
/*    printf("\n"); * JL */
/*    printf("An academic license for the DSSP source code"); * JL */
    printf("((c) W. Kabsch, C. Sander and MPI-MF, 1983, 1985, 1988, 1994, 1995)\n");
    printf("is available in exchange for the following commitments:\n");
    printf("\n");
    printf("I hereby certify that\n");
    printf("\n");
    printf("        (1) I am an academic user at an academic research institution. In\n");
    printf("            using the software, I will respect the interests of the authors\n");
    printf("            and their institutions.\n");
    printf("\n");
    printf("        (2) I will not use the software in commercial activities without\n");
    printf("            a written commercial license agreement; commercial activities\n");
    printf("            include, in particular, work under contract from a commercial\n");
    printf("            company.\n");
    printf("\n");
    printf("        (3) I will not redistribute the software to others outside of my\n");
    printf("            immediate research group. I will suggest to other interested\n");
    printf("            research groups to contact the authors directly.\n");
    printf("\n");
    printf("        (4) I will not alter or suppress the run-time copyright message.\n");
    printf("\n");
    printf("        (5) I will acknowledge the program authors on any publication of\n");
    printf("            scientific results based in part on use of the program and\n");
    printf("            cite the article in which the program was described.\n");
    printf("\n");
    printf("        (6) I will report evidence of program bugs to the authors.\n");
    printf("\n");
    printf("        (7) I will send the source code of any bug corrections and program\n");
    printf("            extensions, major or minor, to the original authors, for free\n");
    printf("            academic use. If I have made major extensions which are incor-\n");
    printf("            porated by the authors, I reserve the right to be appropriately\n");
    printf("            included in any future commercial license agreement.\n");
    printf("\n");
    printf("        (8) I will not extract part of the software, e.g. modules or sub-\n");
    printf("            routines, for use in other contexts without permission by the\n");
    printf("            authors.\n");
    printf("\n");
    printf("        (9) I will not use the program in the context of classified research.\n");
    printf("\n");
/*    printf("To apply for a source code license send an e-mail message to\n"); * JL */
/*    printf("Sander@embl-heidelberg.de or use the URL"); * JL */
/*    printf("http://www.sander.embl-heidelberg.de/dssp/\n"); * JL */
}

static void Usage(char*progname,FILE* fp)
{
  fprintf(fp,"COPYRIGHT\n");
  fprintf(fp,"  W. Kabsch, C. Sander and MPI-MF, 1983, 1985, 1988, 1994 1995\n");
  fprintf(fp,"USAGE\n");
  fprintf(fp,"  %s [-na] [-v] pdb_file [dssp_file]\n",progname);
  fprintf(fp,"  %s [-na] [-v] -- [dssp_file]\n",progname);
  fprintf(fp,"  %s [-h] [-?] [-V]\n",progname);
  fprintf(fp,"OPTIONS\n");
  fprintf(fp,"  -na   Disables the calculation of accessible surface.\n");
  fprintf(fp,"  -c    Classic (old) format.\n");
  fprintf(fp,"  -v    Verbose.\n");
  fprintf(fp,"  --    Read from standard input.\n");
  fprintf(fp,"\n");
  fprintf(fp,"  -h -? Prints a help message.\n");
/*  fprintf(fp,"  -l    Prints the license information.\n"); * JL */
  fprintf(fp,"  -V    Prints version, as in first line of the output.\n");
/*  fprintf(fp,"LICENSED TO\n"); * JL */
/*  fprintf(fp,"  %s\n",GetLicenseString()); * JL */
}

static void Help(char*progname,FILE* fp)
{

}/***/
static void P(char*progname,FILE* fp)
{

}/***/
/*--------------------------------------------------------------------*/

int
main(int argc, char *argv[])
{
  int noaccFlag,verboseFlag,iarg;
  int print_header=0;
  int license = 1;
  noaccFlag=0;
  verboseFlag=0;
  gClassicFormatFlag=0;
  tapein=NULL;
  tapeout=NULL;
  /* license=HasLicense(); * JL */
  if(argc==1)
    {
      Usage(argv[0],stdout);
      exit(-1);
    }
  for(iarg=1;iarg<argc && *argv[iarg]=='-';iarg++)
    {
	if(strcmp("-na",argv[iarg])==0) {
	    noaccFlag=1;
	} else if(strcmp("-c",argv[iarg])==0) {
	    gClassicFormatFlag=1;
	    gDSSP_Version=DSSP_VERSION_OLD;
	} else if(strcmp("-v",argv[iarg])==0) {
	    verboseFlag=1;
	} else if(strcmp("-V",argv[iarg])==0) {
	    print_header=1;
	} else if(strcmp("-l",argv[iarg])==0) {
	    PrintLicense();
	    exit(0);
	} else if(strcmp("-?",argv[iarg])==0 ||
		  strcmp("-h",argv[iarg])==0) {
	    Copyright();
	    /*Usage(argv[0],stdout);*/
	    exit(0);
	} else if(strcmp("--",argv[iarg])==0) {
	    tapein=stdin;
	} else {
	    fprintf(stderr," !!! Unknown argument '%s'!\n\n",argv[iarg]);
	    Usage(argv[0],stderr);
	    exit(-1);
	}
    }
  if(print_header) {
    printf("%s\n",gDSSP_Version);
    exit(0);
  }
  if(tapein==NULL && iarg<argc) {
      if((tapein=fopen(argv[iarg],"r"))==NULL) {
	  fprintf(stderr,"Can't open file <%s> for reading!\n",argv[iarg]);
	  exit(-1);
      }
      iarg++;
  }
  if(iarg<argc) {
      if((tapeout=fopen(argv[iarg],"w"))==NULL) {
	  fprintf(stderr,"Can't open file <%s> for writing!\n",argv[iarg]);
	  exit(-1);
      }
  } else {
      tapeout=stdout;
  }
  if(tapein==NULL || tapeout==NULL) {
      Usage(argv[0],stderr);
      exit(-1);
  }
#ifdef TIMELOCK
  if(time(0)>789397618)
    {
      fprintf(stderr,"*************************************************************************\n");
      fprintf(stderr,"*           This version of dssp has expired !!!!!                      *\n");
      fprintf(stderr,"*************************************************************************\n");
      exit(-1)
    }
#endif
  lchain = 0;
  if(license) {
      Inputcoordinates(&lchain);
      if (!Nochainbreak(1, lchain))
	if(verboseFlag) 
	    fprintf(stderr," !!! Polypeptide chain interrupted !!!\n");
      if(verboseFlag) 
	  fprintf(stderr,"Inputcoordinates done%12ld\n", lchain);
      Flagssbonds();
      if(verboseFlag) 
	  fprintf(stderr,"Flagssbonds done\n");
      Flagchirality();
      if(verboseFlag) 
	  fprintf(stderr,"Flagchirality done\n");
      FlagHydrogenBonds();
      if(verboseFlag) 
	  fprintf(stderr,"Flaghydrogenbonds done\n");
      Flagbridge();
      if(verboseFlag) 
	  fprintf(stderr,"Flagbridge done\n");
      Flagturn();
      if(verboseFlag) 
	  fprintf(stderr,"Flagturn done\n");
      if(noaccFlag==1)
	{
	  if(verboseFlag) fprintf(stderr,"*Accessible surface *NOT* calculated*\n");
	}
      else
	{
      Flagaccess(ORDER);
      if(verboseFlag) 
	  fprintf(stderr,"Flagaccess done\n");
	}
      Printout();
      if(verboseFlag) 
	  fprintf(stderr,"Printout done\n");
      if (tapein != stdin)
	fclose(tapein);
      if (tapeout != stdout)
	fclose(tapeout);
  } else {
      fprintf(stderr,"***This executable is not licensed!!!!***\n");
  }
  exit(0);
}  /* END OF PROGRAM DSSP */

/* sample input file, entry 1PPT from Brookhaven data base, cut after this line.
HEADER    PANCREATIC HORMONE                      16-JAN-81   1PPT      1PPT   3
COMPND    AVIAN PANCREATIC POLYPEPTIDE                                  1PPT   4
SOURCE    TURKEY (MELEAGRIS GALLOPAVO) PANCREAS                         1PPT   5
AUTHOR    T.L.BLUNDELL,J.E.PITTS,I.J.TICKLE,S.P.WOOD                    1PPT   6
ATOM      1  N   GLY     1       2.296  -9.636  18.253  1.00  0.00      1PPT  65
ATOM      2  CA  GLY     1       1.470  -9.017  17.255  1.00  0.00      1PPT  66
ATOM      3  C   GLY     1        .448  -9.983  16.703  1.00  0.00      1PPT  67
ATOM      4  O   GLY     1        .208 -11.066  17.345  1.00  0.00      1PPT  68
ATOM      5  N   PRO     2       -.170  -9.672  15.624  1.00  0.00      1PPT  69
ATOM      6  CA  PRO     2      -1.135 -10.606  14.958  1.00  0.00      1PPT  70
ATOM      7  C   PRO     2       -.376 -11.824  14.490  1.00  0.00      1PPT  71
ATOM      8  O   PRO     2        .776 -11.860  14.075  1.00  0.00      1PPT  72
ATOM      9  CB  PRO     2      -1.717  -9.829  13.776  1.00  0.00      1PPT  73
ATOM     10  CG  PRO     2       -.817  -8.685  13.546  1.00  0.00      1PPT  74
ATOM     11  CD  PRO     2        .108  -8.454  14.780  1.00  0.00      1PPT  75
ATOM     12  N   SER     3      -1.184 -12.918  14.566  1.00  0.00      1PPT  76
ATOM     13  CA  SER     3       -.626 -14.187  14.053  1.00  0.00      1PPT  77
ATOM     14  C   SER     3       -.642 -14.190  12.493  1.00  0.00      1PPT  78
ATOM     15  O   SER     3      -1.149 -13.332  11.830  1.00  0.00      1PPT  79
ATOM     16  CB  SER     3      -1.360 -15.359  14.573  1.00  0.00      1PPT  80
ATOM     17  OG  SER     3      -2.655 -15.234  14.212  1.00  0.00      1PPT  81
ATOM     18  N   GLN     4        .243 -14.995  11.964  1.00  0.00      1PPT  82
ATOM     19  CA  GLN     4        .489 -14.940  10.481  1.00  0.00      1PPT  83
ATOM     20  C   GLN     4       -.766 -15.384   9.734  1.00  0.00      1PPT  84
ATOM     21  O   GLN     4      -1.330 -16.452  10.019  1.00  0.00      1PPT  85
ATOM     22  CB  GLN     4       1.639 -15.895  10.114  1.00  0.00      1PPT  86
ATOM     23  CG  GLN     4       2.182 -15.697   8.704  1.00  0.00      1PPT  87
ATOM     24  CD  GLN     4       3.315 -16.670   8.366  1.00  0.00      1PPT  88
ATOM     25  OE1 GLN     4       3.718 -16.761   7.207  1.00  0.00      1PPT  89
ATOM     26  NE2 GLN     4       3.864 -17.403   9.317  1.00  0.00      1PPT  90
ATOM     27  N   PRO     5      -1.196 -14.647   8.750  1.00  0.00      1PPT  91
ATOM     28  CA  PRO     5      -2.414 -14.970   8.087  1.00  0.00      1PPT  92
ATOM     29  C   PRO     5      -2.264 -16.297   7.258  1.00  0.00      1PPT  93
ATOM     30  O   PRO     5      -1.184 -16.595   6.819  1.00  0.00      1PPT  94
ATOM     31  CB  PRO     5      -2.798 -13.854   7.153  1.00  0.00      1PPT  95
ATOM     32  CG  PRO     5      -1.809 -12.748   7.438  1.00  0.00      1PPT  96
ATOM     33  CD  PRO     5       -.768 -13.190   8.408  1.00  0.00      1PPT  97
ATOM     34  N   THR     6      -3.381 -16.917   7.174  1.00  0.00      1PPT  98
ATOM     35  CA  THR     6      -3.548 -18.158   6.308  1.00  0.00      1PPT  99
ATOM     36  C   THR     6      -3.745 -17.747   4.861  1.00  0.00      1PPT 100
ATOM     37  O   THR     6      -4.693 -17.045   4.518  1.00  0.00      1PPT 101
ATOM     38  CB  THR     6      -4.752 -18.911   6.884  1.00  0.00      1PPT 102
ATOM     39  OG1 THR     6      -4.040 -19.502   8.074  1.00  0.00      1PPT 103
ATOM     40  CG2 THR     6      -4.799 -20.260   6.058  1.00  0.00      1PPT 104
ATOM     41  N   TYR     7      -2.893 -18.207   3.953  1.00  0.00      1PPT 105
ATOM     42  CA  TYR     7      -3.065 -18.017   2.495  1.00  0.00      1PPT 106
ATOM     43  C   TYR     7      -4.327 -18.738   2.010  1.00  0.00      1PPT 107
ATOM     44  O   TYR     7      -4.536 -19.927   2.291  1.00  0.00      1PPT 108
ATOM     45  CB  TYR     7      -1.828 -18.587   1.791  1.00  0.00      1PPT 109
ATOM     46  CG  TYR     7      -1.913 -18.407    .265  1.00  0.00      1PPT 110
ATOM     47  CD1 TYR     7      -2.029 -17.122   -.283  1.00  0.00      1PPT 111
ATOM     48  CD2 TYR     7      -1.884 -19.519   -.588  1.00  0.00      1PPT 112
ATOM     49  CE1 TYR     7      -2.090 -16.948  -1.671  1.00  0.00      1PPT 113
ATOM     50  CE2 TYR     7      -1.943 -19.344  -1.978  1.00  0.00      1PPT 114
ATOM     51  CZ  TYR     7      -2.039 -18.057  -2.521  1.00  0.00      1PPT 115
ATOM     52  OH  TYR     7      -2.067 -17.876  -3.868  1.00  0.00      1PPT 116
ATOM     53  N   PRO     8      -5.261 -18.068   1.439  1.00  0.00      1PPT 117
ATOM     54  CA  PRO     8      -6.566 -18.626    .996  1.00  0.00      1PPT 118
ATOM     55  C   PRO     8      -6.492 -19.530   -.193  1.00  0.00      1PPT 119
ATOM     56  O   PRO     8      -7.584 -20.240   -.510  1.00  0.00      1PPT 120
ATOM     57  CB  PRO     8      -7.488 -17.397    .798  1.00  0.00      1PPT 121
ATOM     58  CG  PRO     8      -6.583 -16.313    .428  1.00  0.00      1PPT 122
ATOM     59  CD  PRO     8      -5.230 -16.608   1.173  1.00  0.00      1PPT 123
ATOM     60  N   GLY     9      -5.375 -19.740   -.857  1.00  0.00      1PPT 124
ATOM     61  CA  GLY     9      -5.423 -20.730  -1.983  1.00  0.00      1PPT 125
ATOM     62  C   GLY     9      -5.256 -19.861  -3.214  1.00  0.00      1PPT 126
ATOM     63  O   GLY     9      -5.790 -18.783  -3.338  1.00  0.00      1PPT 127
ATOM     64  N   ASP    10      -4.626 -20.463  -4.248  1.00  0.00      1PPT 128
ATOM     65  CA  ASP    10      -4.521 -19.865  -5.611  1.00  0.00      1PPT 129
ATOM     66  C   ASP    10      -5.891 -19.644  -6.232  1.00  0.00      1PPT 130
ATOM     67  O   ASP    10      -6.079 -18.696  -7.006  1.00  0.00      1PPT 131
ATOM     68  CB  ASP    10      -3.697 -20.772  -6.523  1.00  0.00      1PPT 132
ATOM     69  CG  ASP    10      -2.225 -20.895  -6.117  1.00  0.00      1PPT 133
ATOM     70  OD1 ASP    10      -1.521 -21.886  -6.544  1.00  0.00      1PPT 134
ATOM     71  OD2 ASP    10      -1.682 -20.014  -5.347  1.00  0.00      1PPT 135
ATOM     72  N   ASP    11      -6.810 -20.490  -5.917  1.00  0.00      1PPT 136
ATOM     73  CA  ASP    11      -8.106 -20.411  -6.537  1.00  0.00      1PPT 137
ATOM     74  C   ASP    11      -9.141 -19.681  -5.693  1.00  0.00      1PPT 138
ATOM     75  O   ASP    11     -10.273 -19.451  -6.151  1.00  0.00      1PPT 139
ATOM     76  CB  ASP    11      -8.681 -21.809  -6.852  1.00  0.00      1PPT 140
ATOM     77  CG  ASP    11      -7.791 -22.570  -7.829  1.00  0.00      1PPT 141
ATOM     78  OD1 ASP    11      -7.396 -21.995  -8.913  1.00  0.00      1PPT 142
ATOM     79  OD2 ASP    11      -7.431 -23.778  -7.563  1.00  0.00      1PPT 143
ATOM     80  N   ALA    12      -8.612 -18.887  -4.775  1.00  0.00      1PPT 144
ATOM     81  CA  ALA    12      -9.622 -18.132  -4.017  1.00  0.00      1PPT 145
ATOM     82  C   ALA    12     -10.101 -16.933  -4.820  1.00  0.00      1PPT 146
ATOM     83  O   ALA    12      -9.482 -16.473  -5.779  1.00  0.00      1PPT 147
ATOM     84  CB  ALA    12      -8.829 -17.575  -2.793  1.00  0.00      1PPT 148
ATOM     85  N   PRO    13     -11.366 -16.547  -4.687  1.00  0.00      1PPT 149
ATOM     86  CA  PRO    13     -11.981 -15.406  -5.466  1.00  0.00      1PPT 150
ATOM     87  C   PRO    13     -11.187 -14.121  -5.215  1.00  0.00      1PPT 151
ATOM     88  O   PRO    13     -10.522 -13.958  -4.032  1.00  0.00      1PPT 152
ATOM     89  CB  PRO    13     -13.424 -15.243  -4.908  1.00  0.00      1PPT 153
ATOM     90  CG  PRO    13     -13.500 -16.009  -3.659  1.00  0.00      1PPT 154
ATOM     91  CD  PRO    13     -12.236 -16.882  -3.480  1.00  0.00      1PPT 155
ATOM     92  N   VAL    14     -11.180 -13.190  -6.126  1.00  0.00      1PPT 156
ATOM     93  CA  VAL    14     -10.341 -11.949  -5.914  1.00  0.00      1PPT 157
ATOM     94  C   VAL    14     -10.673 -11.235  -4.610  1.00  0.00      1PPT 158
ATOM     95  O   VAL    14      -9.729 -10.720  -3.902  1.00  0.00      1PPT 159
ATOM     96  CB  VAL    14     -10.477 -11.110  -7.162  1.00  0.00      1PPT 160
ATOM     97  CG1 VAL    14      -9.809  -9.750  -7.062  1.00  0.00      1PPT 161
ATOM     98  CG2 VAL    14     -10.013 -11.873  -8.431  1.00  0.00      1PPT 162
ATOM     99  N   GLU    15     -11.842 -11.320  -4.113  1.00  0.00      1PPT 163
ATOM    100  CA  GLU    15     -12.142 -10.553  -2.855  1.00  0.00      1PPT 164
ATOM    101  C   GLU    15     -11.451 -11.174  -1.704  1.00  0.00      1PPT 165
ATOM    102  O   GLU    15     -11.170 -10.380   -.760  1.00  0.00      1PPT 166
ATOM    103  CB  GLU    15     -13.711 -10.553  -2.589  1.00  0.00      1PPT 167
ATOM    104  CG  GLU    15     -14.210 -11.854  -1.955  1.00  0.00      1PPT 168
ATOM    105  CD  GLU    15     -15.561 -11.698  -1.254  1.00  0.00      1PPT 169
ATOM    106  OE1 GLU    15     -16.152 -12.731   -.757  1.00  0.00      1PPT 170
ATOM    107  OE2 GLU    15     -16.108 -10.534  -1.161  1.00  0.00      1PPT 171
ATOM    108  N   ASP    16     -11.319 -12.496  -1.702  1.00  0.00      1PPT 172
ATOM    109  CA  ASP    16     -10.488 -13.190   -.722  1.00  0.00      1PPT 173
ATOM    110  C   ASP    16      -9.000 -12.982   -.958  1.00  0.00      1PPT 174
ATOM    111  O   ASP    16      -8.238 -12.840    .013  1.00  0.00      1PPT 175
ATOM    112  CB  ASP    16     -10.706 -14.670   -.638  1.00  0.00      1PPT 176
ATOM    113  CG  ASP    16     -12.106 -15.022   -.152  1.00  0.00      1PPT 177
ATOM    114  OD1 ASP    16     -12.571 -16.208   -.343  1.00  0.00      1PPT 178
ATOM    115  OD2 ASP    16     -12.821 -14.131    .443  1.00  0.00      1PPT 179
ATOM    116  N   LEU    17      -8.476 -12.788  -2.105  1.00  0.00      1PPT 180
ATOM    117  CA  LEU    17      -7.028 -12.438  -2.126  1.00  0.00      1PPT 181
ATOM    118  C   LEU    17      -6.810 -10.983  -1.717  1.00  0.00      1PPT 182
ATOM    119  O   LEU    17      -5.812 -10.718  -1.159  1.00  0.00      1PPT 183
ATOM    120  CB  LEU    17      -6.647 -12.470  -3.630  1.00  0.00      1PPT 184
ATOM    121  CG  LEU    17      -6.525 -13.931  -4.064  1.00  0.00      1PPT 185
ATOM    122  CD1 LEU    17      -6.250 -13.916  -5.582  1.00  0.00      1PPT 186
ATOM    123  CD2 LEU    17      -5.372 -14.656  -3.263  1.00  0.00      1PPT 187
ATOM    124  N   ILE    18      -7.786 -10.075  -1.877  1.00  0.00      1PPT 188
ATOM    125  CA  ILE    18      -7.676  -8.682  -1.354  1.00  0.00      1PPT 189
ATOM    126  C   ILE    18      -7.789  -8.754    .184  1.00  0.00      1PPT 190
ATOM    127  O   ILE    18      -6.937  -8.096    .818  1.00  0.00      1PPT 191
ATOM    128  CB  ILE    18      -8.822  -7.866  -1.960  1.00  0.00      1PPT 192
ATOM    129  CG1 ILE    18      -8.514  -7.564  -3.452  1.00  0.00      1PPT 193
ATOM    130  CG2 ILE    18      -8.973  -6.595  -1.116  1.00  0.00      1PPT 194
ATOM    131  CD1 ILE    18      -9.640  -6.945  -4.117  1.00  0.00      1PPT 195
ATOM    132  N   ARG    19      -8.698  -9.540    .771  1.00  0.00      1PPT 196
ATOM    133  CA  ARG    19      -8.636  -9.644   2.284  1.00  0.00      1PPT 197
ATOM    134  C   ARG    19      -7.271 -10.180   2.719  1.00  0.00      1PPT 198
ATOM    135  O   ARG    19      -6.742  -9.797   3.773  1.00  0.00      1PPT 199
ATOM    136  CB  ARG    19      -9.736 -10.584   2.604  1.00  0.00      1PPT 200
ATOM    137  CG  ARG    19     -11.010  -9.860   3.115  1.00  0.00      1PPT 201
ATOM    138  CD  ARG    19     -12.209 -10.742   3.604  1.00  0.00      1PPT 202
ATOM    139  NE  ARG    19     -13.244 -10.589   2.620  1.00  0.00      1PPT 203
ATOM    140  CZ  ARG    19     -14.562 -10.256   2.523  1.00  0.00      1PPT 204
ATOM    141  NH1 ARG    19     -15.649  -9.931   3.450  1.00  0.00      1PPT 205
ATOM    142  NH2 ARG    19     -14.877 -10.167   1.331  1.00  0.00      1PPT 206
ATOM    143  N   PHE    20      -6.704 -11.222   2.134  1.00  0.00      1PPT 207
ATOM    144  CA  PHE    20      -5.441 -11.780   2.548  1.00  0.00      1PPT 208
ATOM    145  C   PHE    20      -4.374 -10.694   2.394  1.00  0.00      1PPT 209
ATOM    146  O   PHE    20      -3.489 -10.541   3.246  1.00  0.00      1PPT 210
ATOM    147  CB  PHE    20      -5.101 -12.968   1.646  1.00  0.00      1PPT 211
ATOM    148  CG  PHE    20      -3.710 -13.511   1.946  1.00  0.00      1PPT 212
ATOM    149  CD1 PHE    20      -2.666 -13.334   1.030  1.00  0.00      1PPT 213
ATOM    150  CD2 PHE    20      -3.489 -14.185   3.148  1.00  0.00      1PPT 214
ATOM    151  CE1 PHE    20      -1.392 -13.833   1.327  1.00  0.00      1PPT 215
ATOM    152  CE2 PHE    20      -2.218 -14.680   3.443  1.00  0.00      1PPT 216
ATOM    153  CZ  PHE    20      -1.168 -14.502   2.535  1.00  0.00      1PPT 217
ATOM    154  N   TYR    21      -4.347  -9.918   1.280  1.00  0.00      1PPT 218
ATOM    155  CA  TYR    21      -3.402  -8.885   1.007  1.00  0.00      1PPT 219
ATOM    156  C   TYR    21      -3.411  -7.902   2.181  1.00  0.00      1PPT 220
ATOM    157  O   TYR    21      -2.361  -7.533   2.718  1.00  0.00      1PPT 221
ATOM    158  CB  TYR    21      -3.821  -8.142   -.271  1.00  0.00      1PPT 222
ATOM    159  CG  TYR    21      -3.056  -6.870   -.555  1.00  0.00      1PPT 223
ATOM    160  CD1 TYR    21      -1.730  -6.925   -.993  1.00  0.00      1PPT 224
ATOM    161  CD2 TYR    21      -3.708  -5.654   -.376  1.00  0.00      1PPT 225
ATOM    162  CE1 TYR    21      -1.042  -5.733  -1.246  1.00  0.00      1PPT 226
ATOM    163  CE2 TYR    21      -3.022  -4.467   -.625  1.00  0.00      1PPT 227
ATOM    164  CZ  TYR    21      -1.692  -4.505  -1.058  1.00  0.00      1PPT 228
ATOM    165  OH  TYR    21      -1.035  -3.335  -1.287  1.00  0.00      1PPT 229
ATOM    166  N   ASP    22      -4.642  -7.469   2.585  1.00  0.00      1PPT 230
ATOM    167  CA  ASP    22      -4.728  -6.473   3.688  1.00  0.00      1PPT 231
ATOM    168  C   ASP    22      -4.207  -7.089   4.978  1.00  0.00      1PPT 232
ATOM    169  O   ASP    22      -3.568  -6.404   5.780  1.00  0.00      1PPT 233
ATOM    170  CB  ASP    22      -6.194  -6.063   3.854  1.00  0.00      1PPT 234
ATOM    171  CG  ASP    22      -6.707  -5.158   2.746  1.00  0.00      1PPT 235
ATOM    172  OD1 ASP    22      -7.967  -5.120   2.494  1.00  0.00      1PPT 236
ATOM    173  OD2 ASP    22      -5.890  -4.431   2.071  1.00  0.00      1PPT 237
ATOM    174  N   ASN    23      -4.505  -8.385   5.271  1.00  0.00      1PPT 238
ATOM    175  CA  ASN    23      -4.038  -8.947   6.542  1.00  0.00      1PPT 239
ATOM    176  C   ASN    23      -2.506  -9.142   6.402  1.00  0.00      1PPT 240
ATOM    177  O   ASN    23      -1.837  -8.942   7.451  1.00  0.00      1PPT 241
ATOM    178  CB  ASN    23      -4.716 -10.254   6.855  1.00  0.00      1PPT 242
ATOM    179  CG  ASN    23      -6.186 -10.108   7.257  1.00  0.00      1PPT 243
ATOM    180  OD1 ASN    23      -6.841 -11.277   7.242  1.00  0.00      1PPT 244
ATOM    181  ND2 ASN    23      -6.691  -9.038   7.581  1.00  0.00      1PPT 245
ATOM    182  N   LEU    24      -1.987  -9.505   5.267  1.00  0.00      1PPT 246
ATOM    183  CA  LEU    24       -.451  -9.679   5.213  1.00  0.00      1PPT 247
ATOM    184  C   LEU    24        .173  -8.275   5.430  1.00  0.00      1PPT 248
ATOM    185  O   LEU    24       1.220  -8.264   5.981  1.00  0.00      1PPT 249
ATOM    186  CB  LEU    24       -.202 -10.148   3.805  1.00  0.00      1PPT 250
ATOM    187  CG  LEU    24       1.277 -10.481   3.548  1.00  0.00      1PPT 251
ATOM    188  CD1 LEU    24       1.764 -11.470   4.539  1.00  0.00      1PPT 252
ATOM    189  CD2 LEU    24       1.499 -11.012   2.076  1.00  0.00      1PPT 253
ATOM    190  N   GLN    25       -.410  -7.204   4.911  1.00  0.00      1PPT 254
ATOM    191  CA  GLN    25        .085  -5.877   5.096  1.00  0.00      1PPT 255
ATOM    192  C   GLN    25        .191  -5.541   6.582  1.00  0.00      1PPT 256
ATOM    193  O   GLN    25       1.265  -5.150   7.065  1.00  0.00      1PPT 257
ATOM    194  CB  GLN    25       -.806  -4.832   4.429  1.00  0.00      1PPT 258
ATOM    195  CG  GLN    25       -.281  -3.402   4.489  1.00  0.00      1PPT 259
ATOM    196  CD  GLN    25       -.921  -2.504   3.422  1.00  0.00      1PPT 260
ATOM    197  OE1 GLN    25       -.397  -1.428   3.134  1.00  0.00      1PPT 261
ATOM    198  NE2 GLN    25      -2.028  -2.888   2.811  1.00  0.00      1PPT 262
ATOM    199  N   GLN    26       -.856  -5.749   7.310  1.00  0.00      1PPT 263
ATOM    200  CA  GLN    26       -.856  -5.563   8.757  1.00  0.00      1PPT 264
ATOM    201  C   GLN    26        .308  -6.289   9.398  1.00  0.00      1PPT 265
ATOM    202  O   GLN    26       1.009  -5.921  10.266  1.00  0.00      1PPT 266
ATOM    203  CB  GLN    26      -2.266  -5.866   9.301  1.00  0.00      1PPT 267
ATOM    204  CG  GLN    26      -2.357  -5.747  10.824  1.00  0.00      1PPT 268
ATOM    205  CD  GLN    26      -2.333  -4.297  11.313  1.00  0.00      1PPT 269
ATOM    206  OE1 GLN    26      -2.414  -4.053  12.516  1.00  0.00      1PPT 270
ATOM    207  NE2 GLN    26      -2.225  -3.309  10.444  1.00  0.00      1PPT 271
ATOM    208  N   TYR    27        .366  -7.626   9.157  1.00  0.00      1PPT 272
ATOM    209  CA  TYR    27       1.356  -8.520   9.698  1.00  0.00      1PPT 273
ATOM    210  C   TYR    27       2.759  -8.036   9.359  1.00  0.00      1PPT 274
ATOM    211  O   TYR    27       3.622  -7.942  10.245  1.00  0.00      1PPT 275
ATOM    212  CB  TYR    27       1.111  -9.949   9.116  1.00  0.00      1PPT 276
ATOM    213  CG  TYR    27       2.122 -10.986   9.594  1.00  0.00      1PPT 277
ATOM    214  CD1 TYR    27       3.124 -11.439   8.729  1.00  0.00      1PPT 278
ATOM    215  CD2 TYR    27       2.044 -11.485  10.899  1.00  0.00      1PPT 279
ATOM    216  CE1 TYR    27       4.061 -12.377   9.173  1.00  0.00      1PPT 280
ATOM    217  CE2 TYR    27       2.984 -12.421  11.346  1.00  0.00      1PPT 281
ATOM    218  CZ  TYR    27       3.998 -12.862  10.486  1.00  0.00      1PPT 282
ATOM    219  OH  TYR    27       4.932 -13.751  10.922  1.00  0.00      1PPT 283
ATOM    220  N   LEU    28       3.051  -7.735   8.171  1.00  0.00      1PPT 284
ATOM    221  CA  LEU    28       4.342  -7.279   7.784  1.00  0.00      1PPT 285
ATOM    222  C   LEU    28       4.675  -5.907   8.478  1.00  0.00      1PPT 286
ATOM    223  O   LEU    28       5.887  -5.811   8.927  1.00  0.00      1PPT 287
ATOM    224  CB  LEU    28       4.614  -7.138   6.312  1.00  0.00      1PPT 288
ATOM    225  CG  LEU    28       4.490  -8.512   5.562  1.00  0.00      1PPT 289
ATOM    226  CD1 LEU    28       4.429  -8.253   3.969  1.00  0.00      1PPT 290
ATOM    227  CD2 LEU    28       5.761  -9.362   5.838  1.00  0.00      1PPT 291
ATOM    228  N   ASN    29       3.737  -4.993   8.530  1.00  0.00      1PPT 292
ATOM    229  CA  ASN    29       4.064  -3.735   9.328  1.00  0.00      1PPT 293
ATOM    230  C   ASN    29       4.502  -4.080  10.751  1.00  0.00      1PPT 294
ATOM    231  O   ASN    29       5.252  -3.321  11.381  1.00  0.00      1PPT 295
ATOM    232  CB  ASN    29       2.748  -2.894   9.387  1.00  0.00      1PPT 296
ATOM    233  CG  ASN    29       2.581  -2.132   8.083  1.00  0.00      1PPT 297
ATOM    234  OD1 ASN    29       1.565  -1.465   7.896  1.00  0.00      1PPT 298
ATOM    235  ND2 ASN    29       3.539  -2.200   7.175  1.00  0.00      1PPT 299
ATOM    236  N   VAL    30       3.812  -5.044  11.456  1.00  0.00      1PPT 300
ATOM    237  CA  VAL    30       4.069  -5.352  12.824  1.00  0.00      1PPT 301
ATOM    238  C   VAL    30       5.379  -6.029  12.954  1.00  0.00      1PPT 302
ATOM    239  O   VAL    30       6.270  -5.680  13.834  1.00  0.00      1PPT 303
ATOM    240  CB  VAL    30       2.961  -6.250  13.460  1.00  0.00      1PPT 304
ATOM    241  CG1 VAL    30       3.526  -6.593  14.954  1.00  0.00      1PPT 305
ATOM    242  CG2 VAL    30       1.683  -5.435  13.674  1.00  0.00      1PPT 306
ATOM    243  N   VAL    31       5.698  -6.965  12.038  1.00  0.00      1PPT 307
ATOM    244  CA  VAL    31       7.001  -7.699  12.162  1.00  0.00      1PPT 308
ATOM    245  C   VAL    31       8.157  -6.782  11.957  1.00  0.00      1PPT 309
ATOM    246  O   VAL    31       9.302  -7.004  12.397  1.00  0.00      1PPT 310
ATOM    247  CB  VAL    31       6.845  -8.810  10.979  1.00  0.00      1PPT 311
ATOM    248  CG1 VAL    31       8.217  -9.440  10.757  1.00  0.00      1PPT 312
ATOM    249  CG2 VAL    31       5.842  -9.913  11.342  1.00  0.00      1PPT 313
ATOM    250  N   THR    32       8.037  -5.850  11.015  1.00  0.00      1PPT 314
ATOM    251  CA  THR    32       9.068  -4.861  10.736  1.00  0.00      1PPT 315
ATOM    252  C   THR    32       8.975  -3.622  11.693  1.00  0.00      1PPT 316
ATOM    253  O   THR    32       9.882  -2.752  11.522  1.00  0.00      1PPT 317
ATOM    254  CB  THR    32       8.833  -4.336   9.342  1.00  0.00      1PPT 318
ATOM    255  OG1 THR    32       7.762  -3.572   9.058  1.00  0.00      1PPT 319
ATOM    256  CG2 THR    32       9.614  -4.848   8.439  1.00  0.00      1PPT 320
ATOM    257  N   ARG    33       8.154  -3.603  12.666  1.00  0.00      1PPT 321
ATOM    258  CA  ARG    33       7.995  -2.499  13.658  1.00  0.00      1PPT 322
ATOM    259  C   ARG    33       7.787  -1.166  12.943  1.00  0.00      1PPT 323
ATOM    260  O   ARG    33       8.043   -.093  13.507  1.00  0.00      1PPT 324
ATOM    261  CB  ARG    33       9.235  -2.379  14.561  1.00  0.00      1PPT 325
ATOM    262  CG  ARG    33       9.711  -3.734  15.065  1.00  0.00      1PPT 326
ATOM    263  CD  ARG    33      10.875  -3.613  16.039  1.00  0.00      1PPT 327
ATOM    264  NE  ARG    33      10.549  -2.811  17.223  1.00  0.00      1PPT 328
ATOM    265  CZ  ARG    33       9.938  -3.306  18.303  1.00  0.00      1PPT 329
ATOM    266  NH1 ARG    33       9.575  -4.596  18.352  1.00  0.00      1PPT 330
ATOM    267  NH2 ARG    33       9.649  -2.584  19.395  1.00  0.00      1PPT 331
ATOM    268  N   HIS    34       7.031  -1.228  11.897  1.00  0.00      1PPT 332
ATOM    269  CA  HIS    34       6.779    .039  11.099  1.00  0.00      1PPT 333
ATOM    270  C   HIS    34       5.289    .163  10.798  1.00  0.00      1PPT 334
ATOM    271  O   HIS    34       4.835   -.137   9.689  1.00  0.00      1PPT 335
ATOM    272  CB  HIS    34       7.587   -.011   9.878  1.00  0.00      1PPT 336
ATOM    273  CG  HIS    34       7.608   1.293   9.098  1.00  0.00      1PPT 337
ATOM    274  ND1 HIS    34       6.953   1.430   7.879  1.00  0.00      1PPT 338
ATOM    275  CD2 HIS    34       8.195   2.486   9.363  1.00  0.00      1PPT 339
ATOM    276  CE1 HIS    34       7.144   2.668   7.454  1.00  0.00      1PPT 340
ATOM    277  NE2 HIS    34       7.884   3.310   8.330  1.00  0.00      1PPT 341
ATOM    278  N   ARG    35       4.524    .544  11.879  1.00  0.00      1PPT 342
ATOM    279  CA  ARG    35       3.108    .616  11.852  1.00  0.00      1PPT 343
ATOM    280  C   ARG    35       2.637   1.882  11.134  1.00  0.00      1PPT 344
ATOM    281  O   ARG    35       1.446   2.229  11.173  1.00  0.00      1PPT 345
ATOM    282  CB  ARG    35       2.550    .636  13.314  1.00  0.00      1PPT 346
ATOM    283  CG  ARG    35       2.848   -.712  13.994  1.00  0.00      1PPT 347
ATOM    284  CD  ARG    35       2.475   -.788  15.476  1.00  0.00      1PPT 348
ATOM    285  NE  ARG    35       3.312  -1.745  16.223  1.00  0.00      1PPT 349
ATOM    286  CZ  ARG    35       2.837  -2.773  16.945  1.00  0.00      1PPT 350
ATOM    287  NH1 ARG    35       1.521  -3.007  17.037  1.00  0.00      1PPT 351
ATOM    288  NH2 ARG    35       3.612  -3.634  17.621  1.00  0.00      1PPT 352
ATOM    289  N   TYR    36       3.365   2.609  10.443  1.00  0.00      1PPT 353
ATOM    290  CA  TYR    36       2.765   3.446   9.296  1.00  0.00      1PPT 354
ATOM    291  C   TYR    36       2.332   2.479   8.197  1.00  0.00      1PPT 355
ATOM    292  O   TYR    36       3.166   1.720   7.671  1.00  0.00      1PPT 356
ATOM    293  CB  TYR    36       4.021   4.330   8.787  1.00  0.00      1PPT 357
ATOM    294  CG  TYR    36       4.734   4.795  10.058  1.00  0.00      1PPT 358
ATOM    295  CD1 TYR    36       5.675   3.963  10.681  1.00  0.00      1PPT 359
ATOM    296  CD2 TYR    36       4.424   6.040  10.616  1.00  0.00      1PPT 360
ATOM    297  CE1 TYR    36       6.332   4.393  11.840  1.00  0.00      1PPT 361
ATOM    298  CE2 TYR    36       5.083   6.471  11.773  1.00  0.00      1PPT 362
ATOM    299  CZ  TYR    36       6.043   5.652  12.379  1.00  0.00      1PPT 363
ATOM    300  OH  TYR    36       6.704   6.088  13.485  1.00  0.00      1PPT 364
ATOM    301  OXT TYR    36       1.276   2.139   7.885  1.00  0.00      1PPT 365
TER     302      TYR    36                                              1PPT 366
END                                                                     1PPT 369
sample input file, entry 1PPT from Brookhaven data base, cut before this line.*/

/* sample output file, from PDB entry 1PPT, cut after this line ---------------
**** SECONDARY STRUCTURE DEFINITION BY THE PROGRAM DSSP, VERSION OCT. 1985 **** MONTH=12 DAY=24 YEAR=2001                      .
REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637                                                              .
HEADER    PANCREATIC HORMONE                      16-JAN-81   1PPT                                                             .
COMPND    AVIAN PANCREATIC POLYPEPTIDE                                                                                         .
SOURCE    TURKEY (MELEAGRIS GALLOPAVO) PANCREAS                                                                                .
AUTHOR    T.L.BLUNDELL,J.E.PITTS,I.J.TICKLE,S.P.WOOD                                                                           .
   36  1  0  0  0 TOTAL NUMBER OF RESIDUES, NUMBER OF CHAINS, NUMBER OF SS-BRIDGES(TOTAL,INTRACHAIN,INTERCHAIN)                .
  3430.0   ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)                                                                         .
   22 61.1   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(J)  , SAME NUMBER PER 100 RESIDUES                              .
    0  0.0   TOTAL NUMBER OF HYDROGEN BONDS IN     PARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .
    0  0.0   TOTAL NUMBER OF HYDROGEN BONDS IN ANTIPARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .
    0  0.0   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-5), SAME NUMBER PER 100 RESIDUES                              .
    0  0.0   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-4), SAME NUMBER PER 100 RESIDUES                              .
    0  0.0   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-3), SAME NUMBER PER 100 RESIDUES                              .
    0  0.0   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-2), SAME NUMBER PER 100 RESIDUES                              .
    0  0.0   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I-1), SAME NUMBER PER 100 RESIDUES                              .
    0  0.0   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+0), SAME NUMBER PER 100 RESIDUES                              .
    0  0.0   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+1), SAME NUMBER PER 100 RESIDUES                              .
    2  5.6   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+2), SAME NUMBER PER 100 RESIDUES                              .
    3  8.3   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+3), SAME NUMBER PER 100 RESIDUES                              .
   16 44.4   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+4), SAME NUMBER PER 100 RESIDUES                              .
    1  2.8   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I+5), SAME NUMBER PER 100 RESIDUES                              .
  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30     *** HISTOGRAMS OF ***           .
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0    RESIDUES PER ALPHA HELIX         .
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    PARALLEL BRIDGES PER LADDER      .
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    ANTIPARALLEL BRIDGES PER LADDER  .
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    LADDERS PER SHEET                .
  #  RESIDUE AA STRUCTURE BP1 BP2  ACC   N-H-->O  O-->H-N  N-H-->O  O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
    1    1   G              0   0  101    0, 0.0  26,-0.1   0, 0.0   2, 0.0   0.000 360.0 360.0 360.0-165.7    1.5   -9.0   17.3
    2    2   P        -     0   0   60    0, 0.0  28, 0.0   0, 0.0   0, 0.0  -0.321 360.0-107.5 -64.4 147.2   -1.1  -10.6   15.0
    3    3   S        -     0   0  106    1,-0.1  24, 0.0  -2, 0.0   0, 0.0  -0.232  48.2 -82.1 -76.1 156.0   -0.6  -14.2   14.1
    4    4   Q        -     0   0  139    1,-0.1  -1,-0.1  23,-0.1  20,-0.1  -0.552  52.4-113.7 -65.3 131.7    0.5  -14.9   10.5
    5    5   P        -     0   0   26    0, 0.0   2,-0.5   0, 0.0  -1,-0.1  -0.279  23.3-114.6 -68.9 150.5   -2.4  -15.0    8.1
    6    6   T        -     0   0  126    1, 0.0   0, 0.0   2, 0.0   0, 0.0  -0.732  28.2-122.6 -80.5 122.2   -3.5  -18.2    6.3
    7    7   Y        -     0   0  121   -2,-0.5   3,-0.1   1,-0.2  10,-0.1  -0.481  18.6-142.4 -64.0 120.9   -3.1  -18.0    2.5
    8    8   P        -     0   0   55    0, 0.0  -1,-0.2   0, 0.0   3,-0.1   0.592  46.7 -93.0 -71.2  -4.7   -6.6  -18.6    1.0
    9    9   G    >   -     0   0   27    1,-0.1   3,-1.3   2,-0.1  -3, 0.0  -0.088  31.8 -90.2 110.8 148.9   -5.4  -20.7   -2.0
   10   10   D  T 3  S+     0   0  128    1,-0.2  -1,-0.1  -3,-0.1   7,-0.1   0.814 125.8  53.1 -61.0 -32.5   -4.5  -19.9   -5.6
   11   11   D  T 3  S+     0   0  165   -3,-0.1  -1,-0.2   2, 0.0  -2,-0.1   0.510  77.9 125.0 -97.6  20.7   -8.1  -20.4   -6.5
   12   12   A  S <  S-     0   0   17   -3,-1.3  -3,-0.1   1,-0.1   5, 0.0  -0.102  71.9 -98.5 -79.6 146.6   -9.6  -18.1   -4.0
   13   13   P    >>  -     0   0   73    0, 0.0   4,-2.2   0, 0.0   3,-1.6  -0.385  38.5-104.7 -57.5 153.1  -12.0  -15.4   -5.5
   14   14   V  H 3> S+     0   0  108    1,-0.3   4,-2.6   2,-0.2   5,-0.2   0.880 120.1  59.6 -55.7 -28.4  -10.3  -11.9   -5.9
   15   15   E  H 3> S+     0   0  111    1,-0.3   4,-1.9   2,-0.2  -1,-0.3   0.824 105.4  49.3 -71.6 -34.4  -12.1  -10.6   -2.9
   16   16   D  H <> S+     0   0   58   -3,-1.6   4,-2.0   2,-0.2  -1,-0.3   0.791 109.0  50.6 -71.8 -30.1  -10.5  -13.2   -0.7
   17   17   L  H  X S+     0   0   62   -4,-2.2   4,-2.1   2,-0.2  -2,-0.2   0.919 110.0  51.7 -76.3 -25.0   -7.0  -12.4   -2.1
   18   18   I  H  X S+     0   0   91   -4,-2.6   4,-2.5  -5,-0.2  -2,-0.2   0.940 110.3  48.7 -71.3 -42.9   -7.7   -8.7   -1.4
   19   19   R  H  X S+     0   0  142   -4,-1.9   4,-2.1   1,-0.2  -1,-0.2   0.885 113.7  46.4 -59.1 -45.9   -8.6   -9.6    2.3
   20   20   F  H  X S+     0   0   39   -4,-2.0   4,-3.2   2,-0.2  -1,-0.2   0.928 110.4  54.1 -60.9 -43.8   -5.4  -11.8    2.5
   21   21   Y  H  X S+     0   0  143   -4,-2.1   4,-2.1   1,-0.2  -2,-0.2   0.943 110.2  46.7 -54.7 -49.0   -3.4   -8.9    1.0
   22   22   D  H  X S+     0   0   79   -4,-2.5   4,-2.0   1,-0.2  -1,-0.2   0.889 113.7  47.4 -64.9 -39.6   -4.7   -6.5    3.7
   23   23   N  H  X S+     0   0   97   -4,-2.1   4,-2.1  -5,-0.2  -1,-0.2   0.905 111.3  50.6 -72.5 -36.0   -4.0   -8.9    6.5
   24   24   L  H  X S+     0   0   58   -4,-3.2   4,-2.6   1,-0.2  -1,-0.2   0.859 108.3  55.0 -65.8 -37.1   -0.5   -9.7    5.2
   25   25   Q  H  X S+     0   0   90   -4,-2.1   4,-2.8   2,-0.2   5,-0.2   0.950 109.1  45.9 -56.8 -51.4    0.1   -5.9    5.1
   26   26   Q  H  X S+     0   0  112   -4,-2.0   4,-1.5   2,-0.2  -1,-0.2   0.903 114.6  48.2 -50.5 -57.3   -0.9   -5.6    8.8
   27   27   Y  H  X S+     0   0   56   -4,-2.1   4,-2.8   1,-0.2   3,-0.3   0.967 111.3  48.8 -55.3 -50.1    1.4   -8.5    9.7
   28   28   L  H  X S+     0   0   90   -4,-2.6   4,-1.9   1,-0.3   6,-0.3   0.904 108.7  53.2 -63.6 -42.9    4.3   -7.3    7.8
   29   29   N  H  <>S+     0   0   33   -4,-2.8   5,-2.1   1,-0.2  -1,-0.3   0.788 112.6  46.5 -54.5 -43.0    4.1   -3.7    9.3
   30   30   V  H ><5S+     0   0   21   -4,-1.5   3,-1.4  -3,-0.3  -2,-0.2   0.923 111.9  48.6 -70.1 -41.8    4.1   -5.4   12.8
   31   31   V  H 3<5S+     0   0   85   -4,-2.8  -2,-0.2   1,-0.3  -1,-0.2   0.821 115.5  44.5 -66.3 -38.9    7.0   -7.7   12.2
   32   32   T  T 3<5S-     0   0   91   -4,-1.9  -1,-0.3  -5,-0.3  -2,-0.2   0.487 109.4-125.4 -84.1   6.3    9.1   -4.9   10.7
   33   33   R  T < 5S+     0   0  215   -3,-1.4   2,-0.5  -4,-0.2  -3,-0.2   0.680  75.0 122.0  52.5  39.6    8.0   -2.5   13.7
   34   34   H      < +     0   0  101   -5,-2.1  -1,-0.1  -6,-0.3  -4,-0.1  -0.695  33.2 173.5-133.3  74.7    6.8    0.0   11.1
   35   35   R              0   0  180   -2,-0.5  -1,-0.1  -5,-0.1  -5,-0.1   0.257 360.0 360.0 -77.6  13.2    3.1    0.6   11.9
   36   36   Y              0   0  224   -7,-0.1  -2, 0.0   0, 0.0   0, 0.0  -0.827 360.0 360.0 -69.4 360.0    2.8    3.4    9.3
sample output file, from PDB entry 1PPT, cut before this line. ----------  */

/*-------------------------------    end of file   ------------------------*/



/* End. */
