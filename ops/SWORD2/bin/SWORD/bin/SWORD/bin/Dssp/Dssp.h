#ifndef _Dssp_h_
#define _Dssp_h_

#ifdef __GNUG__
#  pragma interface
#endif

/***************  ARRAY DIMENSIONING CONSTANTS  ***********************
 NMAX     - MAXIMUM NUMBER OF AMINOACID RESIDUES IN ARRAY CHAIN
 MAXATOM  - MAXIMUM NUMBER OF SIDECHAIN ATOMS IN ARRAY SIDECHAIN
  --------------------------------------------------------------------*/
#define NMAX            20000
#define MAXATOM         133000L


typedef double Vector[4]; /* the 4th component is the radius!! */

typedef enum {
  symbol, turn3, turn4, turn5, bend, chirality, beta1, beta2
} structure;

typedef char Char4[4];
typedef char Char6[6];

typedef struct HydrogenBond {
  long residue, energy;
} HydrogenBond;

typedef HydrogenBond Bonds[2];

typedef struct Backbone {
  Char6 aaident;
  char sheetlabel, aa;
  Char4 threelettercode;
  char ss[(long)beta2 - (long)symbol + 1];
  long partner[(long)beta2 - (long)beta1 + 1];
  long access;
  double alpha, kappa;
  Bonds acceptor, donor;
  Vector boxmin;   /* minimun of a box around the residue */
  Vector boxmax;   /* maximum of a box around the residue */
  Vector h, n, ca, c, o;
  long atompointer, nsideatoms;
} Backbone;

extern Backbone chain[NMAX + 1];
extern long lchain;
extern Vector sidechain[MAXATOM];

#endif
