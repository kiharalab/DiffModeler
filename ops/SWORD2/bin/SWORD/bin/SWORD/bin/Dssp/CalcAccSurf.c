#ifdef __GNUG__
#  pragma implementation
#endif

#include "Vector.h"
/*#include <malloc.h>*/
#include <stdlib.h>

/***************  ARRAY DIMENSIONING CONSTANTS  ***********************
  NFACE,   - NUMBER OF FACES OF POLYHEDRON. THE COORDINATES OF THE CENTRE
  ORDER      OF EACH TRIANGULAR FACE ARE STORED IN ARRAY P, THE AREA
  IS STORED IN ARRAY WP IN PROCEDURE FLAGACCESS. NFACE MUST
  BE OF THE FORM NFACE=20*(4**ORDER), ORDER=0,1,2,...
  THE ACCURACY OF THE SOLVENT ACCESSIBLE SURFACE OF EACH
  AMINOACID RESIDUE IS ONE ANGSTROM**2 FOR ORDER=2,NFACE=320.
  MAXNEIGHBOURS  - MAXIMUM NUMBER OF PROTEIN ATOMS WHICH CAN INTRUDE INTO
  SOLVENT AROUND ANY GIVEN TEST ATOM. THE COORDINATES OF
  THESE ATOMS ARE STORED IN ARRAY X, THEIR RADII IN ARRAY RX
  IN PROCEDURE SURFACE.
  --------------------------------------------------------------------*/

typedef double Vector[3];
typedef double VectorRad[4]; /* the 4th component is the radius!! */


typedef struct Polyeder {
    Vector *p;
    double *wp;
    char *accept;
    int np;
    int iter;
} Polyeder;


#define CAS CAS_


typedef struct CAS {
    Polyeder* polyeder;
    VectorRad atom;
    VectorRad *neighbours;
    int nneighbours;
    int maxneighbours;
    int wantdots;
} CAS;

#include "CalcAccSurf.h"

#define SQ(x) ((x)*(x))

#define MAXNEIGHBOURS         256

/*******************  MATHEMATICAL CONSTANTS  **************************
  YVERTEX, - ARE Y,Z-COMPONENTS OF THE FIRST ICOSAHEDRON VERTEX. THE
  ZVERTEX    X-COMPONENT IS 0.
  EPS      - NUMERICAL TOLERANCE
  --------------------------------------------------------------------*/

#define FOURPI          12.56637
#define YVERTEX         0.8506508
#define ZVERTEX         0.5257311
#define EPS             0.00001

static void Triangle(Polyeder *polyed, double *x1, double *x2, double *x3,
		     long level)
{
    long k, level1;
    double xnorm;
    Vector x4, x5, x6;

    if (level > 0) {
	level1 = level - 1;
	for (k = 0; k <= 2; k++) {
	    x4[k] = x1[k] + x2[k];
	    x5[k] = x2[k] + x3[k];
	    x6[k] = x1[k] + x3[k];
	}
	Norm(x4, &xnorm);
	Norm(x5, &xnorm);
	Norm(x6, &xnorm);
	Triangle(polyed, x1, x4, x6, level1);
	Triangle(polyed, x4, x5, x6, level1);
	Triangle(polyed, x4, x2, x5, level1);
	Triangle(polyed, x5, x3, x6, level1);
	return;
    }
    for (k = 0; k <= 2; k++)
	x6[k] = x1[k] + x2[k] + x3[k];
    Norm(x6, &xnorm);
    for (k=0; k<3 ;k++)
	polyed->p[polyed->np][k]=x6[k];
    Diff(x3, x1, x5);
    Diff(x2, x1, x4);
    Cross(x5, x4, x6);
    Norm(x6, &xnorm);
    polyed->wp[polyed->np] = xnorm / 2.0;
    polyed->np++;
}  /* Triangle */

static void PolyederReset(Polyeder *polyed)
{
    int i;
    for(i=0;i<polyed->np;i++)
	polyed->accept[i]=0;
}

static void PolyederInit(long order, Polyeder *polyed)
{  /* GENERATES ALL 12 VERTICES OF ICOSAHEDRON */
    Vector v[12];
    double a, b;
    long i, j, k, level;

    k = 0;
    a = YVERTEX;
    b = ZVERTEX;
    for (i = 0; i < 2; i++) {
	a = -a;
	for (j = 0; j < 2; j++,k++) {
	    b = -b;
	    v[k][0] = 0.0;
	    v[k][1] = a;
	    v[k][2] = b;
	    k++;
	    v[k][0] = b;
	    v[k][1] = 0.0;
	    v[k][2] = a;
	    k++;
	    v[k][0] = a;
	    v[k][1] = b;
	    v[k][2] = 0.0;
	}
    }
    polyed->np = 0;
    level = order;
    /* GET ALL 20 FACES OF ICOSAHEDRON */
    for (i = 0; i <= 9; i++) {   /* FIND INTEGRATION POINTS */
	for (j = i + 1; j <= 10; j++) {
	    if (Distance(v[i], v[j]) < 1.1) {
		for (k = j + 1; k <= 11; k++) {
		    if ((Distance(v[i], v[k]) < 1.1) & (Distance(v[j], v[k]) < 1.1))
			Triangle(polyed, v[i], v[j], v[k], level);
		}
	    }
	}
    }
    a = 0.0;
    for (i = 0; i < polyed->np; i++)
	a += polyed->wp[i];
    a = FOURPI / a;
    for (i = 0; i < polyed->np; i++)
	polyed->wp[i] *= a;
    PolyederReset(polyed);
}  /* CreatePolyeder EHCIR NAE */

static Polyeder * PolyederCreate(long order)
{  
    int i,np=20;
    Polyeder *polyeder;
    for(i=0;i<order;i++)
	np*=4;
    polyeder=(Polyeder*)malloc(sizeof(Polyeder));
    polyeder->np=np;
    polyeder->iter=np;
    polyeder->p=(Vector*)malloc(sizeof(Vector)*np);
    polyeder->wp=(double*)malloc(sizeof(double)*np);
    polyeder->accept=(char*)malloc(sizeof(char)*np);
    PolyederInit(order,polyeder);
    return polyeder;
}

static void PolyederDelete(Polyeder*polyeder)
{  
    if(polyeder) {
	free(polyeder->p);
	free(polyeder->wp);
	free(polyeder->accept);
	free(polyeder);
    }
}
/*-------------------------------------------------------------------------*/

CAS* CASCreate(int order)
{
    CAS *cas=(CAS*)malloc(sizeof(CAS));
    cas->wantdots=0;
    cas->polyeder=PolyederCreate(order);
    cas->nneighbours=0;
    cas->maxneighbours=MAXNEIGHBOURS;
    cas->neighbours=(VectorRad*)malloc(sizeof(VectorRad)*cas->maxneighbours);
    return cas;
}

void CASDelete(CAS*cas)
{
    if(cas) {
	PolyederDelete(cas->polyeder);
	free(cas->neighbours);
	free(cas);
    }
}

void CASResetAtoms(CAS*cas)
{
    cas->nneighbours=0;
}

void CASSetCenterAtom(CAS*cas, double x, double y, double z, double radius)
{
    cas->atom[0]=x;
    cas->atom[1]=y;
    cas->atom[2]=z;
    cas->atom[3]=radius;
}

void CASAddNeigbourAtom(CAS*cas, 
			double x, double y, double z,
			double radius, double squaredist)
{
    /*VectorRad atom={x,y,z,radius};*/
    VectorRad atom;
    atom[0]=x;
    atom[1]=y;
    atom[2]=z;
    atom[3]=radius;
    if (squaredist <= EPS)
	return;
    if (cas->nneighbours >= cas->maxneighbours) {
	cas->maxneighbours*=2;
	cas->neighbours=(VectorRad*)realloc(cas->neighbours,sizeof(VectorRad)*cas->maxneighbours);
    }
    Diff(atom, cas->atom, cas->neighbours[cas->nneighbours]);
    cas->neighbours[cas->nneighbours][3] = 
	(SQ(cas->atom[3]) - SQ(radius) + squaredist) / (2 * cas->atom[3]);
    if(cas->neighbours[cas->nneighbours][3] < cas->neighbours[0][3]) {
	double tmp;
	int i;
	for(i=0;i<4;i++) {
	    /*put the atom covering most dots as the first in the list... */
	    tmp=cas->neighbours[cas->nneighbours][i];
	    cas->neighbours[cas->nneighbours][i]=cas->neighbours[0][i];
	    cas->neighbours[0][i]=tmp;
	}
    }
    cas->nneighbours++;
}

double CASSurface(CAS*cas)
{
    int i,k;
    double f;
    int lastk;
    Vector *p=cas->polyeder->p;
    int np= cas->polyeder->np;
    VectorRad *x=cas->neighbours;
    int nx=cas->nneighbours;
    if(cas->wantdots)
	PolyederReset(cas->polyeder);
    f = 0.0;
    lastk = 0;
    for (i = 0; i < np; i++) {
	/* here is the major speedup of the algorythm */
	/* Try first the atom, that kicked out the last point */
	if (Dot(p[i], x[lastk]) <= x[lastk][3] ) {
	    /* only if the last atom didn't cover the point */
	    /* search the neighbours */
	    for (k = 0; k < nx; k++)
		if(Dot(p[i], x[k]) > x[k][3])
		    break;
	    if(k<nx) {
		lastk = k; /* atom k coveres point i */
	    } else {
		f += cas->polyeder->wp[i]; /* point i on the surface */
		cas->polyeder->accept[i]=1;
	    }
	}
    }
    /* scale it with the radius of the current atom */
    return SQ(cas->atom[3])*f;
}

void CASNResetNextPoint(CAS*cas)
{
    if(cas->wantdots)
	cas->polyeder->iter=0;
    else
	cas->polyeder->iter=cas->polyeder->np; /* go to the end...*/
}

int CASGetNextPoint(CAS*cas, double *x, double *y, double *z, double *surface)
{
    int i;
    int n=cas->polyeder->np;
    for(i=cas->polyeder->iter;i<n;i++) {
	if(cas->polyeder->accept[i]) {
	    *x=cas->polyeder->p[i][0]*cas->atom[3]+cas->atom[0];
	    *y=cas->polyeder->p[i][1]*cas->atom[3]+cas->atom[1];
	    *z=cas->polyeder->p[i][2]*cas->atom[3]+cas->atom[2];
	    if(surface)
		*surface=SQ(cas->atom[3])*cas->polyeder->wp[i];
	    cas->polyeder->iter=i+1;
	    return 1;
	}
    }
    return 0;
}
