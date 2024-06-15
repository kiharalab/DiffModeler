#ifdef __GNUG__
#  pragma implementation
#endif

#include "DsspCMBI.h"
#include "CalcAccSurf.h"
#include "AccSurf.h"
#include "Contacts.h"
#include "Vector.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


#define RWATER          1.4

static long NeighbourRes[NMAX + 1];
static long LastNeighbourRes;

static void Listentry(CAS* cas, double *xx, double *neighbour)
{
    double sqdist=Distsq(xx, neighbour);
    if(sqdist < (xx[3]+neighbour[3])*(xx[3]+neighbour[3]))
	CASAddNeigbourAtom(cas,
			   neighbour[0],neighbour[1], neighbour[2],
			   neighbour[3],sqdist);
}  /* Listentry */

static inline int InBox(double *v, double *vmin, double *vmax)
{
  if (v[0] - v[3] < vmax[0] && v[1] - v[3] < vmax[1] &&
      v[2] - v[3] < vmax[2] && v[0] + v[3] > vmin[0] &&
      v[1] + v[3] > vmin[1] && v[2] + v[3] > vmin[2])
    return 1;
  else
    return 0;
}


static void Liste(CAS cas, double *atom)
{
    long i, k, FORLIM;
    Backbone *WITH;
    long FORLIM1;

    FORLIM = LastNeighbourRes;
    for (i = 1; i <= FORLIM; i++) {
	WITH = &chain[NeighbourRes[i]];
	if (InBox(atom, WITH->boxmin, WITH->boxmax)) {
	    Listentry(cas, atom, WITH->n);
	    Listentry(cas, atom, WITH->ca);
	    Listentry(cas, atom, WITH->c);
	    Listentry(cas, atom, WITH->o);
	    if (WITH->nsideatoms > 0) {
		FORLIM1 = WITH->nsideatoms;
		for (k = 0; k < FORLIM1; k++)
		    Listentry(cas, atom, sidechain[WITH->atompointer + k]);
	    }
	}
    }
}  /* Liste */



void FindNeighbourRes(double *vmin, double *vmax)
{
  long i, FORLIM;

  LastNeighbourRes = 0;
  FORLIM = lchain;
  /* gcc (2.5.7 and 2.4.5) have problems on some machines with */
  /* -funroll-loops here.*/
  /* e.g on SGI 1tim produces wrong results!!!! */
  for (i = 1; i <= FORLIM; i++) {
      if (chain[i].aa != '!' &&
	  vmin[0] < chain[i].boxmax[0] && vmin[1] < chain[i].boxmax[1] &&
	  vmin[2] < chain[i].boxmax[2] && vmax[0] > chain[i].boxmin[0] &&
	  vmax[1] > chain[i].boxmin[1] && vmax[2] > chain[i].boxmin[2]) {
	  LastNeighbourRes++;
	  NeighbourRes[LastNeighbourRes] = i;
      }
  }
}  /*  FindNeighbourRes*/


static double Surface(CAS* cas,double *xatom)
{
    double f;
    CASResetAtoms(cas);
    CASSetCenterAtom(cas,xatom[0],xatom[1],xatom[2],xatom[3]);
    Liste(cas, xatom);
    f= CASSurface(cas);
    return f;
}  /* Surface */



void Flagaccess(int order)
{
    long i, k;
    double f;
    long FORLIM;
    Backbone *WITH;
    long FORLIM1;
    CAS *cas= CASCreate(order);
    AddToAllAtomRadii(RWATER); /* add the water radius :-) */
    CalcResidueCenters();
    FORLIM = lchain;
    for (i = 1; i <= FORLIM; i++) {
	if (chain[i].aa != '!') {
	    WITH = &chain[i];
	    FindNeighbourRes(WITH->boxmin, WITH->boxmax);
	    f = Surface(cas,WITH->n) + Surface(cas, WITH->ca) +
		Surface(cas, WITH->c) + Surface(cas, WITH->o);
	    if (WITH->nsideatoms > 0) {
		FORLIM1 = WITH->nsideatoms;
		for (k = 0; k < FORLIM1; k++)
		    f += Surface(cas, sidechain[WITH->atompointer + k]);
	    }
	    WITH->access = (long)floor(f + 0.5);
	}
    }
    AddToAllAtomRadii(-RWATER);   /* remove the water :-) */
    CASDelete(cas);
}  /* Flagaccess */

