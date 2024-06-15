#ifdef __GNUG__
#  pragma implementation
#endif

#include "Contacts.h"
#include "DsspCMBI.h"
#include "Vector.h"

/***/
/*--------------------------------------------------------------------*/

static void MinMax(double *v, double *vmin, double *vmax)
{
  long i;

  for (i = 0; i <= 2; i++) {
    if (v[i] - v[3] < vmin[i])
      vmin[i] = v[i] - v[3];
    if (v[i] + v[3] > vmax[i])
      vmax[i] = v[i] + v[3];
  }
}  /* MinMax */


void AddToAllAtomRadii(double r)
{
  long i, j, FORLIM;
  Backbone *WITH;
  long FORLIM1;

  FORLIM = lchain;
  for (i = 1; i <= FORLIM; i++) {
    if (chain[i].aa != '!') {
      WITH = &chain[i];
      WITH->h[3] += r;
      WITH->n[3] += r;
      WITH->ca[3] += r;
      WITH->c[3] += r;
      WITH->o[3] += r;
      if (WITH->nsideatoms > 0) {
	FORLIM1 = WITH->nsideatoms;
	for (j = 0; j < FORLIM1; j++)
	  sidechain[WITH->atompointer + j][3] += r;
      }
    }
  }
}  /* AddToAllAtomRadii */


void CalcResidueCenters()
{
  long i, j, FORLIM;
  Backbone *WITH;
  long FORLIM1;

  FORLIM = lchain;
  for (i = 1; i <= FORLIM; i++) {
    if (chain[i].aa != '!') {
      WITH = &chain[i];
      WITH->boxmax[0] = -1e10;
      WITH->boxmax[1] = -1e10;
      WITH->boxmax[2] = -1e10;
      WITH->boxmin[0] = 1e10;
      WITH->boxmin[1] = 1e10;
      WITH->boxmin[2] = 1e10;
      MinMax(WITH->n, WITH->boxmin, WITH->boxmax);
      MinMax(WITH->ca, WITH->boxmin, WITH->boxmax);
      MinMax(WITH->c, WITH->boxmin, WITH->boxmax);
      MinMax(WITH->o, WITH->boxmin, WITH->boxmax);
      if (WITH->nsideatoms > 0) {
	FORLIM1 = WITH->nsideatoms;
	for (j = 0; j < FORLIM1; j++)
	  MinMax(sidechain[WITH->atompointer + j], WITH->boxmin, WITH->boxmax);
      }
    }
  }
}  /* CalcResidueCenters */



