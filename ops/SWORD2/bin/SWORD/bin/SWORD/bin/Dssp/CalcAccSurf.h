#ifndef _CalcAccSurf_h_
#define _CalcAccSurf_h_

#ifdef __GNUG__
#  pragma interface
#endif

#ifndef CAS
#  define CAS void*
#endif

CAS*   CASCreate(int order);
void   CASDelete(CAS*cas);

void   CASResetAtoms(CAS*cas);
void   CASSetCenterAtom(CAS*cas,
			double x, 
			double y, 
			double z, 
			double radius);
void   CASAddNeigbourAtom(CAS*cas,
			  double x,
			  double y,
			  double z,
			  double radius,
			  double squaredist);
double CASSurface(CAS*cas);

void CASNResetNextPoint(CAS*cas);
int CASGetNextPoint(CAS*cas, double *x, double *y, double *z, double *surface);

#endif
