#include "Vector.h"
#include <math.h>

double Distance(double *u, double *v)
{
    return sqrt(Distsq(u,v));
}

void Norm(double *x, double *xnorm)
{
  /* RETURNS INPUT VECTOR X NORMALIZED TO UNIT LENGTH.
     XNORM IS THE ORIGINAL LENGTH OF X.                         */
  double TEMP, TEMP1, TEMP2;

  TEMP = x[0];
  TEMP1 = x[1];
  TEMP2 = x[2];
  *xnorm = TEMP * TEMP + TEMP1 * TEMP1 + TEMP2 * TEMP2;
  if (*xnorm <= 0.0)
    return;
  *xnorm = sqrt(*xnorm);
  x[0] /= *xnorm;
  x[1] /= *xnorm;
  x[2] /= *xnorm;
}

double VLength(double *u)
{
    return sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
}
