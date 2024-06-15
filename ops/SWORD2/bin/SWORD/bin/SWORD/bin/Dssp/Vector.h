#ifndef _Vector_h_
#define _Vector_h_

#ifdef __GNUG__
#  pragma interface
#endif

#ifndef GCC
  #define inline
#endif

/* a vector has 3 components: x,y,z!!! */
double Distance(double *u, double *v);
double VLength(double *u);

static inline double Distsq(double *u, double *v)
{
    return  (u[0] - v[0]) * (u[0] -v[0])
	  + (u[1] - v[1]) * (u[1] -v[1])
	  + (u[2] - v[2]) * (u[2] -v[2]);
} 


static inline void Diff(double *x, double *y, double *z)
{
    z[0] = x[0] - y[0];
    z[1] = x[1] - y[1];
    z[2] = x[2] - y[2];
}

static inline double Dot(double *x, double *y)
{
    return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

static inline void Cross(double *x, double *y, double *z)
{
    z[0] = x[1] * y[2] - y[1] * x[2];
    z[1] = x[2] * y[0] - y[2] * x[0];
    z[2] = x[0] * y[1] - y[0] * x[1];
}


void Norm(double *x, double *xnorm);


#endif
