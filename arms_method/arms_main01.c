/* example calling program for adaptive rejection sampling - normal density */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "arms.h"

/* ********************************************************************* */

struct norm_parm {
  double mean,sd;
};

/* ********************************************************************* */

double norm(double x, void *norm_data)
/* normal log density */
{
   struct norm_parm *d;
   double y;

   /* cast voided pointer into pointer to struct norm_parm */
   d = norm_data;
   y = -0.5*pow(((x - d->mean)/d->sd),2.0);
   return y;
};

/* ********************************************************************* */

struct normix_parm {
  double mean1,sd1,mean2,sd2,p1;
};

/* ********************************************************************* */

double normix(double x, void *normix_data)
/* normal mixture log density */
{
   struct normix_parm *d;
   double y;

   /* cast voided pointer into pointer to struct normix_parm */
   d = normix_data;
   y = log((d->p1/d->sd1)*exp(-0.5*pow(((x - d->mean1)/d->sd1),2.0))
     + ((1-d->p1)/d->sd2)*exp(-0.5*pow(((x - d->mean2)/d->sd2),2.0)));
   return y;
};

/* ********************************************************************* */

int main(void)
{
  FILE *f;
  int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 4 ;
  int neval,i;
  double xinit[10]={0.0,3.0,17.0,20.0}, xl = -100.0, xr = 100.0;
  double xsamp[100], xcent[10], qcent[10] = {5., 30., 70., 95.};
  unsigned seed = 44;
  double convex = 1.;
  int dometrop = 0;
  double xprev = 0.0;

  /* set up structures for each density function */

  struct norm_parm norm_data;
  struct normix_parm normix_data;

  /* initialise data needed by normal density function */
  norm_data.mean = 10.0;
  norm_data.sd = 5.;

  /* initialise data needed by normal mixture density function */
  normix_data.p1 = .3;
  normix_data.mean1 = 5.;
  normix_data.sd1 = 1.;
  normix_data.mean2 = 10.;
  normix_data.sd2 = 2.5;

  /* initialise random number generator */
  srand(seed);

  /* open a file */
  f = fopen("arms.out01","w+");

  for(i=0;i<100000;i++){
    err = arms(xinit,ninit,&xl,&xr,norm,&norm_data,&convex,
              npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
    if(err>0){
      fprintf(f,"error code = %d\n",err);
      exit(1);
    }
    fprintf(f,"%d  %11.5f %d\n",i,xsamp[0],neval);
  }
  return 0;
}

