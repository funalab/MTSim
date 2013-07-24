#include <stdio.h>
/*
 * filename: mnewt.c
 * Author: Akatsuki Kimura <akkimura@nig.ac.jp>
 *         Akira Funahashi <funa@bio.keio.ac.jp>
 * Last modified: Thu, 25 Jul 2013 01:43:28 +0900
 */

#include "mtsim.h"

// Newton-Raphson Method for Nonlinear Systems of Equations (ref) Press et al. "Numerical Recipes in C"*/
boolean mnewt(int ntrial, double x[], int n, double tolx, double tolf,
    int step_counter, FILE* f_out8, void (*usrfun)(double *, int , double*, double**) ) {
  int k,i,*indx;
  double errx,errf,d,*fvec,**fjac,*p;

  indx = ivector(1,n);
  p = dvector(1,n);
  fvec = dvector(1,n);
  fjac = dmatrix(1,n,1,n);
  for (k=1;k<=ntrial;k++) {
    usrfun(x,n,fvec,fjac); /*** specific function ****/
    errf=0.0;
    for (i=1;i<=n;i++) {
      errf += fabs(fvec[i]);
    }
    if (errf <=tolf) {
      if (step_counter%100 == 0){
        fprintf(f_out8,"%d tolf errx=%5.4f\n",k,errx*(1.0e+6));
      }
      free_return(fvec, fjac, n, p, indx);
      return true;
    }
    for (i=1; i<=n; i++) {
      p[i] = -fvec[i];
    }
    ludcmp(fjac,n,indx,&d);
    lubksb(fjac,n,indx,p);
    errx=0.0;
    for (i=1;i<=n;i++){
      errx += fabs(p[i]);
      x[i] += p[i];
    }
    if (errx <=tolx) {
      if (step_counter%100 == 0){
        fprintf(f_out8,"%d tolx errf=%5.4f\n",k,errf*(1.0e+12));
      }
      free_return(fvec, fjac, n, p, indx);
      return true;
    }
  }
  printf("mnewt: did not converge at %d",k);
  free_return(fvec, fjac, n, p, indx);
  return false;
}

