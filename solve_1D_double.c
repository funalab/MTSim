/*
 * Last modified: Wed, 17 Apr 2013 04:09:55 +0900
 */
#include <stdio.h>
#include <math.h>
#define MAXIT 100
#define PI (3.141592653590) /* Pi */

void FV_solution(double xx, double *f_v, double *fp_v);

double Stokes_function(double ff, double gr, double fb)
{
  double result;
  result = (ff-fb)/gr;
  return result;
}

double FV_function(double ff, double vg, double ko, double fd)
{
  double result;
  result = ko*(exp(-1*fd*ff)-1)+vg;
  return result;
}

double rtsafe(void (*funcd)(double,double *, double *), double x1, double x2, double xacc)
{
  void nrerror(char error_text[]);
  int j;
  double df,dx,dxold,f,fh,fl;
  double temp,xh,xl,rts;

  (*funcd)(x1,&fl,&df);
  (*funcd)(x2,&fh,&df);
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
    printf("Root must be bracketed in rtsafe\n");
  if (fl == 0.0) return x1;
  if (fh == 0.0) return x2;
  if (fl <0.0) {
    xl=x1;
    xh=x2;
  } else {
    xh=x1;
    xl=x2;
  }
  rts=0.5*(x1+x2);
  dxold=fabs(x2-x1);
  dx=dxold;
  (*funcd)(rts,&f,&df);
  for (j=1;j<=MAXIT;j++) {
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)|| (fabs(2.0*f) > fabs(dxold*df))) {
      dxold=dx;
      dx=0.5*(xh-xl);
      rts=xl+dx;
      if (xl == rts) return rts;
    } else {
      dxold=dx;
      dx=f/df;
      temp=rts;
      rts -= dx;
      if (temp == rts) return rts;
    }
    if (fabs(dx) < xacc) return rts;
    (*funcd)(rts,&f,&df);
    if (f<0.0)
      xl=rts;
    else
      xh=rts;
  }
  printf("Maximum number of iteration exceeded in rtsafe");
  return 0.0;
}

double rtsafe_mod(void (*funcd)(double,double *, double *), double x1, double x2, double xacc)
{
  void nrerror(char error_text[]);
  int j;
  double df,dx,dxold,f,fh,fl;
  double temp,xh,xl,rts;

  (*funcd)(x1,&fl,&df);
  (*funcd)(x2,&fh,&df);
  // modification starts
  if (fl > 0.0 && fh > 0.0){
    //   printf("dv is slower than Vbuckle\n");
    return x2;
  }
  if (fl < 0.0 && fh < 0.0){ 
    printf("dv is faster than Vg\n");
    return x1;
  }
// modification ends
  if (fl == 0.0) return x1;
  if (fh == 0.0) return x2;
  //     printf("Root is bracketed in rtsafe_mod\n");

  if (fl <0.0) {
    xl=x1;
    xh=x2;
  } else {
    xh=x1;
    xl=x2;
  }
  rts=0.5*(x1+x2);
  dxold=fabs(x2-x1);
  dx=dxold;
  (*funcd)(rts,&f,&df);
  for (j=1;j<=MAXIT;j++) {
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)|| (fabs(2.0*f) > fabs(dxold*df))) {
      dxold=dx;
      dx=0.5*(xh-xl);
      rts=xl+dx;
      if (xl == rts) return rts;
    } else {
      dxold=dx;
      dx=f/df;
      temp=rts;
      rts -= dx;
      if (temp == rts) return rts;
    }
    if (fabs(dx) < xacc) return rts;
    (*funcd)(rts,&f,&df);
    if (f<0.0)
      xl=rts;
    else
      xh=rts;
  }
  printf("Maximum number of iteration exceeded in rtsafe");
  return 0.0;
}
