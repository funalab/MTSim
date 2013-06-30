/*
 * LU decomposition (ref) Press et al., "Numerical Recipes in C"
 * Last modified: Mon, 01 Jul 2013 02:42:43 +0900
 */
#include "mtsim.h"
#define TINY 1.0e-20;

void ludcmp(double **a, int n, int *indx, double *d)
{
  int i,imax,j,k;
  double big,dum,sum,temp;
  double *vv;

  //    TRACE(("enter ludcmp\n"));
  vv=dvector(1,n);
  *d=1.0;
  /*
	for(iii=1;iii<=n;iii++){
	TRACE(("%d (",indx[iii]));
	for (jjj=1;jjj<=n;jjj++){
	TRACE(("%4.3f ",a[iii][jjj]*1000));}
	TRACE((")\n"));
    
	}
	TRACE(("\n"));
  */
 for (i=1;i<=n;i++){
    big=0.0;
    for (j=1; j<=n;j++){
      //      TRACE(("ludcmp: a[%d][%d]=%f\n",i,j,a[i][j]));
      //      TRACE(("ludcmp: fabs(a[%d][%d])=%f\n",i,j,fabs(a[i][j])));
      if ((temp=fabs(a[i][j])) > big) {big=temp;}}
    if (big == 0.0) {TRACE(("Singular matrix in routine ludcmp"));}
    vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sum=a[i][j];
      for (k=1;k<i;k++) {sum -= a[i][k]*a[k][j];}
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<=n;i++) {
      sum=a[i][j];
      for (k=1;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ((dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++){
	//		TRACE(("%d %d\n",imax,k));
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;   /******* missing line ****/
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY; /* original: a[j][j]=TINY */
    if (j != n) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j]*=dum;
    }
         
    
  }
  free_dvector(vv,1,n);
  //    TRACE(("end ludcmp\n"));
}
