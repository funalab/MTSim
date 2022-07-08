/*
 * LU decomposition (ref) Press et al., "Numerical Recipes in C"
 * Author: Akatsuki Kimura <akkimura@nig.ac.jp>
 *         Akira Funahashi <funa@bio.keio.ac.jp>
 * Last modified: Thu, 25 Jul 2013 02:39:29 +0900
 */
#include "mtsim.h"

void lubksb(double **a, int n, int *indx, double b[])
{
  int i,ii=0,ip,j;
  double sum;

  for (i=1;i<=n;i++){
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  for (i=n;i>=1;i--) {
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}
