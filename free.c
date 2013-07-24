/*
 * Author: Akatsuki Kimura <akkimura@nig.ac.jp>
 *         Akira Funahashi <funa@bio.keio.ac.jp>
 * Last modified: Thu, 25 Jul 2013 02:37:35 +0900
 */
#include "mtsim.h"

void free_return(double *fvec, double **fjac, int n, double *p, int *indx) {
  free_dmatrix(fjac,1,n,1,n);
  free_dvector(fvec,1,n);
  free_dvector(p,1,n);
  free_ivector(indx,1,n);
  return;
}
