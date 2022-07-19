/*
 * filename: fv_solution.c
 * Author: Akatsuki Kimura <akkimura@nig.ac.jp>
 *         Akira Funahashi <funa@bio.keio.ac.jp>
 * Last modified: Thu, 25 Jul 2013 02:59:10 +0900
 */

#include "mtsim.h"

void FV_solution(double xx, double *f_v, double *fp_v, mtGlobal* g) {
  *f_v = (g->k_on*(exp(-1*g->F_dependency*xx)-1) + g->Vg) - ((xx - g->Buckling_backward_sum)/g->Stokes_translation);
  *fp_v = -1*g->k_on*g->F_dependency*exp(-1*g->F_dependency*xx)-1/g->Stokes_translation;
}

