/*
 * Author: Akatsuki Kimura <akkimura@nig.ac.jp>
 *         Akira Funahashi <funa@bio.keio.ac.jp>
 * Last modified: Thu, 25 Jul 2013 02:06:24 +0900
 */

#ifndef __mtglobal__
#define __mtglobal__

/* parameters used in main(), function_* and draw_graphs() */
typedef struct _mtGlobal {
  int N; /* the (maximum) number of MTs per two centrosomes */
  int NN;
  double *L; /* MTs length */
  double Buckling_forward_sum;
  double Buckling_backward_sum;
  double Stokes_rad; /* stokes radius of cytosol [m] */
  double Visco; /* viscosity of cytosol [kg/m sec] */
  double Stokes_rotation;
  double Stokes_translation;
  double **DVecNucCen;
  double Vg;
  double k_on;
  double F_dependency;
  double **u; /* MTs unit vector */
  unsigned char *pushing_phase;
  unsigned char *pulling_phase;
  unsigned char *phase;
  double Fbuckle[3];
  double Fbackward[6];
  double fjac_pull[6][6];
  double BucklingConst;
  double *NumberOfMotor;
  unsigned int phase_transition_count;
  double MotorDensity;
  double MotorMaxVel;
  double MotorStallF;
  int mt_start;
  int mt_end;
} mtGlobal;

#endif
