/*
 * Last modified: Sun, 30 Jun 2013 02:53:44 +0900
 */

#ifndef __mtsim__
#define __mtsim__
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mtgraphics.h"
#include "nrutil.h"

#define false 0
#define true 1

#define DEBUG 0

#define RED     "\x1b[31m"
#define GREEN   "\x1b[32m"
#define YELLOW  "\x1b[33m"
#define BLUE    "\x1b[34m"
#define MAGENTA "\x1b[35m"
#define CYAN    "\x1b[36m"
#define DEFAULT "\x1b[39m"
#define UNDER_LINE   "\x1b[4m"
#define BOLD         "\x1b[1m"
#define FONT_DEFAULT "\x1b[0m"

/* parameters used also in other files  */
#define PI (3.141592653590)  /* Pi */
#define WIN_WIDTH (256)      /* graphical settings */
#define WIN_HEIGHT (256)     /*  graphical settings */
#define SC3 (0.01)           /* graphical settings */
#define UNITTIME (10)        /* for calculation of velocities */
#define ST (24000)           /* total steps (20min/0.05sec=24000)*/
#define laserST (24000)      /* when each centrosome moves independently */
#define metaST (18000)       /* when repression of cortical pulling force at "LET-99 band" is inactivated */
#define dT (0.05)            /* time step [sec] default = 0.05 */
#define LL (5*pow(10,-6))    /* radius of pronucleus [m] */
#define Rad (25*pow(10,-6))  /* long axis of egg [m] */
#define RadS (15*pow(10,-6)) /* short axes of egg [m] */

/* declaration of functions */
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);
double ran1(long *idum);
double Stokes_function(double ff, double gr, double fb);
double FV_function(double ff, double vg, double ko, double fd);
double rtsafe(void (*funcd)(double, double *, double *),double x1, double x2, double xacc);
double rtsafe_mod(void (*funcd)(double, double *, double *),double x1, double x2, double xacc);

/* math_func.c */
void UnitVector(double vector1[3], double vector2[3]);
void AddVector(double vector1[3], double vector2[3], double result[3]);
void SubVector(double vector1[3], double vector2[3], double result[3]);
void OutProdVector(double vector1[3], double vector2[3], double result[3]);
void QuadEqu2(double a[3], double b[2]);
void ProductJacVec(double aa_out[3], double bb_jac[3][3], double cc_in[3]);
void VectorRotation(double out[3], double in[3], double rotAx[3], double rotStd[3], double degree);
void MakeRotationMatrix(double RotationMatrix[3][3], double RotationVector[3], double dt);
double Length(double a[3]);
double InnProdVector(double vector1[3], double vector2[3]);
double Poisson(double np, int motornumber);
void FV_solution(double xx, double *f_v, double *fp_v);
void function_FV3D(double *x, int n, double *fvec, double **fjac);
void function_MotorFV (double *x, int n, double *fvec, double **fjac);
void function_laserMotorFV (double *x, int n, double *fvec, double **fjac);
void mnewt(int ntrial, double x[], int n, double tolx, double tolf);

/* store_graphs.c */
void store_graphs(void);

/* save_logs.c */
void save_logs(int i, int p, int N, FILE* f_out1, FILE* f_out2, FILE* f_out3, FILE* f_out4, FILE* f_out10, FILE* data_for_3D, double Nuc[3], double PVecCen[2][3], double MT[][3]);

/* free.c */
void free_return(double *fvec, double **fjac, int n, double *p, int *indx);

/* draw_graphs.c */
void draw_graphs(int i, mtGraphics *mtg, double PVecCen[2][3], double MT[][3], double Nuc[3], double DistanceFromPP, double min, double max, double sum, int currentN, double OldVel, double NewVel);

/* display_setting.c */
void display_setting(mtGraphics *mtg);

/* color_setting.c */
void color_setting(mtGraphics *mtg, int p);

#endif
