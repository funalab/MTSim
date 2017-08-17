/*
 * Author: Akatsuki Kimura <akkimura@nig.ac.jp>
 *         Akira Funahashi <funa@bio.keio.ac.jp>
 * Last modified: Thu, 25 Jul 2013 03:01:46 +0900
 */

#ifndef __mtsim__
#define __mtsim__
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mtgraphics.h"
#include "nrutil.h"
#include "mtglobal.h"

#define false 0
#define true 1

/* debug print */
#ifdef DEBUG_PRINT
#define DEBUG_PRINT_FLAG 1
#else
#define DEBUG_PRINT_FLAG 0
#endif

#define TRACE(x) do { if (DEBUG_PRINT_FLAG) dbg_printf x; } while (0)

/* boolean */
#define true 1
#define false 0
#ifdef boolean
#undef boolean
#endif
#define boolean int

/* color */
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
/* #define ST (30000)           /\* total steps (25min/0.05sec=30000)*\/ */
#define ST (24000)           /* total steps (20min/0.05sec=30000)*/
#define laserST (30000)      /* when each centrosome moves independently */
#define anaST (24000)      /* when mitotic spindle begins to elongate */
#define metaST (18000)       /* when repression of cortical pulling force at "LET-99 band" is inactivated */
#define dT (0.05)            /* time step [sec] default = 0.05 */
/* #define LL (5*pow(10,-6))    /\* radius of pronucleus [m] *\/ */
/* #define Rad (25*pow(10,-6))  /\* long axis of egg [m] *\/ */
/* #define RadS (15*pow(10,-6)) /\* short axes of egg [m] *\/ */
#define Cir_Rad (19.487570*pow(10,-6)) 

/* declaration of functions */
/* callback.c */
void function_FV3D(double *x, int n, double *fvec, double **fjac, mtGlobal *g);
void function_MotorFV (double *x, int n, double *fvec, double **fjac, mtGlobal *g);
void function_laserMotorFV (double *x, int n, double *fvec, double **fjac, mtGlobal *g);
void function_anaMotorFV (double *x, int n, double *fvec, double **fjac, mtGlobal *g);

/* fv_solution.c */
void FV_solution(double xx, double *f_v, double *fp_v, mtGlobal *g);

/* mnewt.c */
boolean mnewt(int ntrial, double x[], int n, double tolx, double tolf,
    int step_counter, FILE* f_out8, void (*usrfun)(double *, int , double*, double**, mtGlobal*), mtGlobal* g);

/* solve_1D_double.c */
double Stokes_function(double ff, double gr, double fb);
double FV_function(double ff, double vg, double ko, double fd);
double rtsafe(void (*funcd)(double, double *, double *, mtGlobal*),double x1, double x2, double xacc, mtGlobal* g);
double rtsafe_mod(void (*funcd)(double, double *, double *, mtGlobal*),double x1, double x2, double xacc, mtGlobal* g);

/* ludcmp_double.c */
void ludcmp(double **a, int n, int *indx, double *d);

/* lubksb_double.c */
void lubksb(double **a, int n, int *indx, double b[]);

/* ran1_double.c */
double ran1(long *idum);

/* main.c */
void usage(char* myname);

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

/* store_graphs.c */
void store_graphs(void);

/* save_logs.c */
/* void save_logs(int i, int p, int N, FILE* f_out1, FILE* f_out2, FILE* f_out3, FILE* f_out4, FILE* f_out5, FILE* f_out6, FILE* f_out7, FILE* f_out8, FILE* f_out9, FILE* f_out10, FILE* data_for_3D, double Nuc[3], double PVecCen[2][3], double MT[][3]); */
void save_logs(int i, int p, int N, FILE* f_out1, FILE* f_out2, FILE* f_out3, FILE* f_out4, FILE* f_out5, FILE* data_for_3D, double PVecCen[2][3], double MT[][3]);

/* free.c */
void free_return(double *fvec, double **fjac, int n, double *p, int *indx);

/* draw_graphs.c */
/* void draw_graphs(int i, mtGraphics *mtg, mtGlobal *g, double PVecCen[2][3], double MT[][3], double Nuc[3], double DistanceFromPP, double min, double max, double sum, int currentN, double OldVel, double NewVel); */
void draw_graphs(int i, mtGraphics *mtg, mtGlobal *g, double PVecCen[2][3], double MT[][3], double Nuc[3], double Rad, double RadS, double MetaSpindle_L);

/* display_setting.c */
void display_setting(mtGraphics *mtg);

/* color_setting.c */
void color_setting(mtGraphics *mtg, int p);

/* dbg_printf.c */
void dbg_printf(const char *fmt, ...);

#endif
