/*
 * filename: main.c
 * previous version: PPNM_SIM050913_1.cc
 * change from previous version: Refactoring
 * This program simulates centrosome positioning in one-cell embryo
 * Unit meter, kilo-gram, sec
 * to compile: make (see Makefile for detail)
 * Last modified: Thu, 25 Jul 2013 01:48:41 +0900
 */

#include "mtsim.h"
#if defined(_MSC_VER) || defined(__STRICT_ANSI__)
#include "my_getopt.h"
#else
#include <unistd.h>
#endif

// parameters used in functions in this file
static double g_Buckling_forward_sum;
static double g_Buckling_backward_sum;
static double g_Stokes_rad;
static double g_Visco;
static double g_Stokes_rotation;
static double g_Stokes_translation;
static double **g_DVecNucCen;
static double g_Vg;
static double g_k_on;
static double g_F_dependency;
int g_N; /* the (maximum) number of MTs per two centrosomes */
int g_NN;
static double **g_u;
static unsigned char *g_pushing_phase;
static unsigned char *g_pulling_phase;
static unsigned char *g_phase;
static double g_Fbuckle[3];
static double g_Fbackward[6];
static double g_fjac_pull[6][6];
static double g_BucklingConst;
double *g_L;
static double *g_NumberOfMotor;
static unsigned int g_phase_transition_count;
static double g_MotorDensity;
static double g_MotorMaxVel;
static double g_MotorStallF;
static int g_mt_start;
static int g_mt_end;

void usage(char* myname) {
  printf("Usage : %s [option]\n", myname);
  printf(" -h   : Show this message\n");
  printf(" -c   : Check output with fixed random seed (0.0)\n");
  printf(" -v   : Verbose output\n");
  printf(" -t # : specify simulation time (ex. -t 100 )\n");
  printf(" -s # : specify simulation step (ex. -s 100 )\n");
  printf(" -d # : specify simulation delta (ex. -d 0.01 [default:1/4096])\n");
  printf(" -m # : specify model (ex. -m 5 )\n");
  printf("        0: Pushing\n");
  printf("        1: Pulling\n");
  printf("        2: Cortical-Anchoring\n");
  printf("        3: Pulling Supplemental 1\n");
  printf("        4: Pushing Supplemental 2\n");
  printf("        5: Length-Dependent + Cortical-Anchoring\n");
}

void FV_solution(double xx, double *f_v, double *fp_v) {
  *f_v = (g_k_on*(exp(-1*g_F_dependency*xx)-1) + g_Vg) - ((xx - g_Buckling_backward_sum)/g_Stokes_translation);
  *fp_v = -1*g_k_on*g_F_dependency*exp(-1*g_F_dependency*xx)-1/g_Stokes_translation;
}

void function_FV3D(double *x, int n, double *fvec, double **fjac) {
  /* x[1],x[2],x[3] are velocity of Nuc for x,y,z axes, respectively */
  int i,j;
  for (i=1;i<=n;i++) {
    fvec[i] = 0.0;
    for (j=1;j<=n;j++){
      fjac[i][j] = 0.0;
    }
  }
  int mt;
  double A,B,C,D,E;
  for (mt=0; mt<g_N; mt++) {
    if (g_pushing_phase[mt]==2) {
      C = g_u[mt][0]*x[1]+g_u[mt][1]*x[2]+g_u[mt][2]*x[3];
      D = C + g_Vg;
      E = D/g_k_on;
      A = 1-E;
      if (A<=0.0){
        g_pushing_phase[mt] = 1;
        g_phase_transition_count++;
        TRACE(("A<0 at function_FV3D\n"));
        for (i=1; i<=n; i++) {
          fvec[i] -= g_u[mt][i-1]*g_F_dependency*g_BucklingConst/(g_L[mt]*g_L[mt]);
        }
      } else {
        B=log(A);
        for (i=1; i<=n; i++) {
          fvec[i] += g_u[mt][i-1]*B;
          for (j=1; j<=n; j++) {
            fjac[i][j] += (g_u[mt][i-1]*g_u[mt][j-1])/A;
          }
        }
      }
    }
  }
  for (i=1;i<=n;i++){
    fvec[i] = fvec[i]/g_F_dependency + g_Fbuckle[i-1] - g_Stokes_translation*x[i];   /* Fbuckle[0,1,2] */
    for (j=1;j<=n;j++){
      fjac[i][j] = fjac[i][j]/(-1.0*g_F_dependency*g_k_on);
      if (j==i) {
        fjac[i][j] = fjac[i][j] - g_Stokes_translation;}
    }
  }
}

void function_MotorFV(double *x, int n, double *fvec, double **fjac) {
  /* x[1],x[2],x[3] are velocity of Nuc for x,y,z axes, respectively */ 
  int i,j;
  for (i=1;i<=n;i++) {
    fvec[i] = g_Fbackward[i-1];
    for (j=1;j<=n;j++){
      fjac[i][j] = g_fjac_pull[i-1][j-1];
    }
  }
  //  for (j=1; j<=6; j++) {
  //    TRACE(("fjac[%d] %lf %lf %lf %lf %lf %lf\n",j,fjac[j][1],fjac[j][2],fjac[j][3],fjac[j][4],fjac[j][5],fjac[j][6]));
  //  }

  int mt;
  int qq;
  double dv, PullingF;
  for (mt=0; mt<g_N; mt++) {
    if (g_pulling_phase[mt]==1) {
      if (mt<g_NN) {qq=0;} else {qq=1;}
      dv = 0.0;
      for (j=0; j<3; j++) {
        dv += (x[(j+2)%3+4]*g_DVecNucCen[qq][(j+1)%3]-x[(j+1)%3+4]*g_DVecNucCen[qq][(j+2)%3]+x[j+1])*g_u[mt][j];
      }
      PullingF = g_MotorStallF*(1-dv/g_MotorMaxVel);
      for (i=1; i<=n; i++) {
        fvec[i] += g_NumberOfMotor[mt]*g_u[mt][i-1]*PullingF;
      }
    }
  }
  for (i=1;i<=(n/2);i++){
    fvec[i] -= g_Stokes_translation*x[i];   /* Fbackward[0,1,2] */
  }
  for (i=(n/2)+1;i<=n;i++){
    fvec[i] -= g_Stokes_rotation*x[i];
  }
}

/////////////////// the below function is identical to function_motorFV and needs to be modified ////////////////////
void function_laserMotorFV (double *x, int n, double *fvec, double **fjac) {
  /* x[1],x[2],x[3] are velocity of Nuc for x,y,z axes, respectively */
  int i,j;
  for (i=1;i<=n;i++) {
    fvec[i]=0.0;
    for (j=1;j<=n;j++){
      fjac[i][j] = g_fjac_pull[i-1][j-1];
    }
  }
  int mt;
  int qq;
  double dv, PullingF;
  /* double dv, PullingF, dFdv; */
  /* double gradient = g_Stokes_translation; */
  for (mt=g_mt_start; mt<g_mt_end; mt++) {
    if (g_pulling_phase[mt]==1) {
      dv = g_u[mt][0]*x[1]+g_u[mt][1]*x[2]+g_u[mt][2]*x[3];
      if (mt<g_NN) {qq=0;} else {qq=1;}
      dv=0;
      for (j=0; j<3; j++) {
        dv += (x[(j+2)%3+4]*g_DVecNucCen[qq][(j+1)%3]-x[(j+1)%3+4]*g_DVecNucCen[qq][(j+2)%3]+x[j+1])*g_u[mt][j];
      }
      PullingF = g_MotorStallF*(1-dv/g_MotorMaxVel);
      for (i=1; i<=n; i++) {
        fvec[i] += g_NumberOfMotor[mt]*g_u[mt][i-1]*PullingF;
      }
    }
  }
  for (i=1;i<=(n/2);i++){
    fvec[i] = fvec[i]+g_Fbackward[i-1] - g_Stokes_translation*x[i];   /* Fbackward[0,1,2] */
  }
  for (i=(n/2)+1;i<=n;i++){
    fvec[i] = fvec[i] + g_Stokes_rotation*x[i];
  }
}

int main(int argc, char* argv[]) {
  mtGraphics mtg;
  display_setting(&mtg);

  /*  Variables for getopt() */
  int ch;
  extern char *optarg;
  extern int optind;
  char* myname;
  boolean is_check = false;
  boolean is_verbose = false;

  //random number (ref) Press et al., "Numerical Recipes in C"
  double random;
  long tt;
  long *idum;
  idum = (long *) calloc(1, sizeof(long));
  // for DEBUG funa
  //tt = (unsigned)time(NULL);
  tt = 0.0;
  printf("tt=%ld\n",tt);
  *idum = tt*(-1);

  //file handling
  FILE *f_out1 = fopen ("out1.dat","w");
  FILE *f_out2 = fopen ("out2.dat","w");
  FILE *f_out3 = fopen ("out3.dat","w");
  FILE *f_out4 = fopen ("out4.dat","w");
  FILE *f_out10 = fopen("out5.dat","w");
  FILE *f_out5 = fopen ("out_parameter.dat","w");
  FILE *f_out6 = fopen ("out_MTnumber.dat","w");
  FILE *f_out7 = fopen ("out_FVcheck.dat","w");
  FILE *f_out8 = fopen ("out_3Dcheck.dat","w");
  FILE *data_for_3D = fopen ("out_for_3D.dat","w");

  ////////////////////////////////////////////
  // DECLEARATION of Constants and Variables//
  ////////////////////////////////////////////

  /************CONSTANT NUMBER*****************/
  //PARAMETERS
  unsigned char strain=8; // 0:WT, 1:par-2 w/o LET-99, 2:par-2 with LET-99, 3:par-3 w/o LET-99, 4:par-3 with LET-99, 5: let-99, 6:ric-8, 7:let-99;ric-8, 8:gpr-1/2, 9:PosteriorCortexPullingOnly
  // dynamic instability of MT
  g_Vg = 0.328e-6; /* g_Vg: growth velocity of MT [m/sec], standard:0.12 micron/sec */
  double Vs = 0.537e-6; /* Vs: shrinkage velocity of MT [m/sec], standard:0.288 micron/sec*/
  double CatFreq = 0.046; /* fcat: catastrophe frequency of MT [/sec], standard:0.014 */
  double ResFreq = 0.133; /* fres: rescue frequency of MT [/sec], standard:0.044 */
  // drag force of pronucleus
  g_Visco = 1.0; /* viscosity of cytosol [kg/m sec], standard 0.1->1.0 */
  g_Stokes_rad = 10.0e-6;
  // pushing forces
  double EI = 10.0e-24; /* rigidity of MT [Nmm]: 4.6 to 41 (Standard 10) */
  double ConstantA = PI*PI;
  g_BucklingConst = ConstantA * EI;
  g_k_on = g_Vg; /* or g_k_on = g_Vg/0.8 */
  double F_dependency_single = 3.2e+10; /* 0.034e+12 when velocity decreases to 60% with 15pN (condition resembles in yeast). theoretical value 17.7 (for C. elegans)*/
  // pulling forces
  g_MotorStallF = 1.1e-12;
  g_MotorMaxVel = 2.0e-6;
  // length-dependent pulling force
  g_MotorDensity = 0.050e+6; /* D: density of motor on MT [/m]: 50,000 to 400,000 (standard 100,000)*/
  // corical pulling force
  int CortexPullingDuration = 1; /* [timepoint] */
  double CortexPullingFreq_PAR3 = 0.8; /* [/sec] */
  double CortexPullingFreq_PAR2 = 1.2; /*  */
  double CortexPullingFreq_LET99 = 0.4;
  double CortexPullingFreq;
  // other models
  double Fc = 1.0e-12;

  // arrangement of microtubules: calculation of the (maximum) number of MTs
  double starting_degree = 0.0; /* starting degree of aster: for the simulation with single aster (Sup Fig S9) this value should be changed*/
  //double starting_degree = PI/2.0; /* starting degree of aster: for the simulation with single aster (Sup Fig S9) this value should be changed*/
  int MTDivision = 6; /* this number defines the (maximum) number of MTs*/
  int MTDivisionPlus = 0;
  int MTPerPlane[MTDivision+MTDivisionPlus+1];
  int MTTotal[MTDivision+MTDivisionPlus+1];
  int kk; /* plane number */
  int k; /* MT number */
  MTPerPlane[0]=1;
  MTTotal[0]=1;
  for (kk=1; kk<MTDivision+1; kk++) {
    MTPerPlane[kk]=(int)(0.5+sin(PI*kk/(2.0*MTDivision))*4.0*MTDivision);
    MTTotal[kk]=MTTotal[kk-1]+MTPerPlane[kk];
    printf("%d %d %d\n",MTDivision,MTPerPlane[kk],MTTotal[kk]);
  }
  for (kk=1; kk<MTDivisionPlus+1; kk++) {
    MTPerPlane[MTDivision+kk]=(int)(0.5+sin(PI*(MTDivision-kk)/(2.0*MTDivision))*4.0*MTDivision);
    MTTotal[MTDivision+kk]=MTTotal[MTDivision+kk-1]+MTPerPlane[MTDivision+kk];
    printf("%d %d %d\n",MTDivision,MTPerPlane[MTDivision+kk],MTTotal[MTDivision+kk]);
  }
  g_NN = MTTotal[MTDivision+MTDivisionPlus];
  g_N = 2*g_NN;
  printf("NN=%d N=%d\n", g_NN, g_N);
  // allocate memory
  g_u = dmatrix(0,g_N-1,0,5); ///////////////////////////////
  g_pushing_phase = cvector(0,g_N-1);
  g_pulling_phase = cvector(0,g_N-1);
  g_phase = cvector(0,g_N-1);
  g_L = dvector(0,g_N-1);
  g_NumberOfMotor = dvector(0,g_N-1);
  int CortexPullingMode[g_N];
  int CortexMotorNumber[g_N];
  for (k=0;k<g_N;k++) {CortexPullingMode[k]=0;}
  // arrangement of microtubules: direction of each MT
  double Sdelta[g_N], Ssita[g_N];
  Sdelta[0]=0.5*PI; Ssita[0]=0.0;
  kk=0;
  for (k=1;k<g_NN; k++) {
    if (k > MTTotal[kk]-1){kk++;}
    Sdelta[k]=(PI/2.0)*(1.0-(double)kk/MTDivision);
    Ssita[k]=2.0*PI*k/(MTPerPlane[kk]);
  }
  Sdelta[g_NN] = -0.5*PI; Ssita[0]=0.0;
  kk=0;
  for (k=g_NN+1;k<g_N; k++) {
    if ((k-g_NN) > MTTotal[kk]-1){kk++;}
    Sdelta[k]=-(PI/2.0)*(1.0-(double)kk/MTDivision);
    Ssita[k]=2.0*PI*(k-g_NN)/(MTPerPlane[kk]);
  }
  // int MTcoefficient = 100; /* a constant used in simulations with increasing number of MTs (Sup Fig S5) */

  /************VARIABLES**********************/
  int i; /* steps */
  double t; /* time [sec] */
  int p; /* different parameter sets */
  int j,jj; /* x,y,z axes */
  int qq; /* centrosome 0 or 1 */
  unsigned int mode; /* pulling or pushing */
  double VECVEC[3], VECVECVEC[3];
  // microtubules //
  double previousL[g_N];
  double MT[g_N][3]; /* rectangular coordinates of n-th MT */
  double tempMT[g_N][3];
  double MTC[3]; /* for calculation of length to cortex used around L448 */
  double min, max, sum; /* to monitor max, mean, min of MT length */
  int currentN; /* current number of active MTs */
  // first and second centrosome//
  g_DVecNucCen = dmatrix(0,1,0,2);
  double PVecCen[2][3];
  double Nuc[3];  /* position vector of center of the nucleus */
  // forces //
  double ForceC[2][6]; /* force vector on 1st and 2nd centrosome */
  //double PushingForceC[2][3];
  double Bucklingforce_single;
  double pushingF[g_N];
  // translational and rotational movement of the pronucleus //
  double DirectionDetermination[3];
  double UnitDirection[3];
  double DirectCos[g_N];
  double StallSum;
  double Vnuc_buckle;
  double Vmt_buckle;
  double F_soln;
  double Vnuc_soln;
  double Vmt_soln;
  double *tempNucVel; 
  tempNucVel=dvector(1,6);
  double CalculateVel;
  unsigned int cycle_count;
  unsigned char final_cycle_check;
  unsigned int eachMT_PTC[g_N];
  unsigned int pushing_phase_count[5];
  double dv; /* increase of distance between contact point and nucleus devided by time */
  int step_counter;
  ////////////////////////ROTATION
  double Rotation[3]; /* rotational vector */
  double rotationAx[3]; /* rotational axis */
  ////////////////////// Nucleus does not cross over the cortex /////////////////
  double KKK;
  double NuclearDistanceRatio;
  double OldDis, NewDis, OldVel, NewVel;
  double DistanceFromPP;
  // solution of a quadratic equation //
  double AA[3], AAA[3]; /* three coefficients of the equation */
  double BB[2], BBB[2]; /* two solutions of the equation */
  // variables for local movement
  double DeltaNucStandard[3]={0.0,0.0,0.0};
  double DeltaNuc[3], OldNuc[3], OldNucLength, OldNucPosition[3], UnitOldNucPosition[3], DeltaNucAxis1, DeltaNucAxis2;
  double np, cumpoisson;
  int motornumber;
  int laser_centrosome;
  double CenVel[2][4];
  double dvid[3];
  double RotationMatrix[3][3];

  char filefin[40];
  char filefin2[40];

  // thresholds in Newton-Raphson method
  double xacc = 1.0e-17; /* 1e-5 pN */
  double tolx = 1.0e-12; /* 1e-6 um/sec */
  double tolf = 1.0e-17; /* 1e-5 pN */
  boolean did_converge = false;
  void (*usr_func)(double *x, int n, double *fvec, double **fjac); /* for callback function */

  //// DECLARATION of Constants and Variables - FINISHED //

  /* Parse options */
  myname = argv[0];
  while ((ch = getopt(argc, argv, "t:s:d:m:cv")) != -1){
    switch (ch) {
      /*
         case 't':
         sim_time = atof(optarg);
         break;
         case 's':
         step = atoi(optarg);
         break;
         case 'd':
         delta = atof(optarg);
         break;
         */
      case 'c':
        is_check = true;
        break;
      case 'v':
        is_verbose = true;
        break;
      case 'm':
        mode = atoi(optarg);
        break;
      case 'h':
        usage(myname);
        exit(1);
      default:
        usage(myname);
        exit(1);
    }
  }
  argc -= optind;
  argv += optind;

  if(argc > 1){
    usage(myname);
  }

  /*
     printf("Choose a model: pushing(0), pulling(1), cortical-anchoring(2), pulling_sup1(3), pushing_sup2(4) or length-dependent+cortical-anchoring(5)?:");
     scanf("%d",&mode);
     */

  //////////////////////////////////////////
  // Examination of different parameters ///
  //////////////////////////////////////////

  for (p=0; p<1; p++) {
    switch (p) {
      case 0: /*orange*/
        //CortexPullingFreq_PAR3 = 0.7; /* [/sec] */
        //CortexPullingFreq_PAR2 = CortexPullingFreq_PAR3*1.5;
        //MotorDensity = 0.125e+5;
        // CortexPullingFreq_PAR2 = 1.0/4.0; /* [/sec] */
        //  CortexPullingFreq_PAR3 = 0.5/4.0; /*  */
        //CortexPullingFreq_LET99 = 0.5; /*  */
        //	  g_Visco = 0.25;
        // MTcoefficient = 800;
        // Fc = 4.0e-12;
        //	  NumberOfMotorOnCortex_PAR3=0.0;
        break;
      case 1: /*pink*/
        //CortexPullingFreq_PAR3 = 0.75; /* [/sec] */
        //CortexPullingFreq_PAR2 = CortexPullingFreq_PAR3*1.5;
        //  CortexPullingFreq_PAR2 = 1.0/2.0; /* [/sec] */
        //  CortexPullingFreq_PAR3 = 0.5/2.0; /*  */
        //CortexPullingFreq_LET99 = 0.1; /*  */
        //	  g_Visco = 0.5;
        // MTcoefficient = 400;
        // Fc = 2.0e-12;
        //	  NumberOfMotorOnCortex_PAR3=0.2;
        break;
      case 2: /* black */
        //CortexPullingFreq_PAR3 = 0.8; /* [/sec] */
        //CortexPullingFreq_PAR2 = CortexPullingFreq_PAR3*1.5;
        //  MotorDensity = 0.5e+5;
        //  CortexPullingFreq_PAR2 = 1.0; /* [/sec] */
        // CortexPullingFreq_PAR3 = 0.5; /*  */
        //CortexPullingFreq_LET99 = 0.05; /*  */
        //	  g_Visco = 1.0;
        // MTcoefficient = 200;
        // Fc = 1.0e-12; 
        //	  NumberOfMotorOnCortex_PAR3=0.5;
        break;
      case 3: /* blue */
        //CortexPullingFreq_PAR3 = 0.85; /* [/sec] */
        //CortexPullingFreq_PAR2 = CortexPullingFreq_PAR3*1.5;
        //  MotorDensity = 1.0e+5;
        //  CortexPullingFreq_PAR2 = 1.0*2.0; /* [/sec] */
        //  CortexPullingFreq_PAR3 = 0.5*2.0; /*  */
        //CortexPullingFreq_LET99 = 0.025; /*  */
        //	  g_Visco = 2.0;
        // MTcoefficient = 100;
        // Fc = 0.5e-12; 
        //	  NumberOfMotorOnCortex_PAR3=1.0;
        break;
      case 4: /* green */
        //CortexPullingFreq_PAR3 = 0.9; /* [/sec] */
        //CortexPullingFreq_PAR2 = CortexPullingFreq_PAR3*1.5;
        //  MotorDensity = 2.0e+5;
        //  CortexPullingFreq_PAR2 = 1.0*4.0; /* [/sec] */
        //  CortexPullingFreq_PAR3 = 0.5*4.0; /*  */
        //CortexPullingFreq_LET99 = 0.0; /*  */
        //	  g_Visco = 4.0;
        // MTcoefficient = 50; 
        // Fc = 0.25e-12;
        //	  NumberOfMotorOnCortex_PAR3=1.5;
        break;
    }
    g_Stokes_translation = 6.0*PI*g_Stokes_rad*g_Visco;
    g_Stokes_rotation = -8.0*PI*pow(g_Stokes_rad, 3)*g_Visco/1.0; /* the formula may be adjusted for easy rotation */
    g_BucklingConst = ConstantA * EI;

    // OUTPUT parameter LOGs
    fprintf(f_out5,"p=%d\nmode=%d\nMTDivision=%d MT=%d\nVg_um=%5.3lf Vs_um=%5.3lf\nCatFreq=%5.3lf ResFreq=%5.3lf\nFstall_pN=%5.3lf Vmax_um=%5.3f M_mm=%5.3lf\nEI_pNumum=%5.3lf A_um=%5.3lf B_pN=%5.3lf\nH=%5.3lf Stokes_rad=%5.3lf\n\n", p, mode, MTDivision, g_N, g_Vg*pow(10,6), Vs*pow(10,6), CatFreq, ResFreq, g_MotorStallF*pow(10,12), g_MotorMaxVel*1.0e+6, g_MotorDensity*pow(10,-3), EI*pow(10,24), g_k_on*1.0e+6, F_dependency_single*1.0e-12, g_Visco, g_Stokes_rad*1.0e+6);

    // color settings
    /* #include "color_setting.c" */
    color_setting(&mtg, p);

    /**************** INITIALIZATION *********************/
    // center of the nucleus //
    Nuc[0]=Rad-LL,Nuc[1]=0.0,Nuc[2]=0.0; /**************** start point ********************/
    //Nuc[0]=0.0,Nuc[1]=0.0,Nuc[2]=0.0; /**************** start point ********************/
    g_DVecNucCen[0][0]=0.0; g_DVecNucCen[0][1]=0.0; g_DVecNucCen[0][2]=LL;      
    g_DVecNucCen[1][0]=0.0; g_DVecNucCen[1][1]=0.0; g_DVecNucCen[1][2]=(-1.0)*LL;

    // 1st rotation
    rotationAx[0]=0.0;
    rotationAx[1]=1.0;
    rotationAx[2]=0.0;
    MakeRotationMatrix(RotationMatrix, rotationAx, starting_degree);

    ProductJacVec(VECVEC, RotationMatrix, g_DVecNucCen[0]);
    for (j=0; j<3; j++) {g_DVecNucCen[0][j]=VECVEC[j];}
    ProductJacVec(VECVEC, RotationMatrix, g_DVecNucCen[1]);
    for (j=0; j<3; j++) {g_DVecNucCen[1][j]=VECVEC[j];}

    for (j=0; j<3; j++) {
      PVecCen[0][j] = Nuc[j]+g_DVecNucCen[0][j];
      PVecCen[1][j] = Nuc[j]+g_DVecNucCen[1][j];
    }

    // phase and length of MTs
    for (k=0; k<g_N; k++){
      g_L[k] = Vs*dT;  
      random=ran1(idum);
      if (random<(ResFreq/(CatFreq+ResFreq)))
        g_phase[k]=1;
      else 
        g_phase[k]=0;
      // phase[k]=2; /*** add when simulation with increasing number of MTs (Sup Fig S3), start = inactive phase ***/
      g_u[k][0] = cos(Sdelta[k])*cos(Ssita[k]);
      g_u[k][1] = cos(Sdelta[k])*sin(Ssita[k]);
      g_u[k][2] = sin(Sdelta[k]);
      for (j=0;j<3;j++){VECVEC[j]=g_u[k][j];}
      ProductJacVec(VECVECVEC, RotationMatrix, VECVEC);
      for (j=0;j<3;j++){g_u[k][j]=VECVECVEC[j];}
      // calculation of  g_u[k][3], g_u[k][4], g_u[k]5]
      if (k<g_NN) {
        OutProdVector(g_DVecNucCen[0], VECVECVEC, VECVEC);
      } else {
        OutProdVector(g_DVecNucCen[1], VECVECVEC, VECVEC);
      }
      for (j=0;j<3;j++){g_u[k][3+j]=VECVEC[j];}
      g_pushing_phase[k] = 0;
    }
    currentN = g_N;
    // currentN = 0;  /*** add when simulation with increasing number of MTs (Sup Fig S3), start = 0 MT ***/

    ////////////////////////////////////////
    ///// Repetition for each time step ////
    ////////////////////////////////////////

    step_counter = -1;
    for (i=0; i<ST; i++) {
      step_counter++;
      ////////////////////  INITIALIZATION ////////////////////////////
      for (j=0;j<6;j++) {
        for (qq=0; qq<2; qq++){
          ForceC[qq][j]=0.0;
        }
      }
      for (j=0; j<3;j++) {	  
        DirectionDetermination[j]=0.0;
      }
      max = 0;
      sum = 0;
      min = Rad;

      //  activation of MTs /****** add when simulation with increasing number of MTs (Sup Fig S3) ****/
      // 	  if (currentN<N)
      //	    {
      //	      FT = MTcoefficient*pow(i*dT/420,1);
      //	      increaseN = (int)(FT-currentN);
      //	      for (int ccc=0; ccc<increaseN; ccc++)
      //		{
      //		  do
      //		    {
      //		      random=ran1(idum);
      //		    } while (phase[(int)(random*N)]!=2);
      //		  phase[(int)(random*N)]=1;
      //		  currentN++;
      //		}
      //	      if (i%20==0){
      //		fprintf(f_out6,"%d, %d, %d\n",i,i/20,currentN);}
      //	    }


      /////////////////////  growth/shrinkage for each microtubule ////////////////////
      for (k=0; k<g_N; k++) {
        g_NumberOfMotor[k] = 0.0; /* initialization */
        if (k<g_NN){qq=0;}else{qq = 1;} /* acting on centrosome 1 or 2 */
        g_pushing_phase[k] = 0; /* initial assumption = no pushing */
        switch (g_phase[k]) { /* phase1=growing, phase0=shrinking */
          case 2: /* inactive MTs */
            g_pushing_phase[k] = 0;
            g_pulling_phase[k] = 0;
            break;
          case 1: /* growing phase */
            g_pulling_phase[k] = 1;
            previousL[k] = g_L[k];
            g_L[k] += g_Vg*dT;
            // if (k > NN-1) {g_L[k]=Vs*dT;} /****** add when simulation with 1aster (Sup Fig S8) ******/
            for (j=0; j<3; j++) {
              MTC[j] = g_L[k]*g_u[k][j];
              tempMT[k][j] = PVecCen[qq][j] + MTC[j];}
            if (sqrt(pow((tempMT[k][0]/Rad),2)+pow((tempMT[k][1]/RadS),2)+pow((tempMT[k][2]/RadS),2)) > 1) { /* the MT reaches the cortex */
              // calculation of the distance from the centrosome to the cell cortex at the angle //
              AA[0]=pow(MTC[0]*RadS,2)+pow(MTC[1]*Rad,2)+pow(MTC[2]*Rad,2); /* quadratic equation */
              AA[1]=2*(MTC[0]*PVecCen[qq][0]*RadS*RadS+MTC[1]*PVecCen[qq][1]*Rad*Rad+MTC[2]*PVecCen[qq][2]*Rad*Rad);
              AA[2]=pow(RadS*PVecCen[qq][0],2)+pow(Rad*PVecCen[qq][1],2)+pow(Rad*PVecCen[qq][2],2)-pow(Rad*RadS,2);
              QuadEqu2(AA,BB);
              g_L[k] = BB[0]*g_L[k];

              // THE PUSHING MODEL: estimation of an initial guess to be used in Newton-Raphson method 
              if (mode == 0) {
                //  if (qq=0) { /************add if 1aster ******/ 
                g_pushing_phase[k]=1; /* temporary assume buckling */
                for (j=0; j<3; j++) {
                  tempMT[k][j] = PVecCen[qq][j] + g_L[k]*g_u[k][j]; /* coordinate of contact point */
                  DirectionDetermination[j] -= g_u[k][j]/(g_L[k]*g_L[k]); /* to estimate the direction vector of an initial guess to be used in the Newton-Raphson method */
                }
              } else {
                g_pushing_phase[k]=0;
              }

              // THE CORTICAL-ANCHORING MODEL (Sup Fig S6) 
              if (mode == 2) {
                for (j=0; j<6; j++){
                  ForceC[qq][j] += g_u[k][j] * Fc;} /* To examine the cortical-anchoring model with the additional assumptions (Sup Fig S6), conditions described in the figure legend are added here */
              }

              // PUSHING MODELS WITH DIFFERENT FORMULAS (Sup Fig S2)	
              if (mode == 4) {
                for (j=0; j<6; j++){
                  ForceC[qq][j] -= g_u[k][j] * 1.0e-22 * pow(g_L[k],-2);}  /* To examine the pushing forces expressed with different formulas (Sup Fig S2), the formulas shown here are modified as described in the figure legend */
              }

              // LENGTH-DEPENDENT PULLING + CORTEX-PULLING
              if (mode == 5) {
                if (CortexPullingMode[k]==0) {
                  for (j=0; j<3; j++) {
                    tempMT[k][j] = PVecCen[qq][j] + g_L[k]*g_u[k][j];
                  }
                  switch (strain)
                  {
                    case 0:  // WT
                      if (tempMT[k][0]<0) {
                        CortexPullingFreq = CortexPullingFreq_PAR3;
                      } else {
                        if (i<metaST) {
                          CortexPullingFreq = CortexPullingFreq_LET99+(tempMT[k][0]/Rad)*(CortexPullingFreq_PAR2-CortexPullingFreq_LET99);
                        } else {
                          CortexPullingFreq = CortexPullingFreq_PAR2;
                        }			
                      }
                      break;
                    case 1: // par-2 w/o LET-99
                      CortexPullingFreq = CortexPullingFreq_PAR3;
                      break;
                    case 2: // par-2 with LET-99
                      if ((i<metaST)&&(tempMT[k][0]>0)) {
                        CortexPullingFreq = ((Rad-tempMT[k][0])/Rad)*CortexPullingFreq_PAR3;
                      } else {
                        CortexPullingFreq = CortexPullingFreq_PAR3;
                      }						    
                      break;
                    case 3: // par-3 w/o LET-99
                      CortexPullingFreq = CortexPullingFreq_PAR2;			    
                      break;
                    case 4: // par-3 with LET-99
                      if (i<metaST) {
                        CortexPullingFreq = ((fabs(tempMT[k][0]))/Rad)*CortexPullingFreq_PAR2;
                      } else {
                        CortexPullingFreq = CortexPullingFreq_PAR2;
                      }						    			    
                      break;
                    case 5: // let-99
                      if (tempMT[k][0]<0) {
                        CortexPullingFreq = CortexPullingFreq_PAR3;
                      } else {
                        CortexPullingFreq = CortexPullingFreq_PAR2;
                      }						    			    
                      break;
                    case 6: // ric-8
                      if (tempMT[k][0]<0) {
                        CortexPullingFreq = CortexPullingFreq_PAR3;
                      } else {
                        CortexPullingFreq = (tempMT[k][0]/Rad)*CortexPullingFreq_PAR2;
                      }			
                      break;
                    case 7: // let-99; ric-8
                      if (tempMT[k][0]<0) {
                        CortexPullingFreq = CortexPullingFreq_PAR3;
                      } else {
                        CortexPullingFreq = CortexPullingFreq_PAR2;
                      }						    			  
                      break;
                    case 8: // gpr-1/2
                      CortexPullingFreq=0.0;
                      break;
                    case 9:
                      if (tempMT[k][0]>0) {
                        CortexPullingFreq = CortexPullingFreq_PAR2;
                      } else {CortexPullingFreq = 0.0;}						    		
                      break;
                  }
                  random = ran1(idum);
                  np = CortexPullingFreq;
                  cumpoisson = 0.0;
                  motornumber = -1;
                  do {
                    motornumber++;
                    cumpoisson += Poisson(np,motornumber);
                  } while (random > cumpoisson);
                  CortexMotorNumber[k] = motornumber;

                  //			random = ran1(idum);
                  //np = CortexPullingDuration/dT;
                  //cumpoisson = 0.0;
                  //motornumber = -1;
                  //do {
                  //  motornumber++;
                  //  cumpoisson += Poisson(np,motornumber);
                  //} while (random > cumpoisson);
                  //CortexPullingMode[k] = motornumber;

                  CortexPullingMode[k]=CortexPullingDuration; /********** no duration ***************/
                }
                if (CortexPullingMode[k]>0){
                  g_NumberOfMotor[k] += CortexMotorNumber[k];
                  CortexPullingMode[k]--;
                } 
              }

              // if (ran1(idum) < CatFreq*20*dT){phase[k] = 0;} /* if assuming CatFreq increases upon contact with the cortex */
              } else {
                g_pushing_phase[k]=0; // To ensure that pushing_phase of the MT without contacting the cortex is 0 // 
              }
              // switching from growth phase to shrinkage phase // 
              if (ran1(idum) < CatFreq*dT)
                g_phase[k] = 0;
              break;

              case 0: /* shrinkage phase */
              if (g_L[k] > Rad/25) {
                g_L[k] -= Vs*dT;
                g_pulling_phase[k]=1;
              }  /* the lengh of MT exceeds 1micron */
              //  if (k > NN-1) {g_L[k]=Vs*dT;} /************* add if 1aster (Sup Fig. S8)********************/

              // switching from shrinkage phase to growth phase // 
              if (ran1(idum) < ResFreq*dT) g_phase[k] = 1;
              g_pushing_phase[k]=0;
              break;
            }

            // PULLING MODELS WITH DIFFERENT FORMULAS (Sup Fig S1)
            if (mode == 3) { 
              for (j=0;j<6;j++){
                ForceC[qq][j] += g_u[k][j] * 1.0e-7 * pow(g_L[k],1);} /* pulling models with different fomulas (Sup Fig S1), the formulas below are modified as described in the figure legend */  
            }

            ///////////////////////////////////////////
            // LENGTH-DEPENDENT PULLING + CORTEX PULLING
            if ((mode == 1)||(mode == 5)) {
              random = ran1(idum);
              np = g_L[k]*g_MotorDensity;
              cumpoisson = 0.0;
              motornumber = -1;
              do {
                motornumber++;
                cumpoisson += Poisson(np,motornumber);
              } while (random > cumpoisson);
              g_NumberOfMotor[k] += motornumber;
              //TRACE(("%4d %4d N=%3.1lf MT=(%3.2lf %3.2lf %3.2lf %3.2lf %3.2lf %3.2lf)\n",i,k,NumberOfMotor[k],u[k][0],u[k][1],u[k][2],u[k][3],u[k][4],u[k][5]));
              for (j=0;j<6;j++){
                ForceC[qq][j] += g_NumberOfMotor[k] * g_u[k][j] * g_MotorStallF;
              } /* an initial guess */
            }

            // to monitor microtubules profile
            sum += g_L[k];
            if (max < g_L[k]){max = g_L[k];}
            if (min > g_L[k]){min = g_L[k];}
        }

        // for (j=0; j<6; j++) {ForceC[1][j] = 0.0;} /********* add if 1 aster (Sup Fig S8) ******/

        ///////////////// Translational and rotational movement of the pronucleus //////////////////////

        //for calculation of the velocity of pronucleus 
        if (i%UNITTIME==0){
          if (i==0){
            OldDis = LL;
            OldVel = 0.0;
          } else {
            OldDis = NewDis;
            OldVel = NewVel;
          }
        }
        // initialization
        for (j=1; j<=6; j++){
          tempNucVel[j]=0.0;
        }

        // THE PUSHING MODEL-1: estimation of an initial guess to be used in Newton-Raphson method 
        if ((mode==0)&&(Length(DirectionDetermination)!=0.0)) {
          usr_func = function_FV3D;
          if (i%100==0) {fprintf(f_out7,"%d ",i);}
          UnitVector(DirectionDetermination,UnitDirection); /* unit vector of an initial guess of pronuclear migration */
          g_Buckling_forward_sum = 0.0;
          g_Buckling_backward_sum = 0.0;
          StallSum=0.0;
          for (k=0; k<g_N; k++) {
            if (g_pushing_phase[k]!=0) {
              DirectCos[k]=InnProdVector(g_u[k],UnitDirection);
              if (DirectCos[k]<0) {
                g_Buckling_forward_sum -= DirectCos[k]*g_BucklingConst/(g_L[k]*g_L[k]);
                StallSum -= DirectCos[k];
                g_pushing_phase[k]=2; /* temporal assumption that MTs directed toward the opposite direction of pronuclear migration and pushing the cortex do not buckle */
              } else {
                g_Buckling_backward_sum += DirectCos[k]*g_BucklingConst/(g_L[k]*g_L[k]);
                g_pushing_phase[k]=1; /* temporal assumption that MTs directed toward the direction of pronuclear migration and pushing the cortex buckle  */
              }
            }
          }
          if (g_Buckling_forward_sum < g_Buckling_backward_sum) {
            printf("something is wrong: cannot determine buckling direction!\n");
          } else {
            g_F_dependency = F_dependency_single/StallSum;
            Vnuc_buckle = Stokes_function(g_Buckling_forward_sum, g_Stokes_translation, g_Buckling_backward_sum);
            Vmt_buckle = FV_function(g_Buckling_forward_sum, g_Vg, g_k_on, g_F_dependency);
            if (Vnuc_buckle<=Vmt_buckle){  /* buckling is dominant */
              if (i%100==0) {fprintf(f_out7,"BK ");}
              F_soln = g_Buckling_forward_sum;
              Vnuc_soln = Vnuc_buckle;
              Vmt_soln = Vmt_buckle;
              for (k=0;k<g_N;k++) {
                if (g_pushing_phase[k]!=0) g_pushing_phase[k]=1; /* temporal assumption that all MTs buckle */
              }
            } else {  /* FV is dominant */
              if (i%100==0) {
                fprintf(f_out7,"FV ");
              }
              F_soln = rtsafe(FV_solution, 0.0, g_Buckling_forward_sum, xacc);
              Vnuc_soln = Stokes_function(F_soln, g_Stokes_translation, g_Buckling_backward_sum);
              Vmt_soln = FV_function(F_soln, g_Vg, g_k_on, g_F_dependency);
            }
          }
          for (j=1; j<=3; j++){
            tempNucVel[j]=Vnuc_soln*UnitDirection[j-1]; /* an initial guess of the velocity of the pronucleus to be used in Newton-Raphson method*/
            if (i%100==0) {fprintf(f_out7,"%5.4f ",tempNucVel[j]*pow(10,6));}
          }

          g_F_dependency = F_dependency_single;
          for (j=0;j<3;j++){
            g_Fbuckle[j] = 0.0;  /* Fbuckle[0,1,2]*/
          }
          for (k=0; k<g_N; k++) {
            if (g_pushing_phase[k]==1){
              Bucklingforce_single = g_BucklingConst/(g_L[k]*g_L[k]);	  
              for (j=0;j<3;j++){
                g_Fbuckle[j] -= g_u[k][j]*Bucklingforce_single;
              }
            }
          }
          usr_func = function_FV3D;
          /* Newton-Raphson method to revise the initial guess of the velocity of the pronucleus */
          did_converge = mnewt(10,tempNucVel,3,tolx,tolf, step_counter, f_out8, usr_func);

          // THE PUSHING MODEL-2: solve the set of equation using Newton-Raphson method with the revised initial guess
          cycle_count = 0;
          final_cycle_check = 0;
          for (k=0; k<g_N; k++) {eachMT_PTC[k] = 0;}
          do {
            g_phase_transition_count = 0;
            for (j=0;j<3;j++){
              g_Fbuckle[j] = 0.0;  /* Fbuckle[0,1,2]*/
            }
            for (k=0; k<g_N; k++) {
              if (g_pushing_phase[k]!=0){
                Bucklingforce_single = g_BucklingConst/(g_L[k]*g_L[k]);	  
                Vmt_buckle = FV_function(Bucklingforce_single, g_Vg, g_k_on, F_dependency_single);
                dv = -1.0 * (g_u[k][0]*tempNucVel[1]+g_u[k][1]*tempNucVel[2]+g_u[k][2]*tempNucVel[3]);
                if (eachMT_PTC[k]>10) {  /* to avoid examining same condition repeatedly */
                  if ((dv > g_Vg) || (dv < Vmt_buckle)||(g_pushing_phase[k]!=2)) { /* check the length. If not appropriate, do not exit */
                    g_phase_transition_count++;
                    eachMT_PTC[k]++;
                  }
                  g_pushing_phase[k] = 2;
                  if (eachMT_PTC[k]>25){ eachMT_PTC[k] = 0; } /* if phase 2 is not appropriate, go back to usual classification */
                } else {
                  if (dv > g_Vg){
                    if (g_pushing_phase[k]!=3){
                      g_phase_transition_count++;
                      eachMT_PTC[k]++;
                    }
                    if (g_pushing_phase[k]==1){
                      g_pushing_phase[k] = 2;
                    } else {
                      g_pushing_phase[k] = 3;
                    }
                  } else {
                    if (dv > Vmt_buckle) {  /* FV is dominant */		
                      if (g_pushing_phase[k]!=2){
                        g_phase_transition_count++;
                        eachMT_PTC[k]++;
                      }
                      g_pushing_phase[k] = 2;
                    } else {
                      if (g_pushing_phase[k]!=1){
                        g_phase_transition_count++;
                        eachMT_PTC[k]++;
                      }
                      if (g_pushing_phase[k]==3) {
                        g_pushing_phase[k] = 2;
                      } else {
                        g_pushing_phase[k] = 1;
                        for (j=0;j<3;j++){
                          g_Fbuckle[j] -= g_u[k][j]*Bucklingforce_single;
                        }
                      }
                    }    
                  }
                }
              }		
            }	      
            if ((g_phase_transition_count!=0)||(cycle_count==0)) {
              /* Newton-Raphson method */
              did_converge = mnewt(10, tempNucVel, 3, tolx,tolf, step_counter, f_out8, usr_func);
            }
            if (i%100==0) {
              fprintf(f_out8,"%d %d %d\n", i, cycle_count, g_phase_transition_count);
            }
            cycle_count++;
          } while ((!did_converge)||(cycle_count<=1)||((g_phase_transition_count!=0)&&(cycle_count<1000))); /* repeat until the solution satisfies all equations and conditions */

          if (g_phase_transition_count!=0){ /* in case the solution is not obtained within 1000 cycles*/
            printf("exit without convergence at t=%d PTC=%d\n", i, g_phase_transition_count);
            fprintf(f_out7,"exit without convergence at t=%d PTC=%d\n", i, g_phase_transition_count);
            printf("Cen1=(%3.1f,%3.1f,%3.1f) Cen2=(%3.1f,%3.1f,%3.1f)\n",PVecCen[0][0]*1.0e+6,PVecCen[0][1]*1.0e+6,PVecCen[0][2]*1.0e+6,PVecCen[1][0]*1.0e+6,PVecCen[1][1]*1.0e+6, PVecCen[1][2]*1.0e+6);
            for (k=0; k<g_N; k++){
              if (g_pushing_phase[k]!=0){
                printf("%d %d L=%3.1f (%3.2f,%3.2f,%3.2f)\n",k,g_pushing_phase[k],g_L[k]*1000000,g_u[k][0],g_u[k][1],g_u[k][2]);
              }
            }
            break;
          }

          // THE PUSHING MODEL-3: calculation of the force exerted by each MT 
          if (i%100==0) {fprintf(f_out7,"%d ",cycle_count);}
          for (j=1; j<=3; j++){
            if (i%100==0){fprintf(f_out7,"%5.4f ",tempNucVel[j]*pow(10,6));}
          }
          for (j=0;j<3;j++) {
            for (qq=0; qq<2; qq++) {
              ForceC[qq][j] = 0.0;}
            pushing_phase_count[j+1] = 0;
          }
          for (k=0; k<g_N; k++) {
            if (k<g_NN) {qq=0;} else {qq=1;}
            if (g_pushing_phase[k]==3) {
              // g_L[k] = previousL[k] + g_Vg*dT; /* this line can be omitted because of redundancy*/
              pushingF[k] = 0.0;
              pushing_phase_count[3]++;
            } 
            if (g_pushing_phase[k]==2) {
              dv = -1 * (g_u[k][0]*tempNucVel[1]+g_u[k][1]*tempNucVel[2]+g_u[k][2]*tempNucVel[3]);
              // g_L[k] = g_L[k]+dv*dT; /* this line can be omitted because of redundancy*/
              pushingF[k] = -1*log(1-(g_Vg - dv)/g_k_on)/F_dependency_single;
              for (j=0;j<3;j++) {
                ForceC[qq][j] -= pushingF[k]*g_u[k][j];}
              pushing_phase_count[2]++;
            }
            if (g_pushing_phase[k]==1) {
              pushingF[k] = g_BucklingConst/(g_L[k]*g_L[k]);
              // g_L[k] = previousL[k]+dT*FV_function(pushingF[k], g_Vg, g_k_on, F_dependency_single); /* omit this line assuming buckling MTs do not elongate any more*/
              for (j=0;j<3;j++) {
                ForceC[qq][j] -= pushingF[k]*g_u[k][j];}
              pushing_phase_count[1]++;
            }
          }
          for (j=1; j<=3; j++){
            if (i%100==0){
              fprintf(f_out7,"%d ",pushing_phase_count[j]);
            }
          }
          if (i%100==0){
            fprintf(f_out7,"\n");
          }
        }

        ///////////////////////////////////////////////////////////////////////
        // THE PULLING MODEL: solve the set of equation using Newton-Raphson method with the initial guess
        if ((mode == 1)||(mode == 5)) {
          if (i<laserST){  /*************** before ablation *******************/
            for (j=1;j<=6;j++){ /* initial value of tempNucVel[j] */
              TRACE(("%4d Force[%d]: %10lf %10lf\n",i, j, ForceC[0][j-1]*1.0e+12, ForceC[1][j-1]*1.0e+12));
              if (j<=3) {
                tempNucVel[j] = (ForceC[0][j-1]+ForceC[1][j-1])/g_Stokes_translation; 
              } else {
                tempNucVel[j] = (ForceC[0][j-1]+ForceC[1][j-1])/g_Stokes_rotation;
              }
            }
            usr_func = function_MotorFV;
            cycle_count = 0;
            TRACE(("%4d init ",i));
            for (j=1; j<=6; j++) {
              if (j<=3) {
                TRACE(("%4.3lf ",tempNucVel[j]*1.0e+6));
              } else {
                TRACE(("%4.3lf ",tempNucVel[j]*100));
              }
            }
            TRACE(("\n"));
            for (k=0; k<g_N; k++) {eachMT_PTC[k]=0;}
            do {
              g_phase_transition_count = 0;
              for (j=0;j<6;j++){
                g_Fbackward[j] = 0.0;
                for (jj=0; jj<6; jj++) {
                  g_fjac_pull[j][jj] = 0.0;
                }}
              for (k=0; k<g_N; k++) {
                if (g_pulling_phase[k]!=0){
                  if (k<g_NN) {qq = 0;} else {qq = 1;}
                  dv = 0.0;
                  for (j=0; j<3; j++) {
                    dv += (tempNucVel[(j+2)%3+4]*g_DVecNucCen[qq][(j+1)%3]-tempNucVel[(j+1)%3+4]*g_DVecNucCen[qq][(j+2)%3]+tempNucVel[j+1])*g_u[k][j];
                  }
                  if (dv > g_MotorMaxVel){ /* the motors on this MT do not exert forces: pulling_phase[k]=3 */
                    if (g_pulling_phase[k]!=3) {
                      g_phase_transition_count++;
                      eachMT_PTC[k]++;
                    }
                    g_pulling_phase[k] = 3;
                  } else {
                    if (dv < 0) { /* the motors on this MT exert the maximum (stall) force: pulling_phase[k]=2 */
                      if (g_pulling_phase[k]!=2) {
                        g_phase_transition_count++;
                        eachMT_PTC[k]++;
                      }
                      g_pulling_phase[k] = 2;
                      for (j=0;j<6;j++) {
                        g_Fbackward[j] += g_u[k][j] * g_MotorStallF * g_NumberOfMotor[k];
                      }
                    } else {
                      if (g_pulling_phase[k]!=1) { /* the motors on this MT exert force dependent on their velocity: pulling_phase[k]=1 */
                        g_phase_transition_count++;
                        eachMT_PTC[k]++;
                      }
                      g_pulling_phase[k] = 1;
                      for (j=0; j<3; j++) {dvid[j] = g_u[k][j];}
                      dvid[3] = g_DVecNucCen[qq][2]*g_u[k][1] - g_DVecNucCen[qq][1]*g_u[k][2];
                      dvid[4] = g_DVecNucCen[qq][0]*g_u[k][2] - g_DVecNucCen[qq][2]*g_u[k][0];
                      dvid[5] = g_DVecNucCen[qq][1]*g_u[k][0] - g_DVecNucCen[qq][0]*g_u[k][1];
                      for (j=0; j<6; j++) {
                        for (jj=0; jj<6; jj++) {
                          g_fjac_pull[j][jj] -= (g_MotorStallF*g_NumberOfMotor[k]/g_MotorMaxVel)*dvid[jj]*g_u[k][j];
                        }
                      }
                    }
                  }
                }
              }
              for (j=0; j<3; j++) {g_fjac_pull[j][j] -= g_Stokes_translation;}
              for (j=3; j<6; j++) {g_fjac_pull[j][j] -= g_Stokes_rotation;}
              if ((g_phase_transition_count!=0)||(cycle_count==0)) {
                did_converge = mnewt(10, tempNucVel, 6, tolx, tolf, step_counter, f_out8, usr_func);
              }
              if (i%100==0) fprintf(f_out8,"%d %d %d\n", i, cycle_count, g_phase_transition_count);
              cycle_count++;
            } while ((cycle_count<=1)||((g_phase_transition_count!=0)&&(cycle_count<1000))); /* repeat until the solution satisfies all equations and conditions */

            TRACE(("%4d %4d ",i,cycle_count));
            for (j=1; j<=6; j++) {
              if (j<=3) {
                TRACE(("%4.3lf ",tempNucVel[j]*1.0e+6));
              } else {
                TRACE(("%4.3lf ",tempNucVel[j]*100));
              }
            }
            TRACE(("\n"));

            if (g_phase_transition_count!=0){ /* when the solution is not obtained within 1000 cycles*/
              printf("exit without convergence at t=%d PTC=%d\n", i, g_phase_transition_count);
              fprintf(f_out7,"exit without convergence at t=%d PTC=%d\n", i, g_phase_transition_count);
              printf("Cen1=(%3.1f,%3.1f,%3.1f) Cen2=(%3.1f,%3.1f,%3.1f)\n",PVecCen[0][0]*1.0e+6,PVecCen[0][1]*1.0e+6,PVecCen[0][2]*1.0e+6,PVecCen[1][0]*1.0e+6,PVecCen[1][1]*1.0e+6, PVecCen[1][2]*1.0e+6);
              break;
            }
          } else {  /******** laser ablation *****************/
            usr_func = function_laserMotorFV;
            for (laser_centrosome = 0; laser_centrosome<2; laser_centrosome++){
              if (laser_centrosome==0) {
                g_mt_start = 0;
                g_mt_end = g_NN;
              } else {
                g_mt_start = g_NN;
                g_mt_end = g_N;
              }
              for (j=1;j<=3;j++){
                tempNucVel[j] = (ForceC[laser_centrosome][j-1])/g_Stokes_translation; /* forces do not include dT **/
              }
              cycle_count = 0;
              for (k=g_mt_start; k<g_mt_end; k++) {eachMT_PTC[k] = 0;}
              do {
                g_phase_transition_count = 0;
                for (j=0;j<3;j++){
                  g_Fbackward[j] = 0.0;  /* Fbackward[0,1,2]*/
                }
                for (k=g_mt_start; k<g_mt_end; k++) {
                  if (g_pulling_phase[k]!=0){
                    dv = (g_u[k][0]*tempNucVel[1]+g_u[k][1]*tempNucVel[2]+g_u[k][2]*tempNucVel[3]);
                    if (dv > g_MotorMaxVel){ /* the motors on this MT do not exert forces */
                      if (g_pulling_phase[k]!=3) {
                        g_phase_transition_count++;
                        eachMT_PTC[k]++;
                      }
                      g_pulling_phase[k] = 3;
                    } else {
                      if (dv < 0) { /* the motors on this MT exert the maximum (stall) force */
                        if (g_pulling_phase[k]!=2) {
                          g_phase_transition_count++;
                          eachMT_PTC[k]++;
                        }
                        g_pulling_phase[k] = 2;
                        for (j=0;j<3;j++) {
                          g_Fbackward[j] += g_u[k][j] * g_MotorStallF * g_NumberOfMotor[k];
                        }
                      } else {
                        if (g_pulling_phase[k]!=1) { /* the motors on this MT exert force dependent on their velocity */
                          g_phase_transition_count++;
                          eachMT_PTC[k]++;
                        }
                        g_pulling_phase[k] = 1;
                      }
                    }
                  }
                }
                if ((g_phase_transition_count!=0)||(cycle_count==0)) {
                  did_converge = mnewt(10, tempNucVel, 3, tolx, tolf, step_counter, f_out8, usr_func);
                }
                if (i%100==0) {fprintf(f_out8,"%d %d %d\n", i, cycle_count, g_phase_transition_count);}
                cycle_count++;
              } while ((cycle_count<=1)||((g_phase_transition_count!=0)&&(cycle_count<1000))); /* repeat until the solution satisfies all equations and conditions */

              if (g_phase_transition_count!=0){ /* when the solution is not obtained within 1000 cycles*/
                printf("exit without convergence at t=%d PTC=%d\n", i, g_phase_transition_count);
                fprintf(f_out7,"exit without convergence at t=%d PTC=%d\n", i, g_phase_transition_count);
                //  printf("Cen1=(%3.1f,%3.1f,%3.1f) Cen2=(%3.1f,%3.1f,%3.1f)\n",Cen[0][0]*1.0e+6,Cen[0][1]*1.0e+6,Cen[0][2]*1.0e+6,Cen[1][0]*1.0e+6,Cen[1][1]*1.0e+6, Cen[1][2]*1.0e+6);
                break;
              }
              for (j=1; j<=3; j++){
                CenVel[laser_centrosome][j] = tempNucVel[j];
              }
            }
          }

          // THE PULLING MODEL: calculation of the force exerted by each MT 
          //	    if (i%100==0) {fprintf(f_out7,"%d ",cycle_count);}
          //	    for (j=1; j<=6; j++){
          //	      if (i%100==0){fprintf(f_out7,"%5.4f ",tempNucVel[j]*pow(10,6));}
          //	    }
          //	    for (j=0;j<6;j++) {
          //	      for (qq=0; qq<2; qq++) {
          //		ForceC[qq][j]=0.0;}
          //	      pulling_phase_count[j+1]=0;
          //	    }
          //	    for (k=0; k<N; k++) {
          //	      if (k<NN) {qq=0;} else {qq=1;}
          //	      if (pulling_phase[k]==3) {
          //		pullingF[k]=0.0;
          //		pulling_phase_count[3]++;
          //	      } 
          //	      if (pulling_phase[k]==2) {
          //		pullingF[k] = g_MotorStallF;
          //	for (j=0;j<6;j++) {
          //		  ForceC[qq][j] += pullingF[k]*NumberOfMotor[k]*u[k][j];}
          //		pulling_phase_count[2]++;
          //	      }
          //	      if (pulling_phase[k]==1) {
          //		dv = (u[k][0]*tempNucVel[1]+u[k][1]*tempNucVel[2]+u[k][2]*tempNucVel[3]); ////////////////////CHANGE/////
          //		pullingF[k]=g_MotorStallF*(1-dv/g_MotorMaxVel);
          //		for (j=0;j<6;j++) {
          //		  ForceC[qq][j] += pullingF[k]*NumberOfMotor[k]*u[k][j];}
          //		pulling_phase_count[1]++;
          //	      }
          //	    }
          //	    for (j=1; j<=3; j++){
          //	      if (i%100==0){fprintf(f_out7,"%d ",pulling_phase_count[j]);}
          //	    }
          //	    if (i%100==0){fprintf(f_out7,"\n");
        }

        // Monitoring local migration of pronucleus
        if ((i%UNITTIME==0)&&(i<metaST)) {
          if (i/UNITTIME>0) {
            for (j=0;j<3;j++) {
              DeltaNuc[j] = Nuc[j]-OldNuc[j];
              OldNucPosition[j] = DeltaNucStandard[j]-OldNuc[j];
            }
            OldNucLength = sqrt(OldNucPosition[0]*OldNucPosition[0]+OldNucPosition[2]*OldNucPosition[2]);
            if ((OldNucLength>Rad*0.24)&&(OldNucLength<Rad*0.72)) {
              UnitOldNucPosition[0] = OldNucPosition[0]/OldNucLength;
              UnitOldNucPosition[2] = OldNucPosition[2]/OldNucLength;
              DeltaNucAxis1 = DeltaNuc[0]*UnitOldNucPosition[0]+DeltaNuc[2]*UnitOldNucPosition[2];
              DeltaNucAxis2 = DeltaNuc[0]*UnitOldNucPosition[2]-DeltaNuc[2]*UnitOldNucPosition[0];
              //draw graph in window 9
              XFillRectangle(mtg.d,mtg.w9,mtg.gc9,(int)(WIN_WIDTH/2+500*1000000*DeltaNucAxis2/(dT*UNITTIME))-1,(int)(WIN_HEIGHT/2-500*1000000*DeltaNucAxis1/(dT*UNITTIME))-1,2,2);
            }
          }
          for (j=0;j<3;j++) {
            OldNuc[j] = Nuc[j];
          }
        }

        // TRANSLATIONAL MOVEMENT OF THE PRONUCLEUS
        CalculateVel = 0.0;
        for (j=1; j<=3; j++) {
          if (((mode==1)||(mode==5))&&(i>=laserST)) {
            tempNucVel[j] = (CenVel[0][j]+CenVel[1][j])/2;
            PVecCen[0][j-1] += CenVel[0][j]*dT;
            PVecCen[1][j-1] += CenVel[1][j]*dT;
          } else {
            if ((mode!=0)&&(mode!=1)&&(mode!=5)){ /* models other than THE PUSHING and PULLING MODELS */
              tempNucVel[j] = (ForceC[0][j-1]+ForceC[1][j-1])/g_Stokes_translation;
            }
          }
          CalculateVel += tempNucVel[j]*tempNucVel[j];
          Nuc[j-1] += tempNucVel[j]*dT;
        }

        // the pronucleus does not cross over the cell cortex 
        NuclearDistanceRatio = sqrt(pow((Nuc[0]/(Rad-LL)),2)+pow((Nuc[1]/(RadS-LL)),2)+pow((Nuc[2]/(RadS-LL)),2));
        if (NuclearDistanceRatio > 1) { /* when the pronucleus contacts the cell cortex */
          if (mode==0) {TRACE(("nucleus crosses over the cortex!!\n"));}
          else {
            AAA[0] = 0.0;
            AAA[1] = 0.0;
            AAA[2] = 0.0;
            for (j=0; j<1; j++){
              Nuc[j] -= tempNucVel[j+1]*dT; /* step -1 */
              AAA[0] += pow((tempNucVel[j+1]*dT)/(Rad-LL),2); /* quadratic equation */
              AAA[1] += 2*(tempNucVel[j+1]*dT)*Nuc[j]/((Rad-LL)*(Rad-LL));
              AAA[2] += pow(Nuc[j]/(Rad-LL),2);
            }
            for (j=1; j<3; j++){
              Nuc[j] -= tempNucVel[j+1]*dT; /* step -1 */
              AAA[0] += pow((tempNucVel[j+1]*dT)/(RadS-LL),2); /* quadratic equation */
              AAA[1] += 2*(tempNucVel[j+1]*dT)*Nuc[j]/((RadS-LL)*(RadS-LL));
              AAA[2] += pow(Nuc[j]/(RadS-LL),2);
            }
            AAA[2] -= 1;
            TRACE(("QuadEqu2 at L1276\n"));
            QuadEqu2(AAA,BBB);
            for (j=0; j<3; j++){
              Nuc[j] += BBB[0]*tempNucVel[j+1]*dT;
            }
            KKK = 0;
            for (j=0; j<1; j++){			  
              KKK += (-1)*Nuc[j]*(1-BBB[0])*(tempNucVel[j+1]*dT)/((Rad-LL)*(Rad-LL));
            }
            for (j=1; j<3; j++){			  
              KKK += (-1)*Nuc[j]*(1-BBB[0])*(tempNucVel[j+1]*dT)/((RadS-LL)*(RadS-LL));
            }
            KKK = KKK/(pow(Nuc[0],2)/pow((Rad-LL),4)+pow(Nuc[1],2)/pow((RadS-LL),4)+pow(Nuc[2],2)/pow((RadS-LL),4));
            for (j=0; j<1; j++){
              Nuc[j] += (1-BBB[0])*tempNucVel[j+1]*dT+KKK*Nuc[j]/pow((Rad-LL),2);}
            for (j=1; j<3; j++){
              Nuc[j] += (1-BBB[0])*tempNucVel[j+1]*dT+KKK*Nuc[j]/pow((RadS-LL),2);}
          }
        }

        // calculation of velocity of the pronucleus
        if(i%UNITTIME==0){
          NewDis = sqrt((Rad-Nuc[0])*(Rad-Nuc[0])+Nuc[1]*Nuc[1]+Nuc[2]*Nuc[2]);
          NewVel = (NewDis-OldDis)/(dT*UNITTIME);
        }
        DistanceFromPP = sqrt((Rad-Nuc[0])*(Rad-Nuc[0])+Nuc[1]*Nuc[1]+Nuc[2]*Nuc[2]);

        if (((mode==1)||(mode==5))&&(i<laserST)) {
          // ROTATIONAL MOVEMENT OF THE PRONUCLEUS
          Rotation[0] = tempNucVel[4];
          Rotation[1] = tempNucVel[5];
          Rotation[2] = tempNucVel[6];
          MakeRotationMatrix(RotationMatrix, Rotation, dT);
          TRACE(("RotationMatrix\n"));
          for (j=0; j<3; j++) {
            for (jj=0; jj<3; jj++) {
              TRACE(("%4.3lf ",RotationMatrix[j][jj]));
            }
            TRACE(("\n"));
          }
          ProductJacVec(VECVEC, RotationMatrix, g_DVecNucCen[0]);
          for (j=0; j<3; j++) {g_DVecNucCen[0][j] = VECVEC[j];}
          ProductJacVec(VECVEC, RotationMatrix, g_DVecNucCen[1]);
          for (j=0; j<3; j++) {g_DVecNucCen[1][j] = VECVEC[j];}
          for (j=0; j<3; j++) {
            PVecCen[0][j] = Nuc[j] + g_DVecNucCen[0][j];
            PVecCen[1][j] = Nuc[j] + g_DVecNucCen[1][j];
          }

          // phase and length of MTs
          for (k=0; k<g_N; k++){
            if (k<g_NN) {qq = 0;} else {qq = 1;}
            for (j=0;j<3;j++){VECVEC[j] = g_u[k][j];}
            ProductJacVec(VECVECVEC, RotationMatrix, VECVEC);
            for (j=0;j<3;j++){
              g_u[k][j] = VECVECVEC[j];
              MT[k][j] = PVecCen[qq][j] + g_L[k]*g_u[k][j];
            }
            // calculation of  g_u[k][3], g_u[k][4], g_u[k]5]
            OutProdVector(g_DVecNucCen[qq], VECVECVEC, VECVEC);
            for (j=0;j<3;j++){g_u[k][3+j] = VECVEC[j];}
          }

          if (i%100==0){
            fprintf(f_out7,"\n");
            for (j=1; j<=6; j++){
              fprintf(f_out7,"%4.3f ",tempNucVel[j]*(1.0e+6));
            }	      
            fprintf(f_out7,"\n%d\n",i+1);
          }
        } else {
          if (((mode==1)||(mode==5))&&(i>=laserST)) {
            //	      printf("laser!!\n");
            for (k=0; k<g_N; k++){
              if (k<g_NN) {qq = 0;} else {qq = 1;}
              for (j=0; j<3; j++) {
                MT[k][j] = PVecCen[qq][j] + g_L[k]*g_u[k][j];
              }	
            }	    
          } else {
            if (t == 0) {
              TRACE(("L1373 no rotation!\n"));
            }
            for (j=0; j<3; j++) {
              PVecCen[0][j] = Nuc[j] + g_DVecNucCen[0][j];
              PVecCen[1][j] = Nuc[j] + g_DVecNucCen[1][j];
            }
            for (k=0; k<g_N; k++){
              if (k<g_NN) {qq = 0;} else {qq = 1;}
              for (j=0; j<3; j++) {
                MT[k][j] = PVecCen[qq][j] + g_L[k]*g_u[k][j];
              }	
            }
          }
        }
        ///////////////////////// 1 STEP FINISHED ////////////////////////////////////////////
        t += dT;
        //OUTPUTS 
        /* draw graphs in X-window */
        draw_graphs(i, &mtg, PVecCen, MT, Nuc, DistanceFromPP, min, max, sum, currentN, OldVel, NewVel);
        /* save logs in texts */
        save_logs(i, p, g_N, f_out1, f_out2, f_out3, f_out4, f_out10, data_for_3D, Nuc, PVecCen, MT);
        }
        //////////////// examination with single parameter set FINISHED ///////////////////////////
        XStoreName(mtg.d,mtg.w1,"fin");
        sprintf (filefin, "xwd -name \'fin\' -out outFIN%d",p);
        system (filefin);
        sprintf (filefin2, "convert outFIN%d outFIN%d.jpg",p,p);
        system (filefin2);
      }
      ///////////// examination with all parameter sets FINISHED ///////////////////
      fclose(f_out1);
      fclose(f_out2);
      fclose(f_out3);
      fclose(f_out4);
      fclose(f_out10);
      fclose(f_out5);
      fclose(data_for_3D);
      fclose(f_out6);
      fclose(f_out7);
      fclose(f_out8);
      XStoreName(mtg.d,mtg.w1,"fin");
      XStoreName(mtg.d,mtg.w2,"path");
      XStoreName(mtg.d,mtg.w3,"distance");
      XStoreName(mtg.d,mtg.w4,"poledistance");
      XStoreName(mtg.d,mtg.w5,"length");
      XStoreName(mtg.d,mtg.w6,"single");
      XStoreName(mtg.d,mtg.w7,"MTnum");
      XStoreName(mtg.d,mtg.w8,"velocity");
      XStoreName(mtg.d,mtg.w9,"vector");
      XFlush (mtg.d);
      /* save graphs */
      store_graphs();

      free_dmatrix(g_u,0,g_N-1,0,2);
      free_cvector(g_pushing_phase,0,g_N-1);
      free_cvector(g_pulling_phase,0,g_N-1);
      free_cvector(g_phase,0,g_N-1);
      free_dvector(g_L,0,g_N-1);

      printf("to end, press 0:");
      scanf("%d",&mode);
      if (mode == 0){
        getchar();
      }
      return 0;
    }
