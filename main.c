/*
 * filename: main.c
 * previous version: PPNM_SIM050913_1.cc
 * change from previous version: Refactoring
 * This program simulates centrosome positioning in one-cell embryo
 * Unit meter, kilo-gram, sec
 * to compile: make (see Makefile for detail)
 * Last modified: Mon, 01 Jul 2013 03:05:59 +0900
 */

#include "mtsim.h"
#if defined(_MSC_VER) || defined(__STRICT_ANSI__)
#include "my_getopt.h"
#else
#include <unistd.h>
#endif

// parameters used in functions in this file
static double Buckling_forward_sum;
static double Buckling_backward_sum;
static double Stokes_rad;
static double Visco;
static double Stokes_rotation;
static double Stokes_translation;
static double **DVecNucCen;
static double Vg;
static double k_on;
static double F_dependency;
int N; /* the (maximum) number of MTs per two centrosomes */
int NN;
static double **u;
static unsigned char *pushing_phase;
static unsigned char *pulling_phase;
static unsigned char *phase;
static double Fbuckle[3];
static double Fbackward[6];
static double fjac_pull[6][6];
static int step_counter;
static double BucklingConst;
double *L;
static double *NumberOfMotor;
static unsigned int phase_transition_count;
static void (*usrfun)(double *x, int n, double *fvec, double **fjac);
static double MotorDensity;
static double MotorMaxVel;
static double MotorStallF;
static unsigned char mnewtconverge;
static int mt_start;
static int mt_end;
FILE *f_out8;

void FV_solution(double xx, double *f_v, double *fp_v) {
  *f_v = (k_on*(exp(-1*F_dependency*xx)-1)+Vg)-((xx-Buckling_backward_sum)/Stokes_translation);
  *fp_v = -1*k_on*F_dependency*exp(-1*F_dependency*xx)-1/Stokes_translation;
}

void function_FV3D(double *x, int n, double *fvec, double **fjac) {
  /* x[1],x[2],x[3] are velocity of Nuc for x,y,z axes, respectively */
  int i,j;
  for (i=1;i<=n;i++) {
    fvec[i]=0.0;
    for (j=1;j<=n;j++){
      fjac[i][j]=0.0;
    }
  }
  int mt;
  double A,B,C,D,E;
  for (mt=0; mt<N; mt++) {
    if (pushing_phase[mt]==2) {
      C=u[mt][0]*x[1]+u[mt][1]*x[2]+u[mt][2]*x[3];
      D=C+Vg;
      E=D/k_on;
      A=1-E;
      if (A<=0.0){
        pushing_phase[mt]=1;
        phase_transition_count++;
        TRACE(("A<0 at function_FV3D\n"));
        for (i=1; i<=n; i++) {
          fvec[i] -= u[mt][i-1]*F_dependency*BucklingConst/(L[mt]*L[mt]);
        }
      } else {
        B=log(A);
        for (i=1; i<=n; i++) {
          fvec[i] += u[mt][i-1]*B;
          for (j=1; j<=n; j++) {
            fjac[i][j] += (u[mt][i-1]*u[mt][j-1])/A;
          }
        }
      }
    }
  }
  for (i=1;i<=n;i++){
    fvec[i] = fvec[i]/F_dependency+Fbuckle[i-1]-Stokes_translation*x[i];   /* Fbuckle[0,1,2] */
    for (j=1;j<=n;j++){
      fjac[i][j] = fjac[i][j]/(-1.0*F_dependency*k_on);
      if (j==i) {
        fjac[i][j] = fjac[i][j]-Stokes_translation;}
    }
  }
}

void function_MotorFV (double *x, int n, double *fvec, double **fjac) {
  /* x[1],x[2],x[3] are velocity of Nuc for x,y,z axes, respectively */ 
  int i,j;
  for (i=1;i<=n;i++) {
    fvec[i]=Fbackward[i-1];
    for (j=1;j<=n;j++){
      fjac[i][j]=fjac_pull[i-1][j-1];
    }
  }
  //  for (j=1; j<=6; j++) {
  //    TRACE(("fjac[%d] %lf %lf %lf %lf %lf %lf\n",j,fjac[j][1],fjac[j][2],fjac[j][3],fjac[j][4],fjac[j][5],fjac[j][6]));
  //  }

  int mt;
  int qq;
  double dv, PullingF;
  for (mt=0; mt<N; mt++) {
    if (pulling_phase[mt]==1) {
      if (mt<NN) {qq=0;} else {qq=1;}
      dv=0.0;
      for (j=0; j<3; j++) {
        dv += (x[(j+2)%3+4]*DVecNucCen[qq][(j+1)%3]-x[(j+1)%3+4]*DVecNucCen[qq][(j+2)%3]+x[j+1])*u[mt][j];
      }
      PullingF = MotorStallF*(1-dv/MotorMaxVel);
      for (i=1; i<=n; i++) {
        fvec[i] += NumberOfMotor[mt]*u[mt][i-1]*PullingF;
      }
    }
  }
  for (i=1;i<=(n/2);i++){
    fvec[i] -= Stokes_translation*x[i];   /* Fbackward[0,1,2] */
  }
  for (i=(n/2)+1;i<=n;i++){
    fvec[i] -= Stokes_rotation*x[i];
  }
}

/////////////////// the below function is identical to function_motorFV and needs to be modified ////////////////////
void function_laserMotorFV (double *x, int n, double *fvec, double **fjac) {
  /* x[1],x[2],x[3] are velocity of Nuc for x,y,z axes, respectively */
  int i,j;
  for (i=1;i<=n;i++) {
    fvec[i]=0.0;
    for (j=1;j<=n;j++){
      fjac[i][j]=fjac_pull[i-1][j-1];
    }
  }
  int mt;
  int qq;
  double dv, PullingF;
  /* double dv, PullingF, dFdv; */
  /* double gradient = Stokes_translation; */
  for (mt=mt_start; mt<mt_end; mt++) {
    if (pulling_phase[mt]==1) {
      dv = u[mt][0]*x[1]+u[mt][1]*x[2]+u[mt][2]*x[3];
      if (mt<NN) {qq=0;} else {qq=1;}
      dv=0;
      for (j=0; j<3; j++) {
        dv += (x[(j+2)%3+4]*DVecNucCen[qq][(j+1)%3]-x[(j+1)%3+4]*DVecNucCen[qq][(j+2)%3]+x[j+1])*u[mt][j];
      }
      PullingF = MotorStallF*(1-dv/MotorMaxVel);
      for (i=1; i<=n; i++) {
        fvec[i] += NumberOfMotor[mt]*u[mt][i-1]*PullingF;
      }
    }
  }
  for (i=1;i<=(n/2);i++){
    fvec[i] = fvec[i]+Fbackward[i-1]-Stokes_translation*x[i];   /* Fbackward[0,1,2] */
  }
  for (i=(n/2)+1;i<=n;i++){
    fvec[i] = fvec[i]+Stokes_rotation*x[i];
  }
}

// Newton-Raphson Method for Nonlinear Systems of Equations (ref) Press et al. "Numerical Recipes in C"*/
void mnewt(int ntrial, double x[], int n, double tolx, double tolf) {
  void lubksb(double **a, int n, int *indx, double b[]);
  void ludcmp(double **a, int n, int *indx, double *d);
  int k,i,*indx;
  double errx,errf,d,*fvec,**fjac,*p;

  indx=ivector(1,n);
  p=dvector(1,n);
  fvec=dvector(1,n);
  fjac=dmatrix(1,n,1,n);
  mnewtconverge =0;
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
      return;
    }
    for (i=1; i<=n; i++) p[i] = -fvec[i];
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
      return;
    }
  }
  printf("mnewt: did not converge at %d",k);
  mnewtconverge++;
  free_return(fvec, fjac, n, p, indx);
  return;
}

int main(int argc, char* argv[]) {
  mtGraphics mtg;
  display_setting(&mtg);

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
  FILE *f_out1;
  f_out1 = fopen ("out1.dat","w");
  FILE *f_out2;
  f_out2 = fopen ("out2.dat","w");
  FILE *f_out3;
  f_out3 = fopen ("out3.dat","w");
  FILE *f_out4;
  f_out4 = fopen ("out4.dat","w");
  FILE *f_out10;
  f_out10 = fopen ("out5.dat","w");
  FILE *f_out5;
  f_out5 = fopen ("out_parameter.dat","w");
  FILE *data_for_3D;
  data_for_3D = fopen ("out_for_3D.dat","w");
  FILE *f_out6;
  f_out6 = fopen ("out_MTnumber.dat","w");
  FILE *f_out7;
  f_out7 = fopen ("out_FVcheck.dat","w");
  f_out8 = fopen ("out_3Dcheck.dat","w");

  ////////////////////////////////////////////
  // DECLEARATION of Constants and Variables//
  ////////////////////////////////////////////

  /************CONSTANT NUMBER*****************/
  //PARAMETERS
  unsigned char strain=8; // 0:WT, 1:par-2 w/o LET-99, 2:par-2 with LET-99, 3:par-3 w/o LET-99, 4:par-3 with LET-99, 5: let-99, 6:ric-8, 7:let-99;ric-8, 8:gpr-1/2, 9:PosteriorCortexPullingOnly
  // dynamic instability of MT
  Vg = 0.328e-6; /* Vg: growth velocity of MT [m/sec], standard:0.12 micron/sec */
  double Vs = 0.537e-6; /* Vs: shrinkage velocity of MT [m/sec], standard:0.288 micron/sec*/
  double CatFreq = 0.046; /* fcat: catastrophe frequency of MT [/sec], standard:0.014 */
  double ResFreq = 0.133; /* fres: rescue frequency of MT [/sec], standard:0.044 */
  // drag force of pronucleus
  Visco = 1.0; /* viscosity of cytosol [kg/m sec], standard 0.1->1.0 */
  Stokes_rad = 10.0e-6;
  // pushing forces
  double EI = 10.0e-24; /* rigidity of MT [Nmm]: 4.6 to 41 (Standard 10) */
  double ConstantA = PI*PI;
  BucklingConst = ConstantA * EI;
  k_on = Vg; /* or k_on = Vg/0.8 */
  double F_dependency_single = 3.2e+10; /* 0.034e+12 when velocity decreases to 60% with 15pN (condition resembles in yeast). theoretical value 17.7 (for C. elegans)*/
  // pulling forces
  MotorStallF = 1.1e-12;
  MotorMaxVel = 2.0e-6;
  // length-dependent pulling force
  MotorDensity = 0.050e+6; /* D: density of motor on MT [/m]: 50,000 to 400,000 (standard 100,000)*/
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
  NN = MTTotal[MTDivision+MTDivisionPlus];
  N = 2*NN;
  printf("NN=%d N=%d\n",NN,N);
  // allocate memory
  u=dmatrix(0,N-1,0,5); ///////////////////////////////
  pushing_phase=cvector(0,N-1);
  pulling_phase=cvector(0,N-1);
  phase=cvector(0,N-1);
  L=dvector(0,N-1);
  NumberOfMotor=dvector(0,N-1);
  int CortexPullingMode[N];
  int CortexMotorNumber[N];
  for (k=0;k<N;k++) {CortexPullingMode[k]=0;}
  // arrangement of microtubules: direction of each MT
  double Sdelta[N], Ssita[N];
  Sdelta[0]=0.5*PI; Ssita[0]=0.0;
  kk=0;
  for (k=1;k<NN; k++) {
    if (k > MTTotal[kk]-1){kk++;}
    Sdelta[k]=(PI/2.0)*(1.0-(double)kk/MTDivision);
    Ssita[k]=2.0*PI*k/(MTPerPlane[kk]);
  }
  Sdelta[NN] = -0.5*PI; Ssita[0]=0.0;
  kk=0;
  for (k=NN+1;k<N; k++) {
    if ((k-NN) > MTTotal[kk]-1){kk++;}
    Sdelta[k]=-(PI/2.0)*(1.0-(double)kk/MTDivision);
    Ssita[k]=2.0*PI*(k-NN)/(MTPerPlane[kk]);
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
  double previousL[N];
  double MT[N][3]; /* rectangular coordinates of n-th MT */
  double tempMT[N][3];
  double MTC[3]; /* for calculation of length to cortex used around L448 */
  double min, max, sum; /* to monitor max, mean, min of MT length */
  int currentN; /* current number of active MTs */
  // first and second centrosome//
  DVecNucCen=dmatrix(0,1,0,2);
  double PVecCen[2][3];
  double Nuc[3];  /* position vector of center of the nucleus */
  // forces //
  double ForceC[2][6]; /* force vector on 1st and 2nd centrosome */
  //double PushingForceC[2][3];
  double Bucklingforce_single;
  double pushingF[N];
  // translational and rotational movement of the pronucleus //
  double DirectionDetermination[3];
  double UnitDirection[3];
  double DirectCos[N];
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
  unsigned int eachMT_PTC[N];
  unsigned int pushing_phase_count[5];
  double dv; /* increase of distance between contact point and nucleus devided by time */
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

  //// DECLARATION of Constants and Variables - FINISHED //

  printf("Choose a model: pushing(0), pulling(1), cortical-anchoring(2), pulling_sup1(3), pushing_sup2(4) or length-dependent+cortical-anchoring(5)?:");
  scanf("%d",&mode);

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
        //	  Visco = 0.25;
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
        //	  Visco = 0.5;
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
        //	  Visco = 1.0;
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
        //	  Visco = 2.0;
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
        //	  Visco = 4.0;
        // MTcoefficient = 50; 
        // Fc = 0.25e-12;
        //	  NumberOfMotorOnCortex_PAR3=1.5;
        break;
    }
    Stokes_translation = 6.0*PI*Stokes_rad*Visco;
    Stokes_rotation = -8.0*PI*pow(Stokes_rad,3)*Visco/1.0; /* the formula may be adjusted for easy rotation */
    BucklingConst = ConstantA * EI;

    // OUTPUT parameter LOGs
    fprintf(f_out5,"p=%d\nmode=%d\nMTDivision=%d MT=%d\nVg_um=%5.3lf Vs_um=%5.3lf\nCatFreq=%5.3lf ResFreq=%5.3lf\nFstall_pN=%5.3lf Vmax_um=%5.3f M_mm=%5.3lf\nEI_pNumum=%5.3lf A_um=%5.3lf B_pN=%5.3lf\nH=%5.3lf Stokes_rad=%5.3lf\n\n",p,mode,MTDivision,N,Vg*pow(10,6),Vs*pow(10,6),CatFreq,ResFreq,MotorStallF*pow(10,12),MotorMaxVel*1.0e+6,MotorDensity*pow(10,-3),EI*pow(10,24),k_on*1.0e+6, F_dependency_single*1.0e-12,Visco,Stokes_rad*1.0e+6);

    // color settings
/* #include "color_setting.c" */
    color_setting(&mtg, p);

    /**************** INITIALIZATION *********************/
    // center of the nucleus //
    Nuc[0]=Rad-LL,Nuc[1]=0.0,Nuc[2]=0.0; /**************** start point ********************/
    //Nuc[0]=0.0,Nuc[1]=0.0,Nuc[2]=0.0; /**************** start point ********************/
    DVecNucCen[0][0]=0.0; DVecNucCen[0][1]=0.0; DVecNucCen[0][2]=LL;      
    DVecNucCen[1][0]=0.0; DVecNucCen[1][1]=0.0; DVecNucCen[1][2]=(-1.0)*LL;

    // 1st rotation
    rotationAx[0]=0.0;
    rotationAx[1]=1.0;
    rotationAx[2]=0.0;
    MakeRotationMatrix(RotationMatrix, rotationAx, starting_degree);

    ProductJacVec(VECVEC, RotationMatrix, DVecNucCen[0]);
    for (j=0; j<3; j++) {DVecNucCen[0][j]=VECVEC[j];}
    ProductJacVec(VECVEC, RotationMatrix, DVecNucCen[1]);
    for (j=0; j<3; j++) {DVecNucCen[1][j]=VECVEC[j];}

    for (j=0; j<3; j++) {
      PVecCen[0][j] = Nuc[j]+DVecNucCen[0][j];
      PVecCen[1][j] = Nuc[j]+DVecNucCen[1][j];
    }

    // phase and length of MTs
    for (k=0; k<N; k++){
      L[k]=Vs*dT;  
      random=ran1(idum);
      if (random<(ResFreq/(CatFreq+ResFreq))){phase[k]=1;}
      else {phase[k]=0;}
      // phase[k]=2; /*** add when simulation with increasing number of MTs (Sup Fig S3), start = inactive phase ***/
      u[k][0] = cos(Sdelta[k])*cos(Ssita[k]);
      u[k][1] = cos(Sdelta[k])*sin(Ssita[k]);
      u[k][2] = sin(Sdelta[k]);
      for (j=0;j<3;j++){VECVEC[j]=u[k][j];}
      ProductJacVec(VECVECVEC, RotationMatrix, VECVEC);
      for (j=0;j<3;j++){u[k][j]=VECVECVEC[j];}
      // calculation of  u[k][3], u[k][4], u[k]5]
      if (k<NN) {
        OutProdVector(DVecNucCen[0], VECVECVEC, VECVEC);
      } else {
        OutProdVector(DVecNucCen[1], VECVECVEC, VECVEC);
      }
      for (j=0;j<3;j++){u[k][3+j]=VECVEC[j];}
      pushing_phase[k] = 0;
    }
    currentN = N;
    // currentN = 0;  /*** add when simulation with increasing number of MTs (Sup Fig S3), start = 0 MT ***/

    ////////////////////////////////////////
    ///// Repetition for each time step ////
    ////////////////////////////////////////

    step_counter = -1;
    for (i=0; i<ST; i++)
    {
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
      for (k=0; k<N; k++)
      {
        NumberOfMotor[k]=0.0; /* initialization */
        if (k<NN){qq=0;}else{qq=1;} /* acting on centrosome 1 or 2 */
        pushing_phase[k]=0; /* initial assumption = no pushing */
        switch (phase[k]) /* phase1=growing, phase0=shrinking */
        {
          case 2: /* inactive MTs */
            pushing_phase[k]=0;
            pulling_phase[k]=0;
            break;

          case 1: /* growing phase */
            pulling_phase[k]=1;
            previousL[k]=L[k];
            L[k]+=Vg*dT;
            // if (k > NN-1) {L[k]=Vs*dT;} /****** add when simulation with 1aster (Sup Fig S8) ******/
            for (j=0; j<3; j++) {
              MTC[j] = L[k]*u[k][j];
              tempMT[k][j] = PVecCen[qq][j] + MTC[j];}
            if (sqrt(pow((tempMT[k][0]/Rad),2)+pow((tempMT[k][1]/RadS),2)+pow((tempMT[k][2]/RadS),2)) > 1) { /* the MT reaches the cortex */
              // calculation of the distance from the centrosome to the cell cortex at the angle //
              AA[0]=pow(MTC[0]*RadS,2)+pow(MTC[1]*Rad,2)+pow(MTC[2]*Rad,2); /* quadratic equation */
              AA[1]=2*(MTC[0]*PVecCen[qq][0]*RadS*RadS+MTC[1]*PVecCen[qq][1]*Rad*Rad+MTC[2]*PVecCen[qq][2]*Rad*Rad);
              AA[2]=pow(RadS*PVecCen[qq][0],2)+pow(Rad*PVecCen[qq][1],2)+pow(Rad*PVecCen[qq][2],2)-pow(Rad*RadS,2);
              QuadEqu2(AA,BB);
              L[k] = BB[0]*L[k];

              // THE PUSHING MODEL: estimation of an initial guess to be used in Newton-Raphson method 
              if (mode == 0) {
                //  if (qq=0) { /************add if 1aster ******/ 
                pushing_phase[k]=1; /* temporary assume buckling */
                for (j=0; j<3; j++) {
                  tempMT[k][j] = PVecCen[qq][j] + L[k]*u[k][j]; /* coordinate of contact point */
                  DirectionDetermination[j] -= u[k][j]/(L[k]*L[k]); /* to estimate the direction vector of an initial guess to be used in the Newton-Raphson method */
                }
              } else {
                pushing_phase[k]=0;
              }

              // THE CORTICAL-ANCHORING MODEL (Sup Fig S6) 
              if (mode == 2) {
                for (j=0; j<6; j++){
                  ForceC[qq][j] += u[k][j] * Fc;} /* To examine the cortical-anchoring model with the additional assumptions (Sup Fig S6), conditions described in the figure legend are added here */
              }

              // PUSHING MODELS WITH DIFFERENT FORMULAS (Sup Fig S2)	
              if (mode == 4) {
                for (j=0; j<6; j++){
                  ForceC[qq][j] -= u[k][j] * 1.0e-22 * pow(L[k],-2);}  /* To examine the pushing forces expressed with different formulas (Sup Fig S2), the formulas shown here are modified as described in the figure legend */
              }

              // LENGTH-DEPENDENT PULLING + CORTEX-PULLING
              if (mode == 5) {
                if (CortexPullingMode[k]==0) {
                  for (j=0; j<3; j++) {
                    tempMT[k][j] = PVecCen[qq][j] + L[k]*u[k][j];
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
                  NumberOfMotor[k] += CortexMotorNumber[k];
                  CortexPullingMode[k]--;
                } 
              }

              // if (ran1(idum) < CatFreq*20*dT){phase[k] = 0;} /* if assuming CatFreq increases upon contact with the cortex */
              } else {
                pushing_phase[k]=0; // To ensure that pushing_phase of the MT without contacting the cortex is 0 // 
              }
              // switching from growth phase to shrinkage phase // 
              if (ran1(idum) < CatFreq*dT){phase[k] = 0;}
              break;

              case 0: /* shrinkage phase */
              if (L[k] > Rad/25) {L[k] -= Vs*dT;pulling_phase[k]=1;}  /* the lengh of MT exceeds 1micron */
              //  if (k > NN-1) {L[k]=Vs*dT;} /************* add if 1aster (Sup Fig. S8)********************/

              // switching from shrinkage phase to growth phase // 
              if (ran1(idum) < ResFreq*dT){phase[k] = 1;}
              pushing_phase[k]=0;
              break;
            }

            // PULLING MODELS WITH DIFFERENT FORMULAS (Sup Fig S1)
            if (mode == 3) { 
              for (j=0;j<6;j++){
                ForceC[qq][j] += u[k][j] * 1.0e-7 * pow(L[k],1);} /* pulling models with different fomulas (Sup Fig S1), the formulas below are modified as described in the figure legend */  
            }

            ///////////////////////////////////////////
            // LENGTH-DEPENDENT PULLING + CORTEX PULLING
            if ((mode == 1)||(mode == 5)) {
              random = ran1(idum);
              np = L[k]*MotorDensity;
              cumpoisson = 0.0;
              motornumber = -1;
              do {
                motornumber++;
                cumpoisson += Poisson(np,motornumber);
              } while (random > cumpoisson);
              NumberOfMotor[k] += motornumber;
              //TRACE(("%4d %4d N=%3.1lf MT=(%3.2lf %3.2lf %3.2lf %3.2lf %3.2lf %3.2lf)\n",i,k,NumberOfMotor[k],u[k][0],u[k][1],u[k][2],u[k][3],u[k][4],u[k][5]));
              for (j=0;j<6;j++){
                ForceC[qq][j] += NumberOfMotor[k] * u[k][j] * MotorStallF;} /* an initial guess */
            }

            // to monitor microtubules profile
            sum += L[k];
            if (max < L[k]){max = L[k];}
            if (min > L[k]){min = L[k];}
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
          usrfun = function_FV3D;
          if (i%100==0) {fprintf(f_out7,"%d ",i);}
          UnitVector(DirectionDetermination,UnitDirection); /* unit vector of an initial guess of pronuclear migration */
          Buckling_forward_sum=0.0;
          Buckling_backward_sum=0.0;
          StallSum=0.0;
          for (k=0; k<N; k++) {
            if (pushing_phase[k]!=0) {
              DirectCos[k]=InnProdVector(u[k],UnitDirection);
              if (DirectCos[k]<0) {
                Buckling_forward_sum -= DirectCos[k]*BucklingConst/(L[k]*L[k]);
                StallSum -= DirectCos[k];
                pushing_phase[k]=2; /* temporal assumption that MTs directed toward the opposite direction of pronuclear migration and pushing the cortex do not buckle */
              } else {
                Buckling_backward_sum += DirectCos[k]*BucklingConst/(L[k]*L[k]);
                pushing_phase[k]=1; /* temporal assumption that MTs directed toward the direction of pronuclear migration and pushing the cortex buckle  */
              }
            }
          }
          if (Buckling_forward_sum < Buckling_backward_sum) {
            printf("something is wrong: cannot determine buckling direction!\n");
          } else {
            F_dependency=F_dependency_single/StallSum;
            Vnuc_buckle=Stokes_function(Buckling_forward_sum,Stokes_translation,Buckling_backward_sum);
            Vmt_buckle=FV_function(Buckling_forward_sum,Vg,k_on,F_dependency);
            if (Vnuc_buckle<=Vmt_buckle){  /* buckling is dominant */
              if (i%100==0) {fprintf(f_out7,"BK ");}
              F_soln = Buckling_forward_sum;
              Vnuc_soln = Vnuc_buckle;
              Vmt_soln = Vmt_buckle;
              for (k=0;k<N;k++) {
                if (pushing_phase[k]!=0) {
                  pushing_phase[k]=1;} /* temporal assumption that all MTs buckle */
              }
            } else {  /* FV is dominant */
              if (i%100==0) {
                fprintf(f_out7,"FV ");
              }
              F_soln = rtsafe(FV_solution,0.0,Buckling_forward_sum,xacc);
              Vnuc_soln = Stokes_function(F_soln,Stokes_translation,Buckling_backward_sum);
              Vmt_soln = FV_function(F_soln,Vg,k_on,F_dependency);
            }
          }
          for (j=1; j<=3; j++){
            tempNucVel[j]=Vnuc_soln*UnitDirection[j-1]; /* an initial guess of the velocity of the pronucleus to be used in Newton-Raphson method*/
            if (i%100==0) {fprintf(f_out7,"%5.4f ",tempNucVel[j]*pow(10,6));}
          }

          F_dependency=F_dependency_single;
          for (j=0;j<3;j++){
            Fbuckle[j]=0.0;  /* Fbuckle[0,1,2]*/
          }
          for (k=0; k<N; k++) {
            if (pushing_phase[k]==1){
              Bucklingforce_single=BucklingConst/(L[k]*L[k]);	  
              for (j=0;j<3;j++){
                Fbuckle[j] -= u[k][j]*Bucklingforce_single;
              }
            }
          }
          usrfun = function_FV3D;
          mnewt(10,tempNucVel,3,tolx,tolf); /* Newton-Raphson method to revise the initial guess of the velocity of the pronucleus */

          // THE PUSHING MODEL-2: solve the set of equation using Newton-Raphson method with the revised initial guess
          cycle_count=0;
          final_cycle_check=0;
          for (k=0; k<N; k++) {eachMT_PTC[k]=0;}
          do {
            phase_transition_count=0;
            for (j=0;j<3;j++){
              Fbuckle[j]=0.0;  /* Fbuckle[0,1,2]*/
            }
            for (k=0; k<N; k++) {
              if (pushing_phase[k]!=0){
                Bucklingforce_single=BucklingConst/(L[k]*L[k]);	  
                Vmt_buckle=FV_function(Bucklingforce_single,Vg,k_on,F_dependency_single);
                dv = -1.0 * (u[k][0]*tempNucVel[1]+u[k][1]*tempNucVel[2]+u[k][2]*tempNucVel[3]);
                if (eachMT_PTC[k]>10) {  /* to avoid examining same condition repeatedly */
                  if ((dv > Vg) || (dv < Vmt_buckle)||(pushing_phase[k]!=2)) { /* check the length. If not appropriate, do not exit */
                    phase_transition_count++;eachMT_PTC[k]++;
                  }
                  pushing_phase[k]=2;
                  if (eachMT_PTC[k]>25){eachMT_PTC[k]=0;} /* if phase 2 is not appropriate, go back to usual classification */
                } else {
                  if (dv > Vg){
                    if (pushing_phase[k]!=3){
                      phase_transition_count++;
                      eachMT_PTC[k]++;
                    }
                    if (pushing_phase[k]==1){
                      pushing_phase[k]=2;
                    } else {
                      pushing_phase[k]=3;
                    }
                  } else {
                    if (dv > Vmt_buckle) {  /* FV is dominant */		
                      if (pushing_phase[k]!=2){phase_transition_count++;eachMT_PTC[k]++;
                      }
                      pushing_phase[k]=2;
                    } else {
                      if (pushing_phase[k]!=1){
                        phase_transition_count++;eachMT_PTC[k]++;
                      }
                      if (pushing_phase[k]==3) {
                        pushing_phase[k]=2;
                      } else {
                        pushing_phase[k]=1;
                        for (j=0;j<3;j++){
                          Fbuckle[j] -= u[k][j]*Bucklingforce_single;
                        }
                      }
                    }    
                  }
                }
              }		
            }	      
            if ((phase_transition_count!=0)||(cycle_count==0)) {
              mnewt(10,tempNucVel,3,tolx,tolf);  /* Newton-Raphson method */
            }
            if (i%100==0) {fprintf(f_out8,"%d %d %d\n",i,cycle_count,phase_transition_count);}
            cycle_count++;
          } while ((mnewtconverge!=0)||(cycle_count<=1)||((phase_transition_count!=0)&&(cycle_count<1000))); /* repeat until the solution satisfies all equations and conditions */

          if (phase_transition_count!=0){ /* in case the solution is not obtained within 1000 cycles*/
            printf("exit without convergence at t=%d PTC=%d\n",i,phase_transition_count);
            fprintf(f_out7,"exit without convergence at t=%d PTC=%d\n",i,phase_transition_count);
            printf("Cen1=(%3.1f,%3.1f,%3.1f) Cen2=(%3.1f,%3.1f,%3.1f)\n",PVecCen[0][0]*1.0e+6,PVecCen[0][1]*1.0e+6,PVecCen[0][2]*1.0e+6,PVecCen[1][0]*1.0e+6,PVecCen[1][1]*1.0e+6, PVecCen[1][2]*1.0e+6);
            for (k=0; k<N; k++){
              if (pushing_phase[k]!=0){
                printf("%d %d L=%3.1f (%3.2f,%3.2f,%3.2f)\n",k,pushing_phase[k],L[k]*1000000,u[k][0],u[k][1],u[k][2]);
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
              ForceC[qq][j]=0.0;}
            pushing_phase_count[j+1]=0;
          }
          for (k=0; k<N; k++) {
            if (k<NN) {qq=0;} else {qq=1;}
            if (pushing_phase[k]==3) {
              // L[k]=previousL[k]+Vg*dT; /* this line can be omitted because of redundancy*/
              pushingF[k]=0.0;
              pushing_phase_count[3]++;
            } 
            if (pushing_phase[k]==2) {
              dv = -1 * (u[k][0]*tempNucVel[1]+u[k][1]*tempNucVel[2]+u[k][2]*tempNucVel[3]);
              // L[k]=L[k]+dv*dT; /* this line can be omitted because of redundancy*/
              pushingF[k] = -1*log(1-(Vg-dv)/k_on)/F_dependency_single;
              for (j=0;j<3;j++) {
                ForceC[qq][j] -= pushingF[k]*u[k][j];}
              pushing_phase_count[2]++;
            }
            if (pushing_phase[k]==1) {
              pushingF[k]=BucklingConst/(L[k]*L[k]);
              // L[k]=previousL[k]+dT*FV_function(pushingF[k],Vg,k_on,F_dependency_single); /* omit this line assuming buckling MTs do not elongate any more*/
              for (j=0;j<3;j++) {
                ForceC[qq][j] -= pushingF[k]*u[k][j];}
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
                tempNucVel[j]=(ForceC[0][j-1]+ForceC[1][j-1])/Stokes_translation; 
              } else {
                tempNucVel[j]=(ForceC[0][j-1]+ForceC[1][j-1])/Stokes_rotation;
              }
            }
            usrfun = function_MotorFV;
            cycle_count=0;
            TRACE(("%4d init ",i));
            for (j=1; j<=6; j++) {
              if (j<=3) {
                TRACE(("%4.3lf ",tempNucVel[j]*1.0e+6));
              } else {
                TRACE(("%4.3lf ",tempNucVel[j]*100));
              }
            }
            TRACE(("\n"));
            for (k=0; k<N; k++) {eachMT_PTC[k]=0;}
            do {
              phase_transition_count=0;
              for (j=0;j<6;j++){
                Fbackward[j]=0.0;
                for (jj=0; jj<6; jj++) {
                  fjac_pull[j][jj] = 0.0;
                }}
              for (k=0; k<N; k++) {
                if (pulling_phase[k]!=0){
                  if (k<NN) {qq=0;} else {qq=1;}
                  dv=0.0;
                  for (j=0; j<3; j++) {
                    dv += (tempNucVel[(j+2)%3+4]*DVecNucCen[qq][(j+1)%3]-tempNucVel[(j+1)%3+4]*DVecNucCen[qq][(j+2)%3]+tempNucVel[j+1])*u[k][j];
                  }
                  if (dv > MotorMaxVel){ /* the motors on this MT do not exert forces: pulling_phase[k]=3 */
                    if (pulling_phase[k]!=3) {
                      phase_transition_count++;
                      eachMT_PTC[k]++;
                    }
                    pulling_phase[k]=3;
                  } else {
                    if (dv < 0) { /* the motors on this MT exert the maximum (stall) force: pulling_phase[k]=2 */
                      if (pulling_phase[k]!=2) {
                        phase_transition_count++;
                        eachMT_PTC[k]++;
                      }
                      pulling_phase[k]=2;
                      for (j=0;j<6;j++) {
                        Fbackward[j] += u[k][j] * MotorStallF * NumberOfMotor[k];
                      }
                    } else {
                      if (pulling_phase[k]!=1) { /* the motors on this MT exert force dependent on their velocity: pulling_phase[k]=1 */
                        phase_transition_count++;
                        eachMT_PTC[k]++;
                      }
                      pulling_phase[k]=1;
                      for (j=0; j<3; j++) {dvid[j]=u[k][j];}
                      dvid[3]=DVecNucCen[qq][2]*u[k][1]-DVecNucCen[qq][1]*u[k][2];
                      dvid[4]=DVecNucCen[qq][0]*u[k][2]-DVecNucCen[qq][2]*u[k][0];
                      dvid[5]=DVecNucCen[qq][1]*u[k][0]-DVecNucCen[qq][0]*u[k][1];
                      for (j=0; j<6; j++) {
                        for (jj=0; jj<6; jj++) {
                          fjac_pull[j][jj] -= (MotorStallF*NumberOfMotor[k]/MotorMaxVel)*dvid[jj]*u[k][j];
                        }
                      }
                    }
                  }
                }
              }
              for (j=0; j<3; j++) {fjac_pull[j][j] -= Stokes_translation;}
              for (j=3; j<6; j++) {fjac_pull[j][j] -= Stokes_rotation;}
              if ((phase_transition_count!=0)||(cycle_count==0)) {
                mnewt(10,tempNucVel,6,tolx,tolf);
              }
              if (i%100==0) {fprintf(f_out8,"%d %d %d\n",i,cycle_count,phase_transition_count);}
              cycle_count++;
            } while ((cycle_count<=1)||((phase_transition_count!=0)&&(cycle_count<1000))); /* repeat until the solution satisfies all equations and conditions */

            TRACE(("%4d %4d ",i,cycle_count));
            for (j=1; j<=6; j++) {
              if (j<=3) {
                TRACE(("%4.3lf ",tempNucVel[j]*1.0e+6));
              } else {
                TRACE(("%4.3lf ",tempNucVel[j]*100));
              }
            }
            TRACE(("\n"));

            if (phase_transition_count!=0){ /* when the solution is not obtained within 1000 cycles*/
              printf("exit without convergence at t=%d PTC=%d\n",i,phase_transition_count);
              fprintf(f_out7,"exit without convergence at t=%d PTC=%d\n",i,phase_transition_count);
              printf("Cen1=(%3.1f,%3.1f,%3.1f) Cen2=(%3.1f,%3.1f,%3.1f)\n",PVecCen[0][0]*1.0e+6,PVecCen[0][1]*1.0e+6,PVecCen[0][2]*1.0e+6,PVecCen[1][0]*1.0e+6,PVecCen[1][1]*1.0e+6, PVecCen[1][2]*1.0e+6);
              break;
            }
          } else {  /******** laser ablation *****************/
            usrfun = function_laserMotorFV;
            for (laser_centrosome = 0; laser_centrosome<2; laser_centrosome++){
              if (laser_centrosome==0) {
                mt_start=0;
                mt_end=NN;
              } else {
                mt_start=NN;
                mt_end=N;
              }
              for (j=1;j<=3;j++){
                tempNucVel[j]=(ForceC[laser_centrosome][j-1])/Stokes_translation; /* forces do not include dT **/
              }
              cycle_count=0;
              for (k=mt_start; k<mt_end; k++) {eachMT_PTC[k]=0;}
              do {
                phase_transition_count=0;
                for (j=0;j<3;j++){
                  Fbackward[j]=0.0;  /* Fbackward[0,1,2]*/
                }
                for (k=mt_start; k<mt_end; k++) {
                  if (pulling_phase[k]!=0){
                    dv = (u[k][0]*tempNucVel[1]+u[k][1]*tempNucVel[2]+u[k][2]*tempNucVel[3]);
                    if (dv > MotorMaxVel){ /* the motors on this MT do not exert forces */
                      if (pulling_phase[k]!=3) {
                        phase_transition_count++;
                        eachMT_PTC[k]++;
                      }
                      pulling_phase[k]=3;
                    } else {
                      if (dv < 0) { /* the motors on this MT exert the maximum (stall) force */
                        if (pulling_phase[k]!=2) {
                          phase_transition_count++;
                          eachMT_PTC[k]++;
                        }
                        pulling_phase[k]=2;
                        for (j=0;j<3;j++) {
                          Fbackward[j] += u[k][j] * MotorStallF * NumberOfMotor[k];
                        }
                      } else {
                        if (pulling_phase[k]!=1) { /* the motors on this MT exert force dependent on their velocity */
                          phase_transition_count++;
                          eachMT_PTC[k]++;
                        }
                        pulling_phase[k]=1;
                      }
                    }
                  }
                }
                if ((phase_transition_count!=0)||(cycle_count==0)) {
                  mnewt(10,tempNucVel,3,tolx,tolf);
                }
                if (i%100==0) {fprintf(f_out8,"%d %d %d\n",i,cycle_count,phase_transition_count);}
                cycle_count++;
              } while ((cycle_count<=1)||((phase_transition_count!=0)&&(cycle_count<1000))); /* repeat until the solution satisfies all equations and conditions */

              if (phase_transition_count!=0){ /* when the solution is not obtained within 1000 cycles*/
                printf("exit without convergence at t=%d PTC=%d\n",i,phase_transition_count);
                fprintf(f_out7,"exit without convergence at t=%d PTC=%d\n",i,phase_transition_count);
                //  printf("Cen1=(%3.1f,%3.1f,%3.1f) Cen2=(%3.1f,%3.1f,%3.1f)\n",Cen[0][0]*1.0e+6,Cen[0][1]*1.0e+6,Cen[0][2]*1.0e+6,Cen[1][0]*1.0e+6,Cen[1][1]*1.0e+6, Cen[1][2]*1.0e+6);
                break;
              }
              for (j=1; j<=3; j++){
                CenVel[laser_centrosome][j]=tempNucVel[j];
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
          //		pullingF[k] = MotorStallF;
          //	for (j=0;j<6;j++) {
          //		  ForceC[qq][j] += pullingF[k]*NumberOfMotor[k]*u[k][j];}
          //		pulling_phase_count[2]++;
          //	      }
          //	      if (pulling_phase[k]==1) {
          //		dv = (u[k][0]*tempNucVel[1]+u[k][1]*tempNucVel[2]+u[k][2]*tempNucVel[3]); ////////////////////CHANGE/////
          //		pullingF[k]=MotorStallF*(1-dv/MotorMaxVel);
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
              DeltaNuc[j]=Nuc[j]-OldNuc[j];
              OldNucPosition[j]=DeltaNucStandard[j]-OldNuc[j];
            }
            OldNucLength=sqrt(OldNucPosition[0]*OldNucPosition[0]+OldNucPosition[2]*OldNucPosition[2]);
            if ((OldNucLength>Rad*0.24)&&(OldNucLength<Rad*0.72)) {
              UnitOldNucPosition[0]=OldNucPosition[0]/OldNucLength;
              UnitOldNucPosition[2]=OldNucPosition[2]/OldNucLength;
              DeltaNucAxis1=DeltaNuc[0]*UnitOldNucPosition[0]+DeltaNuc[2]*UnitOldNucPosition[2];
              DeltaNucAxis2=DeltaNuc[0]*UnitOldNucPosition[2]-DeltaNuc[2]*UnitOldNucPosition[0];
              //draw graph in window 9
              XFillRectangle(mtg.d,mtg.w9,mtg.gc9,(int)(WIN_WIDTH/2+500*1000000*DeltaNucAxis2/(dT*UNITTIME))-1,(int)(WIN_HEIGHT/2-500*1000000*DeltaNucAxis1/(dT*UNITTIME))-1,2,2);
            }
          }
          for (j=0;j<3;j++) {
            OldNuc[j]=Nuc[j];
          }
        }

        // TRANSLATIONAL MOVEMENT OF THE PRONUCLEUS
        CalculateVel=0.0;
        for (j=1; j<=3; j++) {
          if (((mode==1)||(mode==5))&&(i>=laserST)) {
            tempNucVel[j]=(CenVel[0][j]+CenVel[1][j])/2;
            PVecCen[0][j-1] += CenVel[0][j]*dT;
            PVecCen[1][j-1] += CenVel[1][j]*dT;
          } else {
            if ((mode!=0)&&(mode!=1)&&(mode!=5)){ /* models other than THE PUSHING and PULLING MODELS */
              tempNucVel[j]=(ForceC[0][j-1]+ForceC[1][j-1])/Stokes_translation;
            }
          }
          CalculateVel+=tempNucVel[j]*tempNucVel[j];
          Nuc[j-1] += tempNucVel[j]*dT;
        }

        // the pronucleus does not cross over the cell cortex 
        NuclearDistanceRatio = sqrt(pow((Nuc[0]/(Rad-LL)),2)+pow((Nuc[1]/(RadS-LL)),2)+pow((Nuc[2]/(RadS-LL)),2));
        if (NuclearDistanceRatio > 1) { /* when the pronucleus contacts the cell cortex */
          if (mode==0) {TRACE(("nucleus crosses over the cortex!!\n"));}
          else {
            AAA[0]=0.0;
            AAA[1]=0.0;
            AAA[2]=0.0;
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
            AAA[2]-=1;
            TRACE(("QuadEqu2 at L1276\n"));
            QuadEqu2(AAA,BBB);
            for (j=0; j<3; j++){
              Nuc[j] += BBB[0]*tempNucVel[j+1]*dT;
            }
            KKK=0;
            for (j=0; j<1; j++){			  
              KKK+= (-1)*Nuc[j]*(1-BBB[0])*(tempNucVel[j+1]*dT)/((Rad-LL)*(Rad-LL));
            }
            for (j=1; j<3; j++){			  
              KKK+= (-1)*Nuc[j]*(1-BBB[0])*(tempNucVel[j+1]*dT)/((RadS-LL)*(RadS-LL));
            }
            KKK=KKK/(pow(Nuc[0],2)/pow((Rad-LL),4)+pow(Nuc[1],2)/pow((RadS-LL),4)+pow(Nuc[2],2)/pow((RadS-LL),4));
            for (j=0; j<1; j++){
              Nuc[j]+=(1-BBB[0])*tempNucVel[j+1]*dT+KKK*Nuc[j]/pow((Rad-LL),2);}
            for (j=1; j<3; j++){
              Nuc[j]+=(1-BBB[0])*tempNucVel[j+1]*dT+KKK*Nuc[j]/pow((RadS-LL),2);}
          }
        }

        // calculation of velocity of the pronucleus
        if(i%UNITTIME==0){
          NewDis = sqrt((Rad-Nuc[0])*(Rad-Nuc[0])+Nuc[1]*Nuc[1]+Nuc[2]*Nuc[2]);
          NewVel = (NewDis-OldDis)/(dT*UNITTIME);
        }
        DistanceFromPP=sqrt((Rad-Nuc[0])*(Rad-Nuc[0])+Nuc[1]*Nuc[1]+Nuc[2]*Nuc[2]);

        if (((mode==1)||(mode==5))&&(i<laserST)) {
          // ROTATIONAL MOVEMENT OF THE PRONUCLEUS
          Rotation[0]=tempNucVel[4];
          Rotation[1]=tempNucVel[5];
          Rotation[2]=tempNucVel[6];
          MakeRotationMatrix(RotationMatrix, Rotation, dT);
          TRACE(("RotationMatrix\n"));
          for (j=0; j<3; j++) {
            for (jj=0; jj<3; jj++) {
              TRACE(("%4.3lf ",RotationMatrix[j][jj]));
            }
            TRACE(("\n"));
          }
          ProductJacVec(VECVEC, RotationMatrix, DVecNucCen[0]);
          for (j=0; j<3; j++) {DVecNucCen[0][j]=VECVEC[j];}
          ProductJacVec(VECVEC, RotationMatrix, DVecNucCen[1]);
          for (j=0; j<3; j++) {DVecNucCen[1][j]=VECVEC[j];}
          for (j=0; j<3; j++) {
            PVecCen[0][j] = Nuc[j]+DVecNucCen[0][j];
            PVecCen[1][j] = Nuc[j]+DVecNucCen[1][j];
          }

          // phase and length of MTs
          for (k=0; k<N; k++){
            if (k<NN) {qq=0;} else {qq=1;}
            for (j=0;j<3;j++){VECVEC[j]=u[k][j];}
            ProductJacVec(VECVECVEC, RotationMatrix, VECVEC);
            for (j=0;j<3;j++){
              u[k][j]=VECVECVEC[j];
              MT[k][j] = PVecCen[qq][j] + L[k]*u[k][j];
            }
            // calculation of  u[k][3], u[k][4], u[k]5]
            OutProdVector(DVecNucCen[qq], VECVECVEC, VECVEC);
            for (j=0;j<3;j++){u[k][3+j]=VECVEC[j];}
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
            for (k=0; k<N; k++){
              if (k<NN) {qq=0;} else {qq=1;}
              for (j=0; j<3; j++) {
                MT[k][j] = PVecCen[qq][j] + L[k]*u[k][j];
              }	
            }	    
          } else {
            if (t == 0) {
              TRACE(("L1373 no rotation!\n"));
            }
            for (j=0; j<3; j++) {
              PVecCen[0][j] = Nuc[j]+DVecNucCen[0][j];
              PVecCen[1][j] = Nuc[j]+DVecNucCen[1][j];
            }
            for (k=0; k<N; k++){
              if (k<NN) {qq=0;} else {qq=1;}
              for (j=0; j<3; j++) {
                MT[k][j] = PVecCen[qq][j] + L[k]*u[k][j];
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
        save_logs(i, p, N, f_out1, f_out2, f_out3, f_out4, f_out10, data_for_3D, Nuc, PVecCen, MT);
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

      free_dmatrix(u,0,N-1,0,2);
      free_cvector(pushing_phase,0,N-1);
      free_cvector(pulling_phase,0,N-1);
      free_cvector(phase,0,N-1);
      free_dvector(L,0,N-1);

      printf("to end, press 0:");
      scanf("%d",&mode);
      if (mode ==0){
        getchar();
      }
      return 0;
    }
