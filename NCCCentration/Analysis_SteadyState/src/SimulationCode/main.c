/*
 * filename: main.c
 * previous version: PPNM_SIM050913_1.cc
 * change from previous version: Refactoring
 * This program simulates centrosome positioning in one-cell embryo
 * Unit meter, kilo-gram, sec
 * to compile: make (see Makefile for detail)
 * Author: Akatsuki Kimura <akkimura@nig.ac.jp>
 *         Akira Funahashi <funa@bio.keio.ac.jp>
 * Last modified: Thu, 25 Jul 2013 03:16:52 +0900
 */

#include "mtsim.h"
#if defined(_MSC_VER) || defined(__STRICT_ANSI__)
#include "my_getopt.h"
#else
#include <unistd.h>
#endif

int main(int argc, char* argv[]) {
  mtGlobal g;      /* parameters which was global variable in previous code */
  mtGraphics mtg;
  display_setting(&mtg);

  /*  Variables for getopt() */
  int ch;
  extern char *optarg;
  extern int optind;
  char* myname;
  boolean is_check = false;
  boolean is_verbose = false;
  unsigned int strain=8; // 0:WT, 1:par-2 w/o LET-99, 2:par-2 with LET-99, 3:par-3 w/o LET-99, 4:par-3 with LET-99, 5: let-99, 6:ric-8, 7:let-99;ric-8, 8:gpr-1/2, 9:PosteriorCortexPullingOnly

  /* Parse options */
  myname = argv[0];
  while ((ch = getopt(argc, argv, "s:cvh")) != -1){
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
      case 's':
        strain = atoi(optarg);
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

  if(strain > 9){
    usage(myname);
    exit(1);
  }

  //file handling
  FILE *f_out1 = fopen("out1.dat","w");
  FILE *f_out2 = fopen("out2.dat","w");
  FILE *f_out3 = fopen("out3.dat","w");
  FILE *f_out4 = fopen("out4.dat","w");
  FILE *f_out5 = fopen("out5.dat","w");
  /* FILE *f_out_param = fopen("out_parameter.dat","w"); */
  /* FILE *f_out_mtnum = fopen("out_MTnumber.dat","w"); */
  /* FILE *f_out_for3d = fopen("out_for_3D.dat","w"); */
  FILE *f_out_fvcheck = fopen("out_FVcheck.dat","w");
  FILE *f_out_3dcheck = fopen("out_3Dcheck.dat","w");
  FILE *f_out_vector;
  char vector_filename[256];
  sprintf(vector_filename, "./vectordata/centering_vector_strain%d_rad%.0lf_rads%.0lf.dat",strain,Rad*pow(10,6),RadS*pow(10,6));
  f_out_vector = fopen(vector_filename, "w");

  ////////////////////////////////////////////
  // DECLEARATION of Constants and Variables//
  ////////////////////////////////////////////

  /************CONSTANT NUMBER*****************/
  //PARAMETERS
  // drag force of pronucleus
  g.Visco = 1.0; /* viscosity of cytosol [kg/m sec], standard 0.1->1.0 */
  g.Stokes_rad = 10.0e-6;
  // pulling forces
  g.MotorStallF = 1.1e-12;
  g.MotorMaxVel = 2.0e-6;
  // length-dependent pulling force
  g.MotorDensity = 0.10e+6; /* D: density of motor on MT [/m]: 50,000 to 400,000 (standard 100,000)*/
  // corical pulling force
  int CortexPullingDuration = 1; /* [timepoint] */
  double CortexPullingFreq_PAR3 = 0.25; /* the probability of a dynein being in its active state */
  double CortexPullingFreq_PAR2 = 1.0;
  double CortexPullingFreq_LET99 = 0.1;
  double CortexPullingFreq;

  // arrangement of microtubules: calculation of the (maximum) number of MTs
  double starting_degree; /* starting degree of aster: for the simulation with single aster (Sup Fig S9) this value should be changed*/
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
  g.NN = MTTotal[MTDivision+MTDivisionPlus];
  g.N = 2*g.NN;
  printf("NN=%d N=%d\n", g.NN, g.N);
  // allocate memory
  g.u = dmatrix(0,g.N-1,0,5); ///////////////////////////////
  g.pushing_phase = cvector(0,g.N-1);
  g.pulling_phase = cvector(0,g.N-1);
  g.phase = cvector(0,g.N-1);
  g.L = dvector(0,g.N-1);
  g.NumberOfMotor = dvector(0,g.N-1);
  int CortexPullingMode[g.N];
  g.CortexMotorNumber = dvector(0,g.N-1);
  for (k=0;k<g.N;k++) {CortexPullingMode[k]=0;g.CortexMotorNumber[k]=0;}
  // arrangement of microtubules: direction of each MT
  double Sdelta[g.N], Ssita[g.N];
  Sdelta[0]=0.5*PI; Ssita[0]=0.0;
  kk=0;
  for (k=1;k<g.NN; k++) {
    if (k > MTTotal[kk]-1){kk++;}
    Sdelta[k]=(PI/2.0)*(1.0-(double)kk/MTDivision);
    Ssita[k]=2.0*PI*k/(MTPerPlane[kk]);
  }
  Sdelta[g.NN] = -0.5*PI; Ssita[0]=0.0;
  kk=0;
  for (k=g.NN+1;k<g.N; k++) {
    if ((k-g.NN) > MTTotal[kk]-1){kk++;}
    Sdelta[k]=-(PI/2.0)*(1.0-(double)kk/MTDivision);
    Ssita[k]=2.0*PI*(k-g.NN)/(MTPerPlane[kk]);
  }
  // int MTcoefficient = 100; /* a constant used in simulations with increasing number of MTs (Sup Fig S5) */

  g.CortexMotorDensity = 16*pow(MTDivision,2)/(pow(PI,3)*pow(Rad,2)); /* [/m^2] */
  /* g.CortexMotorDensity = pow(MTDivision,2)/(pow(PI,3)*pow(Rad,2))/2; /\* [/m^2] *\/ */

  /************VARIABLES**********************/
  int i; /* steps */
  /* double t=0.0; /\* time [sec] *\/ */
  int j,jj; /* x,y,z axes */
  int qq; /* centrosome 0 or 1 */
  double VECVEC[3], VECVECVEC[3];
  // microtubules //
  double previousL[g.N];
  /* double MT[g.N][3]; /\* rectangular coordinates of n-th MT *\/ */
  double tempMT[g.N][3];
  double MTC[3]; /* for calculation of length to cortex used around L448 */
  double min, max, sum; /* to monitor max, mean, min of MT length */
  int currentN; /* current number of active MTs */
  // first and second centrosome//
  g.DVecNucCen = dmatrix(0,1,0,2);
  double PVecCen[2][3];
  double Nuc[3];  /* position vector of center of the nucleus */
  // forces //
  double ForceC[2][6]; /* force vector on 1st and 2nd centrosome */
  // translational and rotational movement of the pronucleus //
  double DirectionDetermination[3];
  double *tempNucVel; 
  tempNucVel=dvector(1,6);
  double CalculateVel;
  unsigned int cycle_count;
  unsigned int eachMT_PTC[g.N];
  double dv; /* increase of distance between contact point and nucleus devided by time */
  int step_counter;
  ////////////////////////ROTATION
  /* double Rotation[3]; /\* rotational vector *\/ */
  double rotationAx[3]; /* rotational axis */
  ////////////////////// Nucleus does not cross over the cortex /////////////////
  /* double KKK; */
  double NuclearDistanceRatio;
  double OldDis, NewDis, OldVel, NewVel;
  double DistanceFromPP;
  // solution of a quadratic equation //
  double AA[3]; /* three coefficients of the equation */
  double BB[2]; /* two solutions of the equation */
  // variables for local movement
  double DeltaNucStandard[3]={0.0,0.0,0.0};
  double DeltaNuc[3], OldNuc[3], OldNucLength, OldNucPosition[3], UnitOldNucPosition[3], DeltaNucAxis1, DeltaNucAxis2;
  double dvid[6];
  double RotationMatrix[3][3];

  char filefin[40];
  char filefin2[40];

  // thresholds in Newton-Raphson method
  double tolx = 1.0e-12; /* 1e-6 um/sec */
  double tolf = 1.0e-17; /* 1e-5 pN */
  boolean did_converge = false;
  void (*usr_func)(double*, int, double*, double**, mtGlobal*); /* for callback function */

  //// DECLARATION of Constants and Variables - FINISHED //

  //////////////////////////////////////////
  // Examination of different parameters ///
  //////////////////////////////////////////

  g.Stokes_translation = 6.0*PI*g.Stokes_rad*g.Visco;
  g.Stokes_rotation = 8.0*PI*pow(g.Stokes_rad, 3)*g.Visco/1.0; /* the formula may be adjusted for easy rotation */

  /**************** INITIALIZATION *********************/
  // center of the nucleus //
  Nuc[0]=Rad-LL,Nuc[1]=0.0,Nuc[2]=0.0; /**************** start point ********************/
  g.DVecNucCen[0][0]=0.0; g.DVecNucCen[0][1]=0.0; g.DVecNucCen[0][2]=LL;      
  g.DVecNucCen[1][0]=0.0; g.DVecNucCen[1][1]=0.0; g.DVecNucCen[1][2]=(-1.0)*LL;

  // 1st rotation
  starting_degree=0.0;
  rotationAx[0]=0.0;
  rotationAx[1]=1.0;
  rotationAx[2]=0.0;
  MakeRotationMatrix(RotationMatrix, rotationAx, starting_degree);

  ProductJacVec(VECVEC, RotationMatrix, g.DVecNucCen[0]);
  for (j=0; j<3; j++) {g.DVecNucCen[0][j]=VECVEC[j];}
  ProductJacVec(VECVEC, RotationMatrix, g.DVecNucCen[1]);
  for (j=0; j<3; j++) {g.DVecNucCen[1][j]=VECVEC[j];}

  for (j=0; j<3; j++) {
    PVecCen[0][j] = Nuc[j]+g.DVecNucCen[0][j];
    PVecCen[1][j] = Nuc[j]+g.DVecNucCen[1][j];
  }

  // phase and length of MTs
  for (k=0; k<g.N; k++){
    g.L[k] = 0.537e-6*dT;  
    g.phase[k]=1;
    g.u[k][0] = cos(Sdelta[k])*cos(Ssita[k]);
    g.u[k][1] = cos(Sdelta[k])*sin(Ssita[k]);
    g.u[k][2] = sin(Sdelta[k]);
    for (j=0;j<3;j++){
      VECVEC[j]=g.u[k][j];
    }
    ProductJacVec(VECVECVEC, RotationMatrix, VECVEC);
    for (j=0;j<3;j++){
      g.u[k][j]=VECVECVEC[j];
    }
    // calculation of  g.u[k][3], g.u[k][4], g.u[k]5]
    if (k<g.NN) {
      OutProdVector(g.DVecNucCen[0], VECVECVEC, VECVEC);
    } else {
      OutProdVector(g.DVecNucCen[1], VECVECVEC, VECVEC);
    }
    for (j=0;j<3;j++){g.u[k][3+j]=VECVEC[j];}
    g.pushing_phase[k] = 0;
    if (k<g.NN){qq=0;}else{qq = 1;} /* acting on centrosome 1 or 2 */
    for (j=0; j<3; j++) {
      MTC[j] = g.L[k]*g.u[k][j];
      tempMT[k][j] = PVecCen[qq][j] + MTC[j];}
    // calculation of the distance from the centrosome to the cell cortex at the angle //
    AA[0]=pow(MTC[0]*RadS,2)+pow(MTC[1]*Rad,2)+pow(MTC[2]*Rad,2); /* quadratic equation */
    AA[1]=2*(MTC[0]*PVecCen[qq][0]*RadS*RadS+MTC[1]*PVecCen[qq][1]*Rad*Rad+MTC[2]*PVecCen[qq][2]*Rad*Rad);
    AA[2]=pow(RadS*PVecCen[qq][0],2)+pow(Rad*PVecCen[qq][1],2)+pow(Rad*PVecCen[qq][2],2)-pow(Rad*RadS,2);
    QuadEqu2(AA,BB);
    g.L[k] = BB[0]*g.L[k];
  }
  currentN = g.N;
  // currentN = 0;  /*** add when simulation with increasing number of MTs (Sup Fig S3), start = 0 MT ***/

  double phi, xmax, xmin, dx;
  int ix;
  xmax = (Rad-LL)*pow(10,6);
  xmin = (-Rad+LL)*pow(10,6);
  
  ////////////////////////////////////////
  ///// Repetition for each time step ////
  ////////////////////////////////////////

  step_counter = -1;
  for (ix=xmin; ix<=xmax; ix++) {
    for (i=0; i<=180; i++) {
      step_counter++;

      dx = ix * pow(10,-6);
      phi = i*PI/180;
      // center of the nucleus //
      Nuc[0]=dx,Nuc[1]=0.0,Nuc[2]=0.0; /**************** start point ********************/
      //Nuc[0]=0.0,Nuc[1]=0.0,Nuc[2]=0.0; /**************** start point ********************/
      g.DVecNucCen[0][0]=0.0; g.DVecNucCen[0][1]=0.0; g.DVecNucCen[0][2]=LL;      
      g.DVecNucCen[1][0]=0.0; g.DVecNucCen[1][1]=0.0; g.DVecNucCen[1][2]=(-1.0)*LL;

      starting_degree=phi;
      rotationAx[0]=0.0;
      rotationAx[1]=1.0;
      rotationAx[2]=0.0;
      MakeRotationMatrix(RotationMatrix, rotationAx, starting_degree);

      ProductJacVec(VECVEC, RotationMatrix, g.DVecNucCen[0]);
      for (j=0; j<3; j++) {g.DVecNucCen[0][j]=VECVEC[j];}
      ProductJacVec(VECVEC, RotationMatrix, g.DVecNucCen[1]);
      for (j=0; j<3; j++) {g.DVecNucCen[1][j]=VECVEC[j];}

      for (j=0; j<3; j++) {
        PVecCen[0][j] = Nuc[j]+g.DVecNucCen[0][j];
        PVecCen[1][j] = Nuc[j]+g.DVecNucCen[1][j];
      }

      // phase and length of MTs
      for (k=0; k<g.N; k++){
        g.L[k] = 0.537e-6*dT;  
        g.phase[k]=1;
        g.u[k][0] = cos(Sdelta[k])*cos(Ssita[k]);
        g.u[k][1] = cos(Sdelta[k])*sin(Ssita[k]);
        g.u[k][2] = sin(Sdelta[k]);
        for (j=0;j<3;j++){
          VECVEC[j]=g.u[k][j];
        }
        ProductJacVec(VECVECVEC, RotationMatrix, VECVEC);
        for (j=0;j<3;j++){
          g.u[k][j]=VECVECVEC[j];
        }
        // calculation of  g.u[k][3], g.u[k][4], g.u[k]5]
        if (k<g.NN) {
          OutProdVector(g.DVecNucCen[0], VECVECVEC, VECVEC);
        } else {
          OutProdVector(g.DVecNucCen[1], VECVECVEC, VECVEC);
        }
        for (j=0;j<3;j++){g.u[k][3+j]=VECVEC[j];}
        g.pushing_phase[k] = 0;
        if (k<g.NN){qq=0;}else{qq = 1;} /* acting on centrosome 1 or 2 */
        for (j=0; j<3; j++) {
          MTC[j] = g.L[k]*g.u[k][j];
          tempMT[k][j] = PVecCen[qq][j] + MTC[j];}
        // calculation of the distance from the centrosome to the cell cortex at the angle //
        AA[0]=pow(MTC[0]*RadS,2)+pow(MTC[1]*Rad,2)+pow(MTC[2]*Rad,2); /* quadratic equation */
        AA[1]=2*(MTC[0]*PVecCen[qq][0]*RadS*RadS+MTC[1]*PVecCen[qq][1]*Rad*Rad+MTC[2]*PVecCen[qq][2]*Rad*Rad);
        AA[2]=pow(RadS*PVecCen[qq][0],2)+pow(Rad*PVecCen[qq][1],2)+pow(Rad*PVecCen[qq][2],2)-pow(Rad*RadS,2);
        QuadEqu2(AA,BB);
        g.L[k] = BB[0]*g.L[k];
      }
    
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

      /////////////////////  growth/shrinkage for each microtubule ////////////////////
      for (k=0; k<g.N; k++) {
        g.NumberOfMotor[k] = 0.0; /* initialization */
        if (k<g.NN){qq=0;}else{qq = 1;} /* acting on centrosome 1 or 2 */
        g.pushing_phase[k] = 0; /* initial assumption = no pushing */
        switch (g.phase[k]) { /* phase1=growing, phase0=shrinking */
          case 2: /* inactive MTs */
            g.pushing_phase[k] = 0;
            g.pulling_phase[k] = 0;
            break;
          case 1: /* growing phase */
            g.pulling_phase[k] = 1;
            previousL[k] = g.L[k];
            g.L[k] += g.Vg*dT;
            break;
        }
        g.pushing_phase[k]=1;
        
        // LENGTH-DEPENDENT PULLING + CORTEX-PULLING
        if (CortexPullingMode[k]==0) {
          for (j=0; j<3; j++) {
            tempMT[k][j] = PVecCen[qq][j] + g.L[k]*g.u[k][j];
          }
          switch (strain)
          {
            case 0:  // WT
              if (tempMT[k][0]<0) {
                CortexPullingFreq = CortexPullingFreq_PAR3;
              } else {
                CortexPullingFreq = CortexPullingFreq_LET99+(tempMT[k][0]/Rad)*(CortexPullingFreq_PAR2-CortexPullingFreq_LET99);
              }
              break;
            case 1: // par-2 w/o LET-99
              CortexPullingFreq = CortexPullingFreq_PAR3;
              break;
            case 2: // par-2 with LET-99
              /* if ((i<metaST)&&(tempMT[k][0]>0)) { */
              if (tempMT[k][0]>0) {
                CortexPullingFreq = ((Rad-tempMT[k][0])/Rad)*CortexPullingFreq_PAR3;
              } else {
                CortexPullingFreq = CortexPullingFreq_PAR3;
              }
              break;
            case 3: // par-3 w/o LET-99
              CortexPullingFreq = CortexPullingFreq_PAR2;
              break;
            case 4: // par-3 with LET-99
              CortexPullingFreq = ((fabs(tempMT[k][0]))/Rad)*CortexPullingFreq_PAR2;
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
          if (1<=pow((PI*g.L[k])/(4*MTDivision),2)*PI*g.CortexMotorDensity*CortexPullingFreq) {
            g.CortexMotorNumber[k] = 1;}
          else {
            g.CortexMotorNumber[k] = pow((PI*g.L[k])/(4*MTDivision),2)*PI*g.CortexMotorDensity*CortexPullingFreq;}
          
          CortexPullingMode[k]=CortexPullingDuration; /********** no duration ***************/
        }
        if (CortexPullingMode[k]>0){
          g.NumberOfMotor[k] += g.CortexMotorNumber[k];
          CortexPullingMode[k]--;
        }
        
        g.pushing_phase[k]=0; // To ensure that pushing_phase of the MT without contacting the cortex is 0 //
        g.CortexMotorNumber[k] = 0;
        
        ///////////////////////////////////////////
        // LENGTH-DEPENDENT PULLING + CORTEX PULLING
        g.NumberOfMotor[k] += g.L[k]*g.MotorDensity;
        //TRACE(("%4d %4d N=%3.1lf MT=(%3.2lf %3.2lf %3.2lf %3.2lf %3.2lf %3.2lf)\n",i,k,NumberOfMotor[k],u[k][0],u[k][1],u[k][2],u[k][3],u[k][4],u[k][5]));
        for (j=0;j<6;j++){
          ForceC[qq][j] += g.NumberOfMotor[k] * g.u[k][j] * g.MotorStallF;
        } /* an initial guess */

        // to monitor microtubules profile
        sum += g.L[k];
        if (max < g.L[k]){max = g.L[k];}
        if (min > g.L[k]){min = g.L[k];}
      }

      ///////////////// Translational and rotational movement of the pronucleus //////////////////////

      //for calculation of the velocity of pronucleus 
      if (i%UNITTIME==0){
        /* if (i==0){ */
        OldDis = LL;
        OldVel = 0.0;
      }
      // initialization
      for (j=1; j<=6; j++){
        tempNucVel[j]=0.0;
      }

      ///////////////////////////////////////////////////////////////////////
      // THE PULLING MODEL: solve the set of equation using Newton-Raphson method with the initial guess
      for (j=1;j<=6;j++){ /* initial value of tempNucVel[j] */
        TRACE(("%4d Force[%d]: %10lf %10lf\n",i, j, ForceC[0][j-1]*1.0e+12, ForceC[1][j-1]*1.0e+12));
        if (j<=3) {
          tempNucVel[j] = (ForceC[0][j-1]+ForceC[1][j-1])/g.Stokes_translation; 
        } else {
          tempNucVel[j] = (ForceC[0][j-1]+ForceC[1][j-1])/g.Stokes_rotation;
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
      for (k=0; k<g.N; k++) {eachMT_PTC[k]=0;}
      do {
        g.phase_transition_count = 0;
        for (j=0;j<6;j++){
          g.Fbackward[j] = 0.0;
          for (jj=0; jj<6; jj++) {
            g.fjac_pull[j][jj] = 0.0;
          }}
        for (k=0; k<g.N; k++) {
          if (g.pulling_phase[k]!=0){
            if (k<g.NN) {qq = 0;} else {qq = 1;}
            dv = 0.0;
            for (j=0; j<3; j++) {
              dv += (tempNucVel[(j+1)%3+4]*g.DVecNucCen[qq][(j+2)%3]-tempNucVel[(j+2)%3+4]*g.DVecNucCen[qq][(j+1)%3]+tempNucVel[j+1])*g.u[k][j];//6/3
            }
            if (dv > g.MotorMaxVel){ /* the motors on this MT do not exert forces: pulling_phase[k]=3 */
              if (g.pulling_phase[k]!=3) {
                g.phase_transition_count++;
                eachMT_PTC[k]++;
              }
              g.pulling_phase[k] = 3;
            } else {
              if (dv < 0) { /* the motors on this MT exert the maximum (stall) force: pulling_phase[k]=2 */
                if (g.pulling_phase[k]!=2) {
                  g.phase_transition_count++;
                  eachMT_PTC[k]++;
                }
                g.pulling_phase[k] = 2;
                for (j=0;j<6;j++) {
                  g.Fbackward[j] += g.u[k][j] * g.MotorStallF * g.NumberOfMotor[k];
                }
              } else {
                if (g.pulling_phase[k]!=1) { /* the motors on this MT exert force dependent on their velocity: pulling_phase[k]=1 */
                  g.phase_transition_count++;
                  eachMT_PTC[k]++;
                }
                g.pulling_phase[k] = 1;
                for (j=0; j<3; j++) {dvid[j] = g.u[k][j];}
                dvid[3] = g.DVecNucCen[qq][1]*g.u[k][2] - g.DVecNucCen[qq][2]*g.u[k][1];
                dvid[4] = g.DVecNucCen[qq][2]*g.u[k][0] - g.DVecNucCen[qq][0]*g.u[k][2];
                dvid[5] = g.DVecNucCen[qq][0]*g.u[k][1] - g.DVecNucCen[qq][1]*g.u[k][0];
                for (j=0; j<6; j++) {
                  for (jj=0; jj<6; jj++) {
                    g.fjac_pull[j][jj] -= (g.MotorStallF*g.NumberOfMotor[k]/g.MotorMaxVel)*dvid[jj]*g.u[k][j];
                  }
                }
              }
            }
          }
        }
        for (j=0; j<3; j++) {g.fjac_pull[j][j] -= g.Stokes_translation;}
        for (j=3; j<6; j++) {g.fjac_pull[j][j] -= g.Stokes_rotation;}
        if ((g.phase_transition_count!=0)||(cycle_count==0)) {
          did_converge = mnewt(10, tempNucVel, 6, tolx, tolf, step_counter, f_out_3dcheck, usr_func, &g);
        }
        if (i%100==0) fprintf(f_out_3dcheck,"%d %d %d\n", i, cycle_count, g.phase_transition_count);
        cycle_count++;
      } while ((cycle_count<=1)||((g.phase_transition_count!=0)&&(cycle_count<1000))); /* repeat until the solution satisfies all equations and conditions */

      TRACE(("%4d %4d ",i,cycle_count));
      for (j=1; j<=6; j++) {
        if (j<=3) {
          TRACE(("%4.3lf ",tempNucVel[j]*1.0e+6));
        } else {
          TRACE(("%4.3lf ",tempNucVel[j]*100));
        }
      }
      TRACE(("\n"));
      /* if (i%10==0 && ix%5==0) printf("%d %.3lf %.3lf %.1lf\n", ix, phi*180/PI, tempNucVel[1]*1.0e+7, tempNucVel[5]*1000); */
      if (i%10==0 && ix%5==0) fprintf(f_out_vector, "%.3lf %.3lf %.3lf %.1lf\n", ix*100/(Rad*pow(10,6)), phi*180/PI, 2*tempNucVel[1]*1.0e+7, 3*tempNucVel[5]*10);

      if (g.phase_transition_count!=0){ /* when the solution is not obtained within 1000 cycles*/
        printf("exit without convergence at t=%d PTC=%d\n", i, g.phase_transition_count);
        fprintf(f_out_fvcheck,"exit without convergence at t=%d PTC=%d\n", i, g.phase_transition_count);
        printf("Cen1=(%3.1f,%3.1f,%3.1f) Cen2=(%3.1f,%3.1f,%3.1f)\n",PVecCen[0][0]*1.0e+6,PVecCen[0][1]*1.0e+6,PVecCen[0][2]*1.0e+6,PVecCen[1][0]*1.0e+6,PVecCen[1][1]*1.0e+6, PVecCen[1][2]*1.0e+6);
        break;
      }

      // Monitoring local migration of pronucleus
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
          /* XFillRectangle(mtg.d,mtg.w9,mtg.gc9,(int)(WIN_WIDTH/2+500*1000000*DeltaNucAxis2/(dT*UNITTIME))-1,(int)(WIN_HEIGHT/2-500*1000000*DeltaNucAxis1/(dT*UNITTIME))-1,2,2); */
        }
      }
      for (j=0;j<3;j++) {
        OldNuc[j] = Nuc[j];
      }

      // TRANSLATIONAL MOVEMENT OF THE PRONUCLEUS
      CalculateVel = 0.0;
      for (j=1; j<=3; j++) {
        CalculateVel += tempNucVel[j]*tempNucVel[j];
      }

      // the pronucleus does not cross over the cell cortex 
      NuclearDistanceRatio = sqrt(pow((Nuc[0]/(Rad-LL)),2)+pow((Nuc[1]/(RadS-LL)),2)+pow((Nuc[2]/(RadS-LL)),2));
      // calculation of velocity of the pronucleus
      if(i%UNITTIME==0){
        NewDis = sqrt((Rad-Nuc[0])*(Rad-Nuc[0])+Nuc[1]*Nuc[1]+Nuc[2]*Nuc[2]);
        NewVel = (NewDis-OldDis)/(dT*UNITTIME);
      }
      DistanceFromPP = sqrt((Rad-Nuc[0])*(Rad-Nuc[0])+Nuc[1]*Nuc[1]+Nuc[2]*Nuc[2]);
    }
  }
  //////////////// examination with single parameter set FINISHED ///////////////////////////
  XStoreName(mtg.d,mtg.w1,"fin");
  /* sprintf (filefin, "xwd -name \'fin\' -out outFIN%d",p); */
  system (filefin);
  /* sprintf (filefin2, "convert outFIN%d outFIN%d.jpg",p,p); */
  system (filefin2);
  ///////////// examination with all parameter sets FINISHED ///////////////////
  fclose(f_out1);
  fclose(f_out2);
  fclose(f_out3);
  fclose(f_out4);
  fclose(f_out5);
  /* fclose(f_out_param); */
  /* fclose(f_out_mtnum); */
  /* fclose(f_out_for3d); */
  fclose(f_out_fvcheck);
  fclose(f_out_3dcheck);
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

  free_dmatrix(g.u,0,g.N-1,0,2);
  free_cvector(g.pushing_phase,0,g.N-1);
  free_cvector(g.pulling_phase,0,g.N-1);
  free_cvector(g.phase,0,g.N-1);
  free_dvector(g.L,0,g.N-1);

  printf("to end, hit [return] key:");
  getchar();
  return 0;
}
