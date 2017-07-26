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
  unsigned int model = 0; /* 0:MT_angle_fixed, 1:MT_angle_moveable */

  /* Parse options */
  myname = argv[0];
  while ((ch = getopt(argc, argv, "m:cvh")) != -1){
    switch (ch) {
      case 'c':
        is_check = true;
        break;
      case 'v':
        is_verbose = true;
        break;
      case 'm':
        model = atoi(optarg);
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

  if(model > 2){
    usage(myname);
    exit(1);
  }

  //file handling
  FILE *f_out_fit = fopen("out_fit.csv","w");
  FILE *f_out_param = fopen("out_parameter.dat","w");
  FILE *f_out_mtnum = fopen("out_MTnumber.dat","w");
  FILE *f_out_fvcheck = fopen("out_FVcheck.dat","w");
  FILE *f_out_3dcheck = fopen("out_3Dcheck.dat","w");
  FILE *f_out_for3d = fopen("out_for_3D.dat","w");

  ////////////////////////////////////////////
  // DECLEARATION of Constants and Variables//
  ////////////////////////////////////////////

  /************CONSTANT NUMBER*****************/
  //PARAMETERS
  // drag force of pronucleus
  g.Visco = 1.0; /* viscosity of cytosol [kg/m sec], standard 0.1->1.0 */
  /* g.Stokes_rad = 10.0 * pow(10,-6); */
  g.Stokes_rad = 1.5 * pow(10,-6);
  // pulling force
  double ForceCoef3 = 0.0;
  double ForceCoef2 = 1.5 * pow(10,-2); /* the elongation continues for 4min */
  if (model==0) {
    ForceCoef2 = 0.25 * pow(10,-2);
  }
  else if (model==1) {
    ForceCoef2 = 0.09 * pow(10,-2);
  }
  double ForceCoef1 = 0.0;

  // arrangement of microtubules: calculation of the (maximum) number of MTs
  /* double starting_degree = 0.0; /\* starting degree of aster: for the simulation with single aster (Sup Fig S9) this value should be changed*\/ */
  double starting_degree = PI/2.0; /* starting degree of aster: for the simulation with single aster (Sup Fig S9) this value should be changed*/
  int MTInitAngle_degree = 110; /* the maximum angle of MTs [degree] */
  /* int MTDiv90 = 20; /\* distribution of MT in 90 degree [/degree] *\/ */
  int MTDiv90 = 6; /* distribution of MT in 90 degree [/degree] */
  double MTInitAngle = MTInitAngle_degree * PI / 180;
  double MTAngleDens = MTDiv90/90.0; /* the number of MTs in unit degree [/degree] */
  int MTDivision = (int)(MTInitAngle_degree * MTAngleDens); /* this number defines the (maximum) number of MTs*/
  
  double MetaSpindle_L = 14*pow(10,-6);
  double Rad = 23 * pow(10,-6);  /* long axis of egg [m] */
  double RadS = 14 * pow(10,-6); /* short axes of egg [m] */
  double aspect_ratio = Rad/RadS;
  double Exp_Elongated_L = 0.39 * (2*Rad) + 5.2 * pow(10,-6);

  /************VARIABLES**********************/
  int p=0; /* different parameter sets */
  int j; /* x,y,z axes */
  int qq; /* centrosome 0 or 1 */
  double VECVEC[3], VECVECVEC[3];
  // microtubules //
  /* double MT[g.N][3]; /\* rectangular coordinates of n-th MT *\/ */
  double MTC[3]; /* for calculation of length to cortex used around L448 */
  // first and second centrosome//
  g.DVecNucCen = dmatrix(0,1,0,2);
  double PVecCen[2][3];
  double Nuc[3];  /* position vector of center of the nucleus */
  // forces //
  double ForceC[2][3]; /* force vector on 1st and 2nd centrosome */
  // translational and rotational movement of the pronucleus //
  double DirectionDetermination[3];
  int step_counter;
  ////////////////////////ROTATION
  double rotationAx[3]; /* rotational axis */
  ////////////////////// Nucleus does not cross over the cortex /////////////////
  double AA[3]; /* three coefficients of the equation */
  double BB[2]; /* two solutions of the equation */
  // variables for local movement
  double CenVel[2][4];
  double RotationMatrix[3][3];

  char filefin[40];
  char filefin2[40];

  //// DECLARATION of Constants and Variables - FINISHED //

  //////////////////////////////////////////
  // Examination of different parameters ///
  //////////////////////////////////////////

  /* printf ("#MTDivision,NN,MTInitAngle_degree,spindle_length,fabs(exp-sim),time[sec]\n"); */
  printf ("#MTInitAngle_degree,MTDiv90,MTDivision,NN,spindle_length,exp-sim,time[sec]\n");
  fprintf (f_out_fit,"#MTInitAngle_degree,MTDiv90,MTDivision,NN,spindle_length,exp-sim,time[sec]\n");
  /* for (MTInitAngle_degree=110; MTInitAngle_degree<150; MTInitAngle_degree++) { */
  /* for (MTInitAngle_degree=120; MTInitAngle_degree<135; MTInitAngle_degree++) { */
  for (MTInitAngle_degree=91; MTInitAngle_degree<170; MTInitAngle_degree++) {
    /* for (MTDiv90=2; MTDiv90<30; MTDiv90++) { */
    for (MTDiv90=6; MTDiv90<7; MTDiv90++) {
      MTAngleDens = MTDiv90/90.0;
      MTDivision = (int)(MTInitAngle_degree * MTAngleDens);
      MTInitAngle = MTInitAngle_degree * PI / 180;
      int MTPerPlane[MTDivision+1];
      int MTTotal[MTDivision+1];
      int kk; /* plane number */
      int k; /* MT number */
      MTPerPlane[0]=1;
      MTTotal[0]=1;
      for (kk=1; kk<MTDivision+1; kk++) {
        MTPerPlane[kk]=(int)(2*PI*MTDivision*sin(MTInitAngle*kk/MTDivision)/MTInitAngle);
        MTTotal[kk]=MTTotal[kk-1]+MTPerPlane[kk];
        /* printf("%d %d %d\n",MTDivision,MTPerPlane[kk],MTTotal[kk]); */
      }
      g.NN = MTTotal[MTDivision];
      g.N = 2*g.NN;
      /* printf("NN=%d N=%d\n", g.NN, g.N); */
      // allocate memory
      g.u = dmatrix(0,g.N-1,0,2); ///////////////////////////////
      g.L = dvector(0,g.N-1);
      // arrangement of microtubules: direction of each MT
      double Sdelta[g.N], Ssita[g.N];
      Sdelta[0]=0.5*PI; Ssita[0]=0.0;
      kk=0;
      for (k=1;k<g.NN; k++) {
        if (k > MTTotal[kk]-1){kk++;}
        Sdelta[k]=PI/2.0-MTInitAngle*kk/MTDivision;
        Ssita[k]=2.0*PI*k/(MTPerPlane[kk]);
      }
      Sdelta[g.NN] = -0.5*PI; Ssita[0]=0.0;
      kk=0;
      for (k=g.NN+1;k<g.N; k++) {
        if ((k-g.NN) > MTTotal[kk]-1){kk++;}
        Sdelta[k]=-PI/2.0+MTInitAngle*kk/MTDivision;
        Ssita[k]=2.0*PI*(k-g.NN)/(MTPerPlane[kk]);
      }
      double MT[g.N][3]; /* rectangular coordinates of n-th MT */
      
      g.Stokes_translation = 6.0*PI*g.Stokes_rad*g.Visco;
      /* printf("asp=%lf, Rad=%e, RadS=%e\n",aspect_ratio,Rad,RadS); */
      
      // OUTPUT parameter LOGs
      fprintf(f_out_param,"p=%d\nmodel=%d\nMTDivision=%d MTInitAngle=%lf MT=%d\naspect_ratio=%lf Rad=%lf RadS=%lf\nH=%5.3lf Stokes_rad=%5.3lf\nForceCoef1=%lf Coef2=%lf Coef3=%lf\n\n", p, model, MTDivision, MTInitAngle, g.N, aspect_ratio, Rad*pow(10,6), RadS*pow(10,6), g.Visco, g.Stokes_rad*1.0e+6, ForceCoef1, ForceCoef2, ForceCoef3);
      
      // color settings
      /* #include "color_setting.c" */
      color_setting(&mtg, p);
    
      /**************** INITIALIZATION *********************/
      // center of the nucleus //
      /* Nuc[0]=Rad-LL,Nuc[1]=0.0,Nuc[2]=0.0; /\**************** start point ********************\/ */
      Nuc[0]=0.0,Nuc[1]=0.0,Nuc[2]=0.0; /**************** start point ********************/
      g.DVecNucCen[0][0]=0.0; g.DVecNucCen[0][1]=0.0; g.DVecNucCen[0][2]=MetaSpindle_L/2;      
      g.DVecNucCen[1][0]=0.0; g.DVecNucCen[1][1]=0.0; g.DVecNucCen[1][2]=(-1.0)*MetaSpindle_L/2;

      // 1st rotation
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
    
      // force
      for (j=0;j<3;j++) {
        for (qq=0; qq<2; qq++){
          ForceC[qq][j]=0.0;
        }
      }
      // phase and length of MTs
      for (k=0; k<g.N; k++){
        if (k<g.NN){qq=0;}else{qq=1;}
        g.L[k] = 1.0*pow(10,-9);  
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
        for (j=0; j<3; j++) {
          MTC[j] = g.L[k]*g.u[k][j];}
        // calculation of the distance from the centrosome to the cell cortex at the angle //
        AA[0]=pow(MTC[0]*RadS,2)+pow(MTC[1]*Rad,2)+pow(MTC[2]*Rad,2); /* quadratic equation */
        AA[1]=2*(MTC[0]*PVecCen[qq][0]*RadS*RadS+MTC[1]*PVecCen[qq][1]*Rad*Rad+MTC[2]*PVecCen[qq][2]*Rad*Rad);
        AA[2]=pow(RadS*PVecCen[qq][0],2)+pow(Rad*PVecCen[qq][1],2)+pow(Rad*PVecCen[qq][2],2)-pow(Rad*RadS,2);
        QuadEqu2(AA,BB);
        g.L[k] = BB[0]*g.L[k];
        for (j=0;j<3;j++){
          MT[k][j] = PVecCen[qq][j] + g.L[k]*g.u[k][j];
        }
      }

      ////////////////////////////////////////
      ///// Repetition for each time step ////
      ////////////////////////////////////////
    
      step_counter = -1;
      do {
        step_counter++;
        ////////////////////  INITIALIZATION ////////////////////////////
        for (j=0;j<3;j++) {
          for (qq=0; qq<2; qq++){
            ForceC[qq][j]=0.0;
          }
        }
        for (j=0; j<3;j++) {	  
          DirectionDetermination[j]=0.0;
        }

        if (model==0) /* MT_angle_fixed */
        {
          for (k=0; k<g.N; k++) {
            if (k<g.NN){qq=0;}else{qq=1;}
            for (j=0; j<3; j++) {
              MTC[j] = g.L[k]*g.u[k][j];
            }
            // calculation of the distance from the centrosome to the cell cortex at the angle //
            AA[0]=pow(MTC[0]*RadS,2)+pow(MTC[1]*Rad,2)+pow(MTC[2]*Rad,2); /* quadratic equation */
            AA[1]=2*(MTC[0]*PVecCen[qq][0]*RadS*RadS+MTC[1]*PVecCen[qq][1]*Rad*Rad+MTC[2]*PVecCen[qq][2]*Rad*Rad);
            AA[2]=pow(RadS*PVecCen[qq][0],2)+pow(Rad*PVecCen[qq][1],2)+pow(Rad*PVecCen[qq][2],2)-pow(Rad*RadS,2);
            QuadEqu2(AA,BB);
            g.L[k] = BB[0]*g.L[k];
            for (j=0;j<3;j++){
              MT[k][j] = PVecCen[qq][j] + g.L[k]*g.u[k][j];
              ForceC[qq][j] += (ForceCoef3*pow(g.L[k],3)+ForceCoef2*pow(g.L[k],2)+ForceCoef1*g.L[k]) * g.u[k][j];
            }
          }
        }

        else if (model==1) /* MT_angle_variable */
        {
          for (k=0; k<g.N; k++){
            if (k<g.NN){qq=0;}else{qq=1;}
            g.L[k] = sqrt(pow(MT[k][0]-PVecCen[qq][0],2)+pow(MT[k][1]-PVecCen[qq][1],2)+pow(MT[k][2]-PVecCen[qq][2],2));
            for (j=0; j<3; j++) {
              g.u[k][j] = (MT[k][j] - PVecCen[qq][j]) / g.L[k];
              ForceC[qq][j] += (ForceCoef3*pow(g.L[k],3)+ForceCoef2*pow(g.L[k],2)+ForceCoef1*g.L[k]) * g.u[k][j];
            }
          }
        }

        for (j=0; j<3; j++) {
          CenVel[0][j+1] = ForceC[0][j] / g.Stokes_translation;
          CenVel[1][j+1] = ForceC[1][j] / g.Stokes_translation;
          PVecCen[0][j] += CenVel[0][j+1] * dT;
          PVecCen[1][j] += CenVel[1][j+1] * dT;
        }
      
        ///////////////////////// 1 STEP FINISHED ////////////////////////////////////////////
        //OUTPUTS 
        /* draw graphs in X-window */
        draw_graphs(step_counter, &mtg, &g, PVecCen, MT, Nuc, Rad, RadS, MetaSpindle_L);
        /* save logs in texts */
        /* save_logs(step_counter, p, g.N, f_out1, f_out2, f_out3, f_out4, f_out5, f_out6, f_out_for3d, PVecCen, MT); */
        /* } while ((step_counter<100000000) && (0.1*pow(10,-6)<fabs(CenVel[0][0+1]-CenVel[1][0+1])*10)); */
      } while ((step_counter<100000000) && (0.1*pow(10,-10)<fabs(CenVel[0][0+1]-CenVel[1][0+1])*10));

      if (100000000 <= step_counter) {
        printf("%d,%d,%d,%d,%lf,%lf,not_convergence\n", MTInitAngle_degree, MTDiv90, MTDivision, g.NN, fabs(PVecCen[0][0]-PVecCen[1][0])*pow(10,6), (Exp_Elongated_L-fabs(PVecCen[0][0]-PVecCen[1][0]))*pow(10,6));
        fprintf(f_out_fit,"%d,%d,%d,%d,%lf,%lf,not_convergence\n", MTInitAngle_degree, MTDiv90, MTDivision, g.NN, fabs(PVecCen[0][0]-PVecCen[1][0])*pow(10,6), (Exp_Elongated_L-fabs(PVecCen[0][0]-PVecCen[1][0]))*pow(10,6));
      }
      else {
        printf("%d,%d,%d,%d,%lf,%lf,%.3lf\n", MTInitAngle_degree, MTDiv90, MTDivision, g.NN, fabs(PVecCen[0][0]-PVecCen[1][0])*pow(10,6), (Exp_Elongated_L-fabs(PVecCen[0][0]-PVecCen[1][0]))*pow(10,6), step_counter*dT);
        fprintf(f_out_fit,"%d,%d,%d,%d,%lf,%lf,%.3lf\n", MTInitAngle_degree, MTDiv90, MTDivision, g.NN, fabs(PVecCen[0][0]-PVecCen[1][0])*pow(10,6), (Exp_Elongated_L-fabs(PVecCen[0][0]-PVecCen[1][0]))*pow(10,6), step_counter*dT);
      }
      //////////////// examination with single parameter set FINISHED ///////////////////////////
      /* XStoreName(mtg.d,mtg.w1,"fin"); */
      /* sprintf (filefin, "xwd -name \'fin\' -out outFIN%d",p); */
      /* system (filefin); */
      /* sprintf (filefin2, "convert outFIN%d outFIN%d.jpg",p,p); */
      /* system (filefin2); */

    }
    fprintf(f_out_fit,"\n");
  }
  ///////////// examination with all parameter sets FINISHED ///////////////////
  fclose(f_out_fit);
  fclose(f_out_param);
  fclose(f_out_for3d);
  fclose(f_out_mtnum);
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
  free_dvector(g.L,0,g.N-1);

  printf("to end, hit [return] key:");
  getchar();
  return 0;
}
