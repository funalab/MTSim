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
  FILE *f_out1 = fopen("out01.csv","w");
  FILE *f_out2 = fopen("out02.csv","w");
  FILE *f_out3 = fopen("out03.csv","w");
  FILE *f_out4 = fopen("out04.csv","w");
  FILE *f_out5 = fopen("out05.csv","w");
  FILE *f_out_somesize = fopen("out_somesize.csv","w");
  FILE *f_out_param = fopen("out_parameter.dat","w");
  FILE *f_out_mtnum = fopen("out_MTnumber.dat","w");
  FILE *f_out_fvcheck = fopen("out_FVcheck.dat","w");
  FILE *f_out_3dcheck = fopen("out_3Dcheck.dat","w");
  FILE *f_out_for3d = fopen("out_for_3D.csv","w");

  ////////////////////////////////////////////
  // DECLEARATION of Constants and Variables//
  ////////////////////////////////////////////

  /************CONSTANT NUMBER*****************/
  //PARAMETERS
  // drag force of pronucleus
  g.Visco = 1.0; /* viscosity of cytosol [kg/m sec], standard 0.1->1.0 */
  g.Stokes_rad = 1.5 * pow(10,-6);
  // pulling force
  double ForceCoef3 = 0.0;
  double ForceCoef2 = 1.5 * pow(10,-2); /* attraction coefficient [N/m^2] the elongation continues for 4min */
  double ForceCoef1 = 0.0;

  // arrangement of microtubules: calculation of the (maximum) number of MTs
  /* double starting_degree = 0.0; /\* starting degree of aster: for the simulation with single aster (Sup Fig S9) this value should be changed*\/ */
  double starting_degree = PI/2.0; /* starting degree of aster: for the simulation with single aster (Sup Fig S9) this value should be changed*/
  int MTInitAngle_degree = 115; /* the maximum angle of MTs [degree] */
  int MTDiv90 = 6; /* distribution of MT in 90 degree [num] */
  if (model==0) {
    /* MTInitAngle_degree = 132; */
    /* MTDiv90 = 28; */
    ForceCoef2 = 0.25 * pow(10,-2);
    MTInitAngle_degree = 126;
  }
  else if (model==1) {
    /* MTInitAngle_degree = 103; */
    /* MTDiv90 = 28; */
    ForceCoef2 = 0.09 * pow(10,-2);
    MTInitAngle_degree = 97;
  }
  double MTInitAngle = MTInitAngle_degree * PI / 180; /* the maximum angle of MTs [radian] */
  double MTAngleDens = MTDiv90/90.0; /* the number of MTs in unit degree [num/degree] */
  int MTDivision = (int)(MTInitAngle_degree * MTAngleDens); /* this number defines the (maximum) number of MTs [num]*/
  int MTPerPlane[MTDivision+1];
  int MTTotal[MTDivision+1];
  int p_num; /* plane number */
  int MT_num; /* MT number */
  MTPerPlane[0]=1;
  MTTotal[0]=1;
  for (p_num=1; p_num<MTDivision+1; p_num++) {
    MTPerPlane[p_num]=(int)(2*PI*MTDivision*sin(MTInitAngle*p_num/MTDivision)/MTInitAngle);
    MTTotal[p_num]=MTTotal[p_num-1]+MTPerPlane[p_num];
    /* printf("%d %d %d\n",MTDivision,MTPerPlane[p_num],MTTotal[p_num]); */
  }
  g.NN = MTTotal[MTDivision];
  g.N = 2*g.NN;
  printf("NN=%d N=%d\n", g.NN, g.N);
  // allocate memory
  g.u = dmatrix(0,g.N-1,0,2); ///////////////////////////////
  g.L = dvector(0,g.N-1);
  // arrangement of microtubules: direction of each MT
  double Sdelta[g.N], Stheta[g.N];
  Sdelta[0]=0.5*PI; Stheta[0]=0.0;
  p_num=0;
  for (MT_num=1;MT_num<g.NN; MT_num++) {
    if (MT_num > MTTotal[p_num]-1){p_num++;}
    Sdelta[MT_num]=PI/2.0-MTInitAngle*p_num/MTDivision;
    Stheta[MT_num]=2.0*PI*MT_num/(MTPerPlane[p_num]);
  }
  Sdelta[g.NN] = -0.5*PI; Stheta[0]=0.0;
  p_num=0;
  for (MT_num=g.NN+1;MT_num<g.N; MT_num++) {
    if ((MT_num-g.NN) > MTTotal[p_num]-1){p_num++;}
    Sdelta[MT_num]=-PI/2.0+MTInitAngle*p_num/MTDivision;
    Stheta[MT_num]=2.0*PI*(MT_num-g.NN)/(MTPerPlane[p_num]);
  }
  // int MTcoefficient = 100; /* a constant used in simulations with increasing number of MTs (Sup Fig S5) */
  
  double MetaSpindle_L = 14 * pow(10, -6); /* sea-urchin's MetaSpindle_L : 25 * 10^{-6}*/
  double Rad;  /* long axis of egg [m] */
  double RadS; /* short axes of egg [m] */
  double aspect_ratio = 1.0;
  Rad = Cir_Rad * cbrt(aspect_ratio*aspect_ratio);
  RadS = Cir_Rad / cbrt(aspect_ratio);

  /************VARIABLES**********************/
  int p; /* different parameter sets */
  int axis; /* x,y,z axes */
  int qq; /* centrosome 0 or 1 */
  double VECVEC[3], VECVECVEC[3];
  // microtubules //
  double MT[g.N][3]; /* rectangular coordinates of n-th MT */
  double MTC[3]; /* for calculation of length to cortex used around L448 */
  // first and second centrosome//
  g.DVecNucCen = dmatrix(0,1,0,2);
  double PVecCen[2][3];
  double Nuc[3];  /* position vector of center of the nucleus */
  // forces //
  double ForceC[2][3]; /* force vector on 1st and 2nd centrosome */
  // translational and rotational movement of the pronucleus //
  // double DirectionDetermination[3];
  int step_counter;
  ////////////////////////ROTATION
  double rotationAx[3]; /* rotational axis */
  ////////////////////// Nucleus does not cross over the cortex /////////////////
  double QuaFuncCoe[3]; /* three coefficients of the equation */
  double QuaFuncSol[2]; /* two solutions of the equation */
  // variables for local movement
  double CenVel[2][4]; /* velocity of centrosome */
  double RotationMatrix[3][3];
  // variables of aspect ratio
  double ar_init = 1.0;
  double ar_step = 0.5;

  //// DECLARATION of Constants and Variables - FINISHED //

  //////////////////////////////////////////
  // Examination of different parameters ///
  //////////////////////////////////////////

  printf ("#aspect_ratio,Rad,RadS,spindle_length,time[sec],spindle_length/(Rad*2)\n");
  fprintf (f_out_somesize,"#aspect_ratio,Rad,RadS,spindle_length,time[sec],spindle_length/(Rad*2)\n");
	// Mitotic stages are dispersed 
  for (p=0; p<5; p++) {
    aspect_ratio = ar_init + ar_step * p;

    g.Stokes_translation = 6.0*PI*g.Stokes_rad*g.Visco;
    Rad = Cir_Rad * cbrt(aspect_ratio*aspect_ratio);
    RadS = Cir_Rad / cbrt(aspect_ratio);
    MetaSpindle_L = Rad * 2 * (-0.113 * aspect_ratio + 0.4569);

    // OUTPUT parameter LOGs
    fprintf(f_out_param,"p=%d\nmodel=%d\nMTDiv90=%d MTDivision=%d MTInitAngle_degree=%d N=%d\naspect_ratio=%lf Rad=%lf RadS=%lf\nH=%5.3lf Stokes_rad=%5.3lf\nForceCoef1=%lf Coef2=%lf Coef3=%lf\n\n", p, model, MTDiv90, MTDivision, MTInitAngle_degree, g.N, aspect_ratio, Rad*pow(10,6), RadS*pow(10,6), g.Visco, g.Stokes_rad*1.0e+6, ForceCoef1, ForceCoef2, ForceCoef3);

    // color settings
    /* #include "color_setting.c" */
    color_setting(&mtg, p);

    /**************** INITIALIZATION *********************/
    // center of the nucleus //
    /* Nuc[0]=Rad-LL,Nuc[1]=0.0,Nuc[2]=0.0; /\**************** start point ********************\/ */
    Nuc[0]=0.0,Nuc[1]=0.0,Nuc[2]=0.0; /**************** start point ********************/
    g.DVecNucCen[0][0]=0.0; g.DVecNucCen[0][1]=0.0; g.DVecNucCen[0][2]=MetaSpindle_L/2;      
    g.DVecNucCen[1][0]=0.0; g.DVecNucCen[1][1]=0.0; g.DVecNucCen[1][2]=(-1.0)*MetaSpindle_L/2;

    // 1st rotation (rotation axis is y-axis)
    rotationAx[0]=0.0;
    rotationAx[1]=1.0;
    rotationAx[2]=0.0;
    // Rotation from z-axis to x-axis at -PI/2
    MakeRotationMatrix(RotationMatrix, rotationAx, starting_degree);

    ProductJacVec(VECVEC, RotationMatrix, g.DVecNucCen[0]);
    for (axis=0; axis<3; axis++) {g.DVecNucCen[0][axis]=VECVEC[axis];}
    ProductJacVec(VECVEC, RotationMatrix, g.DVecNucCen[1]);
    for (axis=0; axis<3; axis++) {g.DVecNucCen[1][axis]=VECVEC[axis];}

    for (axis=0; axis<3; axis++) {
      PVecCen[0][axis] = Nuc[axis]+g.DVecNucCen[0][axis];
      PVecCen[1][axis] = Nuc[axis]+g.DVecNucCen[1][axis];
    }
    
    // force
    for (axis=0;axis<3;axis++) {
      for (qq=0; qq<2; qq++){
        ForceC[qq][axis]=0.0;
      }
    }
    // phase and length of MTs
    for (MT_num=0; MT_num<g.N; MT_num++){
      if (MT_num<g.NN){qq=0;}else{qq=1;}
      g.L[MT_num] = 1.0*pow(10,-9);  
      g.u[MT_num][0] = cos(Sdelta[MT_num])*cos(Stheta[MT_num]);
      g.u[MT_num][1] = cos(Sdelta[MT_num])*sin(Stheta[MT_num]);
      g.u[MT_num][2] = sin(Sdelta[MT_num]);
      for (axis=0;axis<3;axis++){
        VECVEC[axis]=g.u[MT_num][axis];
      }
      ProductJacVec(VECVECVEC, RotationMatrix, VECVEC);
      for (axis=0;axis<3;axis++){
        g.u[MT_num][axis]=VECVECVEC[axis];
      }
      for (axis=0; axis<3; axis++) {
        MTC[axis] = g.L[MT_num]*g.u[MT_num][axis];}
      // calculation of the distance from the centrosome to the cell cortex at the angle //
      QuaFuncCoe[0]=pow(MTC[0]*RadS,2)+pow(MTC[1]*Rad,2)+pow(MTC[2]*Rad,2); /* quadratic equation */
      QuaFuncCoe[1]=2*(MTC[0]*PVecCen[qq][0]*RadS*RadS+MTC[1]*PVecCen[qq][1]*Rad*Rad+MTC[2]*PVecCen[qq][2]*Rad*Rad);
      QuaFuncCoe[2]=pow(RadS*PVecCen[qq][0],2)+pow(Rad*PVecCen[qq][1],2)+pow(Rad*PVecCen[qq][2],2)-pow(Rad*RadS,2);
      QuadEqu2(QuaFuncCoe,QuaFuncSol);
      g.L[MT_num] = QuaFuncSol[0]*g.L[MT_num];
      for (axis=0;axis<3;axis++){
        MT[MT_num][axis] = PVecCen[qq][axis] + g.L[MT_num]*g.u[MT_num][axis];
      }
    }

    ////////////////////////////////////////
    ///// Repetition for each time step ////
    ////////////////////////////////////////
    
    step_counter = -1;
    do {
      step_counter++;
      ////////////////////  INITIALIZATION ////////////////////////////
      for (axis=0;axis<3;axis++) {
        for (qq=0; qq<2; qq++){
          ForceC[qq][axis]=0.0;
        }
      }
      // for (axis=0; axis<3;axis++) {	  
      //   DirectionDetermination[axis]=0.0;
      // }

      if (model==0) /* MT_angle_fixed */
      {
        for (MT_num=0; MT_num<g.N; MT_num++) {
          if (MT_num<g.NN){qq=0;}else{qq=1;}
          for (axis=0; axis<3; axis++) {
            MTC[axis] = g.L[MT_num]*g.u[MT_num][axis];
          }
          // calculation of the distance from the centrosome to the cell cortex at the angle //
          QuaFuncCoe[0]=pow(MTC[0]*RadS,2)+pow(MTC[1]*Rad,2)+pow(MTC[2]*Rad,2); /* quadratic equation */
          QuaFuncCoe[1]=2*(MTC[0]*PVecCen[qq][0]*RadS*RadS+MTC[1]*PVecCen[qq][1]*Rad*Rad+MTC[2]*PVecCen[qq][2]*Rad*Rad);
          QuaFuncCoe[2]=pow(RadS*PVecCen[qq][0],2)+pow(Rad*PVecCen[qq][1],2)+pow(Rad*PVecCen[qq][2],2)-pow(Rad*RadS,2);
          QuadEqu2(QuaFuncCoe,QuaFuncSol);
          g.L[MT_num] = QuaFuncSol[0]*g.L[MT_num];
          for (axis=0;axis<3;axis++){
            MT[MT_num][axis] = PVecCen[qq][axis] + g.L[MT_num]*g.u[MT_num][axis];
            ForceC[qq][axis] += (ForceCoef3*pow(g.L[MT_num],3)+ForceCoef2*pow(g.L[MT_num],2)+ForceCoef1*g.L[MT_num]) * g.u[MT_num][axis];
          }
        }
      }

      else if (model==1) /* MT_angle_variable */
      {
        for (MT_num=0; MT_num<g.N; MT_num++){
          if (MT_num<g.NN){qq=0;}else{qq=1;}
          g.L[MT_num] = sqrt(pow(MT[MT_num][0]-PVecCen[qq][0],2)+pow(MT[MT_num][1]-PVecCen[qq][1],2)+pow(MT[MT_num][2]-PVecCen[qq][2],2));
          for (axis=0; axis<3; axis++) {
            g.u[MT_num][axis] = (MT[MT_num][axis] - PVecCen[qq][axis]) / g.L[MT_num];
            ForceC[qq][axis] += (ForceCoef3*pow(g.L[MT_num],3)+ForceCoef2*pow(g.L[MT_num],2)+ForceCoef1*g.L[MT_num]) * g.u[MT_num][axis];
          }
        }
      }

      for (axis=0; axis<3; axis++) {
        CenVel[0][axis+1] = ForceC[0][axis] / g.Stokes_translation;
        CenVel[1][axis+1] = ForceC[1][axis] / g.Stokes_translation;
        PVecCen[0][axis] += CenVel[0][axis+1] * dT;
        PVecCen[1][axis] += CenVel[1][axis+1] * dT;
      }
      
      ///////////////////////// 1 STEP FINISHED ////////////////////////////////////////////
      //OUTPUTS 
      /* draw graphs in X-window */
      draw_graphs(step_counter, &mtg, &g, PVecCen, MT, Nuc, Rad, RadS, MetaSpindle_L);
      /* save logs in texts */
      save_logs(step_counter, p, g.N, f_out1, f_out2, f_out3, f_out4, f_out5, f_out_for3d, PVecCen, MT);
    /* } while ((step_counter<100000000) && (0.1*pow(10,-6)<fabs(CenVel[0][0+1]-CenVel[1][0+1])*10)); */
    } while ((step_counter<100000000) && (0.1*pow(10,-10)<fabs(CenVel[0][0+1]-CenVel[1][0+1])*10));

    if (100000000 <= step_counter) {
      printf("%lf,%.3lf,not_convergence\n", aspect_ratio, fabs(PVecCen[0][0]-PVecCen[1][0])*pow(10,6));
    }
    else {
      printf("%.1lf,%.3lf,%.3lf,%.3lf,%lf,%.3lf\n", aspect_ratio, Rad*pow(10,6), RadS*pow(10,6), fabs(PVecCen[0][0]-PVecCen[1][0])*pow(10,6), step_counter*dT, fabs(PVecCen[0][0]-PVecCen[1][0])/(Rad*2));
      fprintf(f_out_somesize,"%.1lf,%.3lf,%.3lf,%.3lf,%lf,%.3lf\n", aspect_ratio, Rad*pow(10,6), RadS*pow(10,6), fabs(PVecCen[0][0]-PVecCen[1][0])*pow(10,6), step_counter*dT, fabs(PVecCen[0][0]-PVecCen[1][0])/(Rad*2));
    }
    //////////////// examination with single parameter set FINISHED ///////////////////////////
    /* XStoreName(mtg.d,mtg.w1,"fin"); */
    /* sprintf (filefin, "xwd -name \'fin\' -out outFIN%d",p); */
    /* system (filefin); */
    /* sprintf (filefin2, "convert outFIN%d outFIN%d.jpg",p,p); */
    /* system (filefin2); */
  }
  ///////////// examination with all parameter sets FINISHED ///////////////////
  fclose(f_out1);
  fclose(f_out2);
  fclose(f_out3);
  fclose(f_out4);
  fclose(f_out5);
  fclose(f_out_somesize);
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
  // store_graphs();

  free_dmatrix(g.u,0,g.N-1,0,2);
  free_dvector(g.L,0,g.N-1);

  printf("to end, hit [return] key:");
  getchar();
  return 0;
}
