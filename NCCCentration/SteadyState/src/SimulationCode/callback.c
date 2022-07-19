/*
 * filename: callback.c
 * Author: Akatsuki Kimura <akkimura@nig.ac.jp>
 *         Akira Funahashi <funa@bio.keio.ac.jp>
 * Last modified: Thu, 25 Jul 2013 02:57:16 +0900
 */
#include "mtsim.h"

void function_FV3D(double *x, int n, double *fvec, double **fjac, mtGlobal* g) {
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
  for (mt=0; mt<g->N; mt++) {
    if (g->pushing_phase[mt]==2) {
      C = g->u[mt][0]*x[1]+g->u[mt][1]*x[2]+g->u[mt][2]*x[3];
      D = C + g->Vg;
      E = D/g->k_on;
      A = 1-E;
      if (A<=0.0){
        g->pushing_phase[mt] = 1;
        g->phase_transition_count++;
        TRACE(("A<0 at function_FV3D\n"));
        for (i=1; i<=n; i++) {
          fvec[i] -= g->u[mt][i-1]*g->F_dependency*g->BucklingConst/(g->L[mt]*g->L[mt]);
        }
      } else {
        B=log(A);
        for (i=1; i<=n; i++) {
          fvec[i] += g->u[mt][i-1]*B;
          for (j=1; j<=n; j++) {
            fjac[i][j] += (g->u[mt][i-1]*g->u[mt][j-1])/A;
          }
        }
      }
    }
  }
  for (i=1;i<=n;i++){
    fvec[i] = fvec[i]/g->F_dependency + g->Fbuckle[i-1] - g->Stokes_translation*x[i];   /* Fbuckle[0,1,2] */
    for (j=1;j<=n;j++){
      fjac[i][j] = fjac[i][j]/(-1.0*g->F_dependency*g->k_on);
      if (j==i) {
        fjac[i][j] = fjac[i][j] - g->Stokes_translation;
      }
    }
  }
}

void function_MotorFV(double *x, int n, double *fvec, double **fjac, mtGlobal* g) {
  /* x[1],x[2],x[3] are velocity of Nuc for x,y,z axes, respectively */ 
  int i,j;
  for (i=1;i<=n;i++) {
    fvec[i] = g->Fbackward[i-1];
    for (j=1;j<=n;j++){
      fjac[i][j] = g->fjac_pull[i-1][j-1];
    }
  }
  //  for (j=1; j<=6; j++) {
  //    TRACE(("fjac[%d] %lf %lf %lf %lf %lf %lf\n",j,fjac[j][1],fjac[j][2],fjac[j][3],fjac[j][4],fjac[j][5],fjac[j][6]));
  //  }

  int mt;
  int qq;
  double dv, PullingF;
  for (mt=0; mt<g->N; mt++) {
    if (g->pulling_phase[mt]==1) {
      if (mt<g->NN) {qq=0;} else {qq=1;}
      dv = 0.0;
      for (j=0; j<3; j++) {
        dv += (x[(j+1)%3+4]*g->DVecNucCen[qq][(j+2)%3]-x[(j+2)%3+4]*g->DVecNucCen[qq][(j+1)%3]+x[j+1])*g->u[mt][j];
      }
      PullingF = g->MotorStallF*(1-dv/g->MotorMaxVel);
      for (i=1; i<=n; i++) {
        fvec[i] += g->NumberOfMotor[mt]*g->u[mt][i-1]*PullingF;
      }
    }
  }
  for (i=1;i<=(n/2);i++){
    fvec[i] -= g->Stokes_translation*x[i];   /* Fbackward[0,1,2] */
  }
  for (i=(n/2)+1;i<=n;i++){
    fvec[i] -= g->Stokes_rotation*x[i];
  }
}

/////////////////// the below function is identical to function_motorFV and needs to be modified ////////////////////
void function_laserMotorFV (double *x, int n, double *fvec, double **fjac, mtGlobal* g) {
  /* x[1],x[2],x[3] are velocity of Nuc for x,y,z axes, respectively */
  int i,j;
  for (i=1;i<=n;i++) {
    fvec[i]=0.0;
    for (j=1;j<=n;j++){
      fjac[i][j] = g->fjac_pull[i-1][j-1];
    }
  }
  int mt;
  int qq;
  double dv, PullingF;
  /* double dv, PullingF, dFdv; */
  /* double gradient = g->Stokes_translation; */
  for (mt=g->mt_start; mt<g->mt_end; mt++) {
    if (g->pulling_phase[mt]==1) {
      dv = g->u[mt][0]*x[1]+g->u[mt][1]*x[2]+g->u[mt][2]*x[3];
      if (mt<g->NN) {qq=0;} else {qq=1;}
      dv=0;
      for (j=0; j<3; j++) {
        dv += (x[(j+1)%3+4]*g->DVecNucCen[qq][(j+2)%3]-x[(j+2)%3+4]*g->DVecNucCen[qq][(j+1)%3]+x[j+1])*g->u[mt][j];
      }
      PullingF = g->MotorStallF*(1-dv/g->MotorMaxVel);
      for (i=1; i<=n; i++) {
        fvec[i] += g->NumberOfMotor[mt]*g->u[mt][i-1]*PullingF;
      }
    }
  }
  for (i=1;i<=(n/2);i++){
    fvec[i] = fvec[i]+g->Fbackward[i-1] - g->Stokes_translation*x[i];   /* Fbackward[0,1,2] */
  }
  for (i=(n/2)+1;i<=n;i++){
    fvec[i] = fvec[i] + g->Stokes_rotation*x[i];
  }
}

void function_anaMotorFV (double *x, int n, double *fvec, double **fjac, mtGlobal* g) {
  /* x[1],x[2],x[3] are velocity of Nuc for x,y,z axes, respectively */
  int i,j;
  for (i=1;i<=n;i++) {
    fvec[i]=0.0;
    for (j=1;j<=n;j++){
      fjac[i][j] = g->fjac_pull[i-1][j-1];
    }
  }
  int mt;
  int qq;
  double dv, PullingF;
  /* double dv, PullingF, dFdv; */
  /* double gradient = g->Stokes_translation; */
  for (mt=g->mt_start; mt<g->mt_end; mt++) {
    if (g->pulling_phase[mt]==1) {
      dv = g->u[mt][0]*x[1]+g->u[mt][1]*x[2]+g->u[mt][2]*x[3];
      if (mt<g->NN) {qq=0;} else {qq=1;}
      dv=0;
      for (j=0; j<3; j++) {
        dv += (x[(j+1)%3+4]*g->DVecNucCen[qq][(j+2)%3]-x[(j+2)%3+4]*g->DVecNucCen[qq][(j+1)%3]+x[j+1])*g->u[mt][j];
      }
      PullingF = g->MotorStallF*(1-dv/g->MotorMaxVel);
      for (i=1; i<=n; i++) {
        fvec[i] += g->NumberOfMotor[mt]*g->u[mt][i-1]*PullingF;
      }
    }
  }
  for (i=1;i<=(n/2);i++){
    fvec[i] = fvec[i]+g->Fbackward[i-1] - g->Stokes_translation*x[i];   /* Fbackward[0,1,2] */
  }
  for (i=(n/2)+1;i<=n;i++){
    fvec[i] = fvec[i] + g->Stokes_rotation*x[i];
  }
}
