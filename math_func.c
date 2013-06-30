/*
 * filename: function8.c
 * this code is a collection of mathematical functions used in the simulation
 * Last modified: Mon, 01 Jul 2013 02:50:07 +0900
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI (3.141592653590)

// Length of a vector: L is the length of a[3]
double Length(double a[3]) {
  double L;
  L = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
  return L;
}


// Unit vector: vector2 is the unit vector of vector1
void UnitVector(double vector1[3], double vector2[3]) {
  int i;
  double L;
  L = Length(vector1);
  for (i=0; i<3; i++){
    vector2[i] = vector1[i] / L;
  }
  return;
}

void MakeRotationMatrix(double RotationMatrix[3][3], double RotationVector[3], double dt) {
  double L;
  L = Length(RotationVector);
  if ((dt>1.0e-4)&&(L>1.0e-6)) {
    double UnitRotVec[3];
    UnitVector(RotationVector, UnitRotVec);
    double AA, BB, CC;
    AA = cos((-1.0)*L*dt);
    BB = sin((-1.0)*L*dt);
    CC = 1.0-AA;
    RotationMatrix[0][0] = AA + CC*UnitRotVec[0]*UnitRotVec[0];
    RotationMatrix[1][1] = AA + CC*UnitRotVec[1]*UnitRotVec[1];
    RotationMatrix[2][2] = AA + CC*UnitRotVec[2]*UnitRotVec[2];
    RotationMatrix[0][1] = CC*UnitRotVec[0]*UnitRotVec[1] - BB*UnitRotVec[2];
    RotationMatrix[0][2] = CC*UnitRotVec[0]*UnitRotVec[2] + BB*UnitRotVec[1];
    RotationMatrix[1][0] = CC*UnitRotVec[0]*UnitRotVec[1] + BB*UnitRotVec[2];
    RotationMatrix[1][2] = CC*UnitRotVec[1]*UnitRotVec[2] - BB*UnitRotVec[0];
    RotationMatrix[2][0] = CC*UnitRotVec[0]*UnitRotVec[2] - BB*UnitRotVec[1];
    RotationMatrix[2][1] = CC*UnitRotVec[1]*UnitRotVec[2] + BB*UnitRotVec[0];
  } else {
    int ii, jj;
    for (ii=0; ii<3; ii++) {
      for (jj=0; jj<3; jj++) {
        if (ii==jj) {
          RotationMatrix[ii][jj]=1.0;
        } else {
          RotationMatrix[ii][jj]=0.0;
        }
      }
    }
  }
}

// Addition of vectors: vector1 + vector2 = result
void AddVector(double vector1[3], double vector2[3], double result[3]) {
  int i;
  for (i=0; i<3; i++){
    result[i] = vector1[i]+vector2[i];
  }
  return;
}

// Subtraction of vectors: vector1 - vector2 = result
void SubVector(double vector1[3], double vector2[3], double result[3]) {
  int i;
  for (i=0; i<3; i++){
    result[i] = vector1[i]-vector2[i];
  }
  return;
}

// Inner Product of vectors: vector1*vector2 = result
double InnProdVector(double vector1[3], double vector2[3]) {
  double result;
  result = vector1[0]*vector2[0]+vector1[1]*vector2[1]+vector1[2]*vector2[2];
  return result;
}


// Outer Product of vectors: vector1 X vector2 = result
void OutProdVector(double vector1[3], double vector2[3], double result[3]) {
  result[0]=vector1[1]*vector2[2]-vector1[2]*vector2[1];
  result[1]=vector1[2]*vector2[0]-vector1[0]*vector2[2];
  result[2]=vector1[0]*vector2[1]-vector1[1]*vector2[0];
  return;
}

// solution of quadratic equation: a[0]x2+a[1]x+a[2]=0
void QuadEqu2(double a[3], double b[2]) {
  if (a[0]==0){
    if (a[1]==0){printf("no solution!\n");}
    else {b[0]=-1*a[2]/a[1];b[1]=-1*a[2]/a[1];}
  }
  else {
    double q, qq;
    q = a[1]*a[1]-4*a[0]*a[2];
    if (q<0){printf("no solution!\n");}
    else {
      if (a[1]>=0){qq=-1*(a[1]+sqrt(q))/2;}
      else {qq=-1*(a[1]-sqrt(q))/2;}
      if (qq/a[0]>=a[2]/qq)
	{b[0]=qq/a[0];b[1]=a[2]/qq;}
      else {b[1]=qq/a[0];b[0]=a[2]/qq;}
    }
  }
  return;
}


// calculation of new axes
// center of nucleus, centrosome,
void ProductJacVec(double aa_out[3], double bb_jac[3][3], double cc_in[3]) {
  int i, j;
  //  TRACE(("in=(%lf, %lf, %lf)\n",cc_in[0],cc_in[1],cc_in[2]));
  for (i=0; i<3; i++) {
    //    TRACE(("in[%d]=%lf in=%d out[%d]=%lf out=%d\n",i,cc_in[i],cc_in,i,aa_out[i],aa_out));
    aa_out[i] = 0.0;
    //    TRACE(("in[%d]=%lf in=%d out[%d]=%lf out=%d\n",i,cc_in[i],cc_in,i,aa_out[i],aa_out));
    for (j=0; j<3; j++) {
      //      TRACE(("jac[%d][%d]=%lf,in[%d]=%lf\n",i,j,bb_jac[i][j],j,cc_in[j]));
      aa_out[i] += bb_jac[i][j]*cc_in[j];
      //      TRACE(("out[%d]=%lf,jac[%d][%d]=%lf,in[%d]=%lf\n",i,aa_out[i],i,j,bb_jac[i][j],j,cc_in[j]));
    }
  }
  return;
}

void VectorRotation(double out[3], double in[3], double rotAx[3], double rotStd[3], double degree) {
  double jac[3][3];
  double A=cos(degree);
  double B=sin(degree);
  double C=1.0-A;
  int k, j;
  //  TRACE(("in=(%lf, %lf, %lf)\n",in[0],in[1],in[2]));
  for (k=0; k<3; k++) {
    for (j=0; j<3; j++) {
      if (k==j) {jac[k][j]=A+C*rotAx[k]*rotAx[k];}
      else {
	if (((k<j)&&((k+j)%2==0))||((k>j)&&((k+j)%2==1))) {
	  jac[k][j]=C*rotAx[k]*rotAx[j]+B*rotAx[3-k-j];
	} else {
	  jac[k][j]=C*rotAx[k]*rotAx[j]-B*rotAx[3-k-j];
	}
      }
      //      TRACE(("jac[%d][%d]=%lf\n",k,j,jac[k][j]));
    }
  }
  ProductJacVec(out, jac, in);
  //  TRACE(("temp out=(%lf, %lf, %lf)\n",out[0],out[1],out[2]));
  for (k=0; k<3; k++) {
    for (j=0; j<3; j++) {
      jac[k][j] = -jac[k][j];
      if (k==j) {jac[k][j]+=1;}
    }
  }
  double temp[3];
  ProductJacVec(temp, jac, rotStd);
  AddVector(out, temp, out);
  //  TRACE(("final out=(%lf, %lf, %lf)\n",out[0],out[1],out[2]));
  return;
}

double Poisson(double np, int motornumber) {
  double result;
  int i, factorial=1;
  for (i=1; i<=motornumber; i++) {
    factorial*=i;
  }
  result = exp(-np)*pow(np,motornumber)/factorial;
  return result;
}
