/*
 * filename: save_logs.c
 * this code is to store simulation results of the position of MTs and pronucleus in text files
 * Author: Akatsuki Kimura <akkimura@nig.ac.jp>
 *         Akira Funahashi <funa@bio.keio.ac.jp>
 * Last modified: Wed, 20 Jul 2022 01:49:08 +0900
 */
#include "mtsim.h"

void save_logs(int i, int p, int N, FILE* f_out1, FILE* f_out2, FILE* f_out3, FILE* f_out4, FILE* f_out5, FILE* f_out6, FILE* f_out7, FILE* f_out8, FILE* f_out9, FILE* f_out10, FILE* data_for_3D, double Nuc[3], double PVecCen[2][3], double MT[][3]) {
  int m;
  FILE* f_out;
  // OUTPUT LOGSoutput per every single timepoint (20timepoints=1sec)
  /* if (i%20==0) { */
  if (i%100==0) {
    switch (p) {
      case 0:
        f_out = f_out1;
        break;
      case 1:
        f_out = f_out2;
        break;
      case 2:
        f_out = f_out3;
        break;
      case 3:
        f_out = f_out4;
        break;
      case 4:
        f_out = f_out5;
        break;
      case 5:
        f_out = f_out6;
        break;
      case 6:
        f_out = f_out7;
        break;
      case 7:
        f_out = f_out8;
        break;
      case 8:
        f_out = f_out9;
        break;
      case 9:
        f_out = f_out10;
        break;
    }
    if(i==0){fprintf(f_out,"#p,i,time[min],Nuc[0]*10^6,Nuc_distance[per],PVecCen[0][0]*10^6,PVecCen[0][2]*10^6,PVecCen[1][0]*10^6,PVecCen[1][2]*10^6,angle[rad],angel[degree]\n");}
    fprintf(f_out,"%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf" ,p,i,i*dT/60,Nuc[0]*pow(10,6),(Rad-Nuc[0])/(2*Rad)*100,PVecCen[0][0]*pow(10,6),PVecCen[0][2]*pow(10,6),PVecCen[1][0]*pow(10,6),PVecCen[1][2]*pow(10,6),PI/2-fabs(atan(fabs(PVecCen[0][2]-PVecCen[1][2])/fabs(PVecCen[0][0]-PVecCen[1][0]))),(PI/2-fabs(atan(fabs(PVecCen[0][2]-PVecCen[1][2])/fabs(PVecCen[0][0]-PVecCen[1][0]))))*180/PI);
    fprintf(f_out,"\n");
  }

  // OUTPUT DATA FOR 3D animation /////// 
  if (p ==0) {
    fprintf(data_for_3D,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",Nuc[0]*pow(10,6),Nuc[1]*pow(10,6),Nuc[2]*pow(10,6),PVecCen[0][0]*pow(10,6),PVecCen[0][1]*pow(10,6),PVecCen[0][2]*pow(10,6),PVecCen[1][0]*pow(10,6), PVecCen[1][1]*pow(10,6),PVecCen[1][2]*pow(10,6));
    for (m=0; m<N; m++) {
      fprintf(data_for_3D,"%lf %lf %lf\n",MT[m][0]*pow(10,6),MT[m][1]*pow(10,6),MT[m][2]*pow(10,6));
    }
  }

  // OUTPUT 2D SNAPSHOTS ////
  /*	  
        if (i%1000 == 0){
        int llll = (int)(i/200);
        int lll = int(llll/100);
        int ll = int((llll-lll*100)/10);
        int l = int(llll-lll*100-ll*10);
        XStoreName(d,w1,"final");
        XFlush (d);

        char filefin[40];
        sprintf (filefin, "xwd -name \'final\' -out FIN%d%d%d -silent",lll,ll,l);
        system (filefin);

        char filefin2[40];
        sprintf (filefin2, "convert xwd:FIN%d%d%d FIN%d%d%d.jpg",lll,ll,l,lll,ll,l);
        system (filefin2);

        char filefin3[40];
        sprintf (filefin3, "convert xwd:FIN%d%d%d FIN%d%d%d.bmp",lll,ll,l,lll,ll,l);
        system (filefin3);
        }
        */
  //////////////////OUTPUT FINISHED //////////////////////////////
}
