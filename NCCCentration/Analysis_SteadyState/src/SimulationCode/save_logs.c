/*
 * filename: save_logs.c
 * this code is to store simulation results of the position of MTs and pronucleus in text files
 * Author: Akatsuki Kimura <akkimura@nig.ac.jp>
 *         Akira Funahashi <funa@bio.keio.ac.jp>
 * Last modified: Wed, 20 Jul 2022 01:46:48 +0900
 */
#include "mtsim.h"

void save_logs(int i, int p, int N, FILE* f_out1, FILE* f_out2, FILE* f_out3, FILE* f_out4, FILE* f_out10, FILE* data_for_3D, double Nuc[3], double PVecCen[2][3], double MT[][3]) {
  int m;
  // OUTPUT LOGSoutput per every single timepoint (20timepoints=1sec)
  if (i%20==0) {
    switch (p) {
      case 0:
        fprintf(f_out1,"%d %d %d %lf %lf %lf" ,p,(int)(i/20),i,Nuc[0]*pow(10,6),Nuc[1]*pow(10,6),Nuc[2]*pow(10,6));
        fprintf(f_out1," %lf %lf %lf %lf %lf %lf", PVecCen[0][0]*pow(10,6),PVecCen[0][1]*pow(10,6),PVecCen[0][2]*pow(10,6),PVecCen[1][0]*pow(10,6), PVecCen[1][1]*pow(10,6),PVecCen[1][2]*pow(10,6));
        fprintf(f_out1,"\n");
        break;
      case 1:
        fprintf(f_out2,"%d %d %d %lf %lf %lf" ,p,(int)(i/20),i,Nuc[0]*pow(10,6),Nuc[1]*pow(10,6),Nuc[2]*pow(10,6));
        fprintf(f_out2," %lf %lf %lf %lf %lf %lf", PVecCen[0][0]*pow(10,6),PVecCen[0][1]*pow(10,6),PVecCen[0][2]*pow(10,6),PVecCen[1][0]*pow(10,6), PVecCen[1][1]*pow(10,6),PVecCen[1][2]*pow(10,6));
        fprintf(f_out2,"\n");
        break;
      case 2:
        fprintf(f_out3,"%d %d %d %lf %lf %lf" ,p,(int)(i/20),i,Nuc[0]*pow(10,6),Nuc[1]*pow(10,6),Nuc[2]*pow(10,6));
        fprintf(f_out3," %lf %lf %lf %lf %lf %lf", PVecCen[0][0]*pow(10,6),PVecCen[0][1]*pow(10,6),PVecCen[0][2]*pow(10,6),PVecCen[1][0]*pow(10,6), PVecCen[1][1]*pow(10,6),PVecCen[1][2]*pow(10,6));
        fprintf(f_out3,"\n");
        break;
      case 3:
        fprintf(f_out4,"%d %d %d %lf %lf %lf" ,p,(int)(i/20),i,Nuc[0]*pow(10,6),Nuc[1]*pow(10,6),Nuc[2]*pow(10,6));
        fprintf(f_out4," %lf %lf %lf %lf %lf %lf", PVecCen[0][0]*pow(10,6),PVecCen[0][1]*pow(10,6),PVecCen[0][2]*pow(10,6),PVecCen[1][0]*pow(10,6), PVecCen[1][1]*pow(10,6),PVecCen[1][2]*pow(10,6));
        fprintf(f_out4,"\n");
        break;
      case 4:
        fprintf(f_out10,"%d %d %d %lf %lf %lf" ,p,(int)(i/20),i,Nuc[0]*pow(10,6),Nuc[1]*pow(10,6),Nuc[2]*pow(10,6));
        fprintf(f_out10," %lf %lf %lf %lf %lf %lf", PVecCen[0][0]*pow(10,6),PVecCen[0][1]*pow(10,6),PVecCen[0][2]*pow(10,6),PVecCen[1][0]*pow(10,6), PVecCen[1][1]*pow(10,6),PVecCen[1][2]*pow(10,6));
        fprintf(f_out10,"\n");
        break;
    }
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
