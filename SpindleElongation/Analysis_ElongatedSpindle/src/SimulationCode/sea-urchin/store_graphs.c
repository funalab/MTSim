/*
 * filename: store_graphs.c
 * this code is to save graphs
 * Author: Akatsuki Kimura <akkimura@nig.ac.jp>
 *         Akira Funahashi <funa@bio.keio.ac.jp>
 * Last modified: Wed, 20 Jul 2022 01:48:36 +0900
 */
#include "mtsim.h"

void store_graphs(void) {
  char filefin[40];
  sprintf (filefin, "xwd -name \'fin\' -out outFINf");
  system (filefin);
  char filefin2[40];
  sprintf (filefin2, "convert xwd:outFINf outFINf.jpg");
  system (filefin2);

  char filepat[40];
  sprintf (filepat, "xwd -name \'path\' -out outPAT");
  system (filepat);
  char filepat2[40];
  sprintf (filepat2, "convert xwd:outPAT outPAT.jpg");
  system (filepat2);

  char filedis[40];
  sprintf (filedis, "xwd -name \'distance\' -out outDIS");
  system (filedis);
  char filedis2[40];
  sprintf (filedis2, "convert xwd:outDIS outDIS.jpg");
  system (filedis2);

  char filePdis[40];
  sprintf (filePdis, "xwd -name \'poledistance\' -out outPDIS");
  system (filePdis);
  char filePdis2[40];
  sprintf (filePdis2, "convert xwd:outPDIS outPDIS.jpg");
  system (filePdis2);

  char filelen[40];
  sprintf (filelen, "xwd -name \'length\' -out outLEN");
  system (filelen);
  char filelen2[40];
  sprintf (filelen2, "convert xwd:outLEN outLEN.jpg");
  system (filelen2);

  char filesin[40];
  sprintf (filesin, "xwd -name \'single\' -out outSIN");
  system (filesin);
  char filesin2[40];
  sprintf (filesin2, "convert xwd:outSIN outSIN.jpg");
  system (filesin2);

  char fileMTnum[40];
  sprintf (fileMTnum, "xwd -name \'MTnum\' -out outMTnum");
  system (fileMTnum);
  char fileMTnum2[40];
  sprintf (fileMTnum2, "convert xwd:outMTnum outMTnum.jpg");
  system (fileMTnum2);

  char fileVec[40];
  sprintf (fileVec, "xwd -name \'vector\' -out outVEC");
  system (fileVec);
  char fileVec2[40];
  sprintf (fileVec2, "convert xwd:outVEC outVEC.jpg");
  system (fileVec2);
}
