/*
 * filename: draw_graphs.c
 * this code is to plot results of the simulation onto the graphs on the X window system
 * Author: Akatsuki Kimura <akkimura@nig.ac.jp>
 *         Akira Funahashi <funa@bio.keio.ac.jp>
 * Last modified: Thu, 25 Jul 2013 02:37:32 +0900
 */
#include "mtsim.h"

/* void draw_graphs(int i, mtGraphics *mtg, mtGlobal *g, double PVecCen[2][3], double MT[][3], double Nuc[3], double DistanceFromPP, double min, double max, double sum, int currentN, double OldVel, double NewVel) { */
void draw_graphs(int i, mtGraphics *mtg, mtGlobal *g, double PVecCen[2][3], double MT[][3], double Nuc[3], double Rad, double RadS, double MetaSpindle_L) {
  Display *d = mtg->d;
  int Scale1 = mtg->Scale1;
  // draw graph in window 1: distribution of MTs at the moment
  XFillRectangle(d,mtg->pixmap,mtg->gc_clr,0,0,WIN_WIDTH,WIN_HEIGHT);
  int m;
  for (m=0; m<g->N; m++) {
    if (m<g->NN){XDrawLine(d,mtg->pixmap,mtg->gc1,WIN_WIDTH/2+(int)(Scale1*PVecCen[0][0]*pow(10,6)),WIN_HEIGHT/2-(int)(Scale1*PVecCen[0][2]*pow(10,6)),WIN_WIDTH/2+(int)(Scale1*MT[m][0]*pow(10,6)),WIN_HEIGHT/2-(int)(Scale1*MT[m][2]*pow(10,6)));}
    else {XDrawLine(d,mtg->pixmap,mtg->gc1,WIN_WIDTH/2+(int)(Scale1*PVecCen[1][0]*pow(10,6)),WIN_HEIGHT/2-(int)(Scale1*PVecCen[1][2]*pow(10,6)),WIN_WIDTH/2+(int)(Scale1*MT[m][0]*pow(10,6)),WIN_HEIGHT/2-(int)(Scale1*MT[m][2]*pow(10,6)));}
  }
  /* if (i<laserST && i<anaST) { */
  /*   XDrawLine(d,mtg->pixmap,mtg->gc1,WIN_WIDTH/2+(int)(Scale1*PVecCen[0][0]*pow(10,6)),WIN_HEIGHT/2-(int)(Scale1*PVecCen[0][2]*pow(10,6)),WIN_WIDTH/2+(int)(Scale1*PVecCen[1][0]*pow(10,6)),WIN_HEIGHT/2-(int)(Scale1*PVecCen[1][2]*pow(10,6))); */
  /*   XDrawArc(d,mtg->pixmap,mtg->gc1,WIN_WIDTH/2+(int)((Scale1*(Nuc[0]-LL)*pow(10,6))),WIN_HEIGHT/2-(int)((Scale1*(Nuc[2]+LL)*pow(10,6))),(int)(Scale1*2*LL*pow(10,6)),(int)(Scale1*2*LL*pow(10,6)),0,360*64); */
  /* } */
  /* else if (i<laserST && i>=anaST) { */
  /*   XDrawLine(d,mtg->pixmap,mtg->gc1,WIN_WIDTH/2+(int)(Scale1*PVecCen[0][0]*pow(10,6)),WIN_HEIGHT/2-(int)(Scale1*PVecCen[0][2]*pow(10,6)),WIN_WIDTH/2+(int)(Scale1*PVecCen[1][0]*pow(10,6)),WIN_HEIGHT/2-(int)(Scale1*PVecCen[1][2]*pow(10,6))); */
  /* } */
  XDrawArc(d,mtg->pixmap,mtg->gc1,WIN_WIDTH/2-(int)(Scale1*Rad*pow(10,6)), WIN_HEIGHT/2-(int)(Scale1*RadS*pow(10,6)), (int)(2*Scale1*Rad*pow(10,6)), (int)(2*Scale1*RadS*pow(10,6)), 0, 360*64);
  XCopyArea(d,mtg->pixmap,mtg->w1,mtg->gc1,0,0,WIN_WIDTH,WIN_HEIGHT,0,0);


  /* // draw graph in window 2: path of the center of the pronucleus */
  /* XDrawPoint(d,mtg->w2,mtg->gc2,WIN_WIDTH/2+(int)(Scale1*Nuc[0]*pow(10,6)),WIN_HEIGHT/2-(int)(Scale1*Nuc[2]*pow(10,6))); */

  /* // draw graph in window 3: distance of the pronucleus from the center of the egg was plotted against time */
  /* XDrawPoint(d,mtg->w3,mtg->gc3,(int)(SC3*i),WIN_HEIGHT-(int)(2*Scale1*sqrt(Nuc[0]*Nuc[0]+Nuc[1]*Nuc[1]+Nuc[2]*Nuc[2])*pow(10,6))); */

  /* // draw graph in window 4: distance of the pronucleus from the posterior pole of the egg was plotted against time */
  /* XDrawPoint(d,mtg->w4,mtg->gc4,(int)(SC3*i),WIN_HEIGHT-(int)(Scale1*DistanceFromPP*pow(10,6))); */
  /* //Inset */
  /* /\**\/ */
  /* int insetcount=20; */
  /* // standard insetcount=20: 633<i<4000, 31.65<t(sec)<200 */
  /* // if insetcount=100: 3300<i<6667, 165<t(sec)<333.35 */
  /* // if insetcount=144: 4<t(min)<6.7 */
  /* // if insetcount=180: 5966<i<9333, 298.3<t(sec)<466.65 */
  /* int AMP=3; */
  /* double AMP_Y=2; */
  /* /\* if ((SC3*i*AMP>insetcount-1)&&(SC3*i*AMP<100+insetcount)&&(Scale1*(DistanceFromPP-LL)*pow(10,6)*AMP_Y<100)) *\/ */
  /* /\* {XDrawPoint(d,mtg->w4,mtg->gc4,(int)(14+SC3*i*AMP-insetcount),(int)(114-AMP_Y*Scale1*(DistanceFromPP-LL)*pow(10,6)));} *\/ */
  /* if ((SC3*i*AMP>insetcount-1)&&(SC3*i*AMP<100+insetcount)&&(Scale1*(DistanceFromPP-MetaSpindle_L)*pow(10,6)*AMP_Y<100)) */
  /* {XDrawPoint(d,mtg->w4,mtg->gc4,(int)(14+SC3*i*AMP-insetcount),(int)(114-AMP_Y*Scale1*(DistanceFromPP-MetaSpindle_L)*pow(10,6)));} */
  /* /\**\/ */
  /* /\* range of the inset can be alterd, for example... *\/ */
  /* //if ((SC3*i*4>189)&&(SC3*i*4<300)&&(Scale1*(DistanceFromPP-LL)*pow(10,6)*4<110)) */
  /* //{XDrawPoint(d,w4,gc4,(int)(14+SC3*i*4-190),(int)(124-4*Scale1*(DistanceFromPP-LL)*pow(10,6)));} */


  /* // draw graph in window 5: maximum, minimum, and mean lengths of the MTs were plotted against time */
  /* XDrawPoint(d,mtg->w5,mtg->gc5,(int)(SC3*i),WIN_HEIGHT-(int)(Scale1*pow(10,6)*sum/g->N)); */
  /* XDrawPoint(d,mtg->w5,mtg->gc5,(int)(SC3*i),WIN_HEIGHT-(int)(Scale1*pow(10,6)*max)); */
  /* XDrawPoint(d,mtg->w5,mtg->gc5,(int)(SC3*i),WIN_HEIGHT-(int)(Scale1*pow(10,6)*min)); */

  /* // draw graph in window 6: length of a single MT was plotted against time */
  /* XDrawPoint(d,mtg->w6,mtg->gc6,(int)(SC3*i),WIN_HEIGHT-(int)(2*Scale1*pow(10,6)*g->L[0])); */


  /* //draw graph in window 7: the number of MTs was plotted against time */
  /* XDrawPoint(d,mtg->w7,mtg->gc7,(int)(SC3*i),WIN_HEIGHT-(int)((Scale1*currentN)/4.0)); */

  /* //draw graph in window 8: velocity of the pronucleus was plotted against time */
  /* if ((i%UNITTIME==0)&&(i!=0)){ */
  /*   XSetLineAttributes (d,mtg->gc8,1,LineSolid, CapButt,JoinMiter); */
  /*   XDrawLine(d, mtg->w8,mtg->gc8,(int)(SC3*(i-UNITTIME)),WIN_HEIGHT-(int)(OldVel*pow(10,6)*1000), (int)(SC3*i),WIN_HEIGHT-(int)(NewVel*pow(10,6)*1000)); */
  /* } */
  return;
}
