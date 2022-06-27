/*
 * filename: display_setting.c
 * this code is to draw graphs in the X window system
 * Author: Akatsuki Kimura <akkimura@nig.ac.jp>
 *         Akira Funahashi <funa@bio.keio.ac.jp>
 * Last modified: Thu, 25 Jul 2013 02:37:27 +0900
 */
#include "mtsim.h"

void display_setting(mtGraphics *mtg) {
  // display settings
  mtg->d = XOpenDisplay(NULL);
  if (mtg->d == NULL) {
    printf("XOpenDisplay failed.\n");
    exit(1);
  }
  Display *d = mtg->d;
  mtg->cm = DefaultColormap(d,0);
  Colormap cm = mtg->cm;
  XAllocNamedColor(d,cm,"violet", &mtg->c1, &mtg->cc);
  XAllocNamedColor(d,cm,"blue", &mtg->c2, &mtg->cc);
  XAllocNamedColor(d,cm,"green", &mtg->c3, &mtg->cc);
  XAllocNamedColor(d,cm,"orange", &mtg->c4, &mtg->cc);
  XAllocNamedColor(d,cm,"red", &mtg->c5, &mtg->cc);
  mtg->w1 = XCreateSimpleWindow(d, RootWindow(d,0), 5, 5, WIN_WIDTH, WIN_HEIGHT, 2, BlackPixel(d,0), WhitePixel(d,0));
  mtg->w2 = XCreateSimpleWindow(d, RootWindow(d,0), 10+WIN_WIDTH, 5, WIN_WIDTH, WIN_HEIGHT, 2, BlackPixel(d,0), WhitePixel(d,0));
  mtg->w3 = XCreateSimpleWindow(d, RootWindow(d,0), 5, 10+WIN_HEIGHT, WIN_WIDTH, WIN_HEIGHT, 2, BlackPixel(d,0), WhitePixel(d,0));
  mtg->w4 = XCreateSimpleWindow(d, RootWindow(d,0), 10+WIN_WIDTH, 10+WIN_HEIGHT, WIN_WIDTH, WIN_HEIGHT, 2, BlackPixel(d,0), WhitePixel(d,0));
  mtg->w5 = XCreateSimpleWindow(d, RootWindow(d,0), 5, 15+2*WIN_HEIGHT, WIN_WIDTH, WIN_HEIGHT, 2, BlackPixel(d,0), WhitePixel(d,0));
  mtg->w6 = XCreateSimpleWindow(d, RootWindow(d,0), 10+WIN_WIDTH, 15+2*WIN_HEIGHT, WIN_WIDTH, WIN_HEIGHT, 2, BlackPixel(d,0), WhitePixel(d,0));
  mtg->w7 = XCreateSimpleWindow(d, RootWindow(d,0), 15+2*+WIN_WIDTH, 5, WIN_WIDTH, WIN_HEIGHT, 2, BlackPixel(d,0), WhitePixel(d,0));
  mtg->w8 = XCreateSimpleWindow(d, RootWindow(d,0), 15+2*+WIN_WIDTH, 10+WIN_HEIGHT, WIN_WIDTH, WIN_HEIGHT, 2, BlackPixel(d,0), WhitePixel(d,0));
  mtg->w9 = XCreateSimpleWindow(d, RootWindow(d,0), 15+2*+WIN_WIDTH, 15+2*WIN_HEIGHT, WIN_WIDTH, WIN_HEIGHT, 2, BlackPixel(d,0), WhitePixel(d,0));
  mtg->pixmap = XCreatePixmap(d, mtg->w1, WIN_WIDTH, WIN_HEIGHT, DefaultDepth(d,0));
  mtg->att.override_redirect = 1;
  XChangeWindowAttributes(d,mtg->w1,CWOverrideRedirect,&mtg->att);
  XChangeWindowAttributes(d,mtg->w2,CWOverrideRedirect,&mtg->att);
  XChangeWindowAttributes(d,mtg->w3,CWOverrideRedirect,&mtg->att);
  XChangeWindowAttributes(d,mtg->w4,CWOverrideRedirect,&mtg->att);
  XChangeWindowAttributes(d,mtg->w5,CWOverrideRedirect,&mtg->att);
  XChangeWindowAttributes(d,mtg->w6,CWOverrideRedirect,&mtg->att);
  XChangeWindowAttributes(d,mtg->w7,CWOverrideRedirect,&mtg->att);
  XChangeWindowAttributes(d,mtg->w8,CWOverrideRedirect,&mtg->att);
  XChangeWindowAttributes(d,mtg->w9,CWOverrideRedirect,&mtg->att);
  XMapWindow(d,mtg->w1);
  mtg->gc1 = XCreateGC(d,mtg->w1,0,0);
  XMapWindow(d,mtg->w2);
  mtg->gc2 = XCreateGC(d,mtg->w2,0,0);
  XMapWindow(d,mtg->w3);
  mtg->gc3 = XCreateGC(d,mtg->w3,0,0);
  XMapWindow(d,mtg->w4);
  mtg->gc4 = XCreateGC(d,mtg->w4,0,0);
  XMapWindow(d,mtg->w5);
  mtg->gc5 = XCreateGC(d,mtg->w5,0,0);
  XMapWindow(d,mtg->w6);
  mtg->gc6 = XCreateGC(d,mtg->w6,0,0);
  XMapWindow(d,mtg->w7);
  mtg->gc7 = XCreateGC(d,mtg->w7,0,0);
  XMapWindow(d,mtg->w8);
  mtg->gc8 = XCreateGC(d,mtg->w8,0,0);
  XMapWindow(d,mtg->w9);
  mtg->gc9 = XCreateGC(d,mtg->w9,0,0);
  mtg->gc_clr = XCreateGC(d,mtg->w1,0,0);
  XSetBackground(d,mtg->gc_clr,WhitePixel(d,0));
  XSetForeground(d,mtg->gc_clr,WhitePixel(d,0));
  int Scale1 = 4; /* for visualization scale setting 100/Rad[um] */
  mtg->Scale1 = Scale1;

  /******* window2 ********/
  XDrawArc(d,mtg->w2,mtg->gc2,WIN_WIDTH/2-(int)(Scale1*Rad*pow(10,6)), WIN_HEIGHT/2-(int)(Scale1*RadS*pow(10,6)), (int)(2*Scale1*Rad*pow(10,6)), (int)(2*Scale1*RadS*pow(10,6)), 0, 360*64);


  /******* window3 ********/
  XDrawLine(d,mtg->w3,mtg->gc3,0,WIN_HEIGHT-(int)(2*Scale1*25),5,WIN_HEIGHT-(int)(2*Scale1*25));
  XDrawLine(d,mtg->w3,mtg->gc3,0,WIN_HEIGHT-(int)(2*Scale1*20),5,WIN_HEIGHT-(int)(2*Scale1*20)); 
  XDrawLine(d,mtg->w3,mtg->gc3,0,WIN_HEIGHT-(int)(2*Scale1*10),5,WIN_HEIGHT-(int)(2*Scale1*10));
  XDrawLine(d,mtg->w3,mtg->gc3,(int)(SC3*anaST),WIN_HEIGHT,(int)(SC3*anaST),WIN_HEIGHT-5);
  XDrawLine(d,mtg->w3,mtg->gc3,(int)(SC3*anaST*1/4),WIN_HEIGHT,(int)(SC3*anaST*1/4),WIN_HEIGHT-5);
  XDrawLine(d,mtg->w3,mtg->gc3,(int)(SC3*anaST*2/4),WIN_HEIGHT,(int)(SC3*anaST*2/4),WIN_HEIGHT-5);
  XDrawLine(d,mtg->w3,mtg->gc3,(int)(SC3*anaST*3/4),WIN_HEIGHT,(int)(SC3*anaST*3/4),WIN_HEIGHT-5);


  /********** window4 **********/
  XDrawLine(d,mtg->w4,mtg->gc4,0,WIN_HEIGHT-(int)(Scale1*50),5,WIN_HEIGHT-(int)(Scale1*50));
  XDrawLine(d,mtg->w4,mtg->gc4,0,WIN_HEIGHT-(int)(Scale1*40),5,WIN_HEIGHT-(int)(Scale1*40)); 
  XDrawLine(d,mtg->w4,mtg->gc4,0,WIN_HEIGHT-(int)(Scale1*30),5,WIN_HEIGHT-(int)(Scale1*30)); 
  XDrawLine(d,mtg->w4,mtg->gc4,0,WIN_HEIGHT-(int)(Scale1*20),5,WIN_HEIGHT-(int)(Scale1*20)); 
  XDrawLine(d,mtg->w4,mtg->gc4,0,WIN_HEIGHT-(int)(Scale1*10),5,WIN_HEIGHT-(int)(Scale1*10)); 
  XDrawLine(d,mtg->w4,mtg->gc4,(int)(SC3*anaST),WIN_HEIGHT,(int)(SC3*anaST),WIN_HEIGHT-5);
  XDrawLine(d,mtg->w4,mtg->gc4,(int)(SC3*anaST/2),WIN_HEIGHT,(int)(SC3*anaST/2),WIN_HEIGHT-5);
  XDrawLine(d,mtg->w4,mtg->gc4,(int)(SC3*anaST*3/4),WIN_HEIGHT,(int)(SC3*anaST*3/4),WIN_HEIGHT-5);
  XDrawLine(d,mtg->w4,mtg->gc4,(int)(SC3*anaST/4),WIN_HEIGHT,(int)(SC3*anaST/4),WIN_HEIGHT-5);
  /* inset */
  /**/
  XDrawRectangle(d,mtg->w4,mtg->gc4,14,14,100,100);
  //XDrawLine(d,w4,gc4,14,34,19,34);
  XDrawLine(d,mtg->w4,mtg->gc4,14,64,19,64); 
  //XDrawLine(d,w4,gc4,14,124,19,124); 
  XDrawLine(d,mtg->w4,mtg->gc4,64,109,64,114);
  XDrawLine(d,mtg->w4,mtg->gc4,114,109,114,114);
  /**/
  XSetLineAttributes(d,mtg->gc4,1,LineOnOffDash, CapButt, JoinMiter);
  XDrawLine(d,mtg->w4,mtg->gc4,0,WIN_HEIGHT-(int)(Scale1*25),(int)(SC3*anaST),WIN_HEIGHT-(int)(Scale1*25));
  //XDrawLine(d,w4,gc4,14,124,124,124);
  XSetLineAttributes(d,mtg->gc4,1,LineSolid, CapButt,JoinMiter);


  /******* window5 ********/
  XDrawLine(d,mtg->w5,mtg->gc5,0,WIN_HEIGHT-(int)(Scale1*50),5,WIN_HEIGHT-(int)(Scale1*50));   
  XDrawLine(d,mtg->w5,mtg->gc5,0,WIN_HEIGHT-(int)(Scale1*40),5,WIN_HEIGHT-(int)(Scale1*40));   
  XDrawLine(d,mtg->w5,mtg->gc5,0,WIN_HEIGHT-(int)(Scale1*30),5,WIN_HEIGHT-(int)(Scale1*30));   
  XDrawLine(d,mtg->w5,mtg->gc5,0,WIN_HEIGHT-(int)(Scale1*20),5,WIN_HEIGHT-(int)(Scale1*20));   
  XDrawLine(d,mtg->w5,mtg->gc5,0,WIN_HEIGHT-(int)(Scale1*10),5,WIN_HEIGHT-(int)(Scale1*10)); 
  XDrawLine(d,mtg->w5,mtg->gc5,(int)(SC3*anaST),WIN_HEIGHT,(int)(SC3*anaST),WIN_HEIGHT-5);
  XDrawLine(d,mtg->w5,mtg->gc5,(int)(SC3*anaST/4),WIN_HEIGHT,(int)(SC3*anaST/4),WIN_HEIGHT-5);
  XDrawLine(d,mtg->w5,mtg->gc5,(int)(SC3*anaST/2),WIN_HEIGHT,(int)(SC3*anaST/2),WIN_HEIGHT-5);
  XDrawLine(d,mtg->w5,mtg->gc5,(int)(SC3*anaST*3/4),WIN_HEIGHT,(int)(SC3*anaST*3/4),WIN_HEIGHT-5);


  /******* window6 ********/
  XDrawLine(d,mtg->w6,mtg->gc6,0,WIN_HEIGHT-(int)(Scale1*50),5,WIN_HEIGHT-(int)(Scale1*50));   
  XDrawLine(d,mtg->w6,mtg->gc6,0,WIN_HEIGHT-(int)(Scale1*40),5,WIN_HEIGHT-(int)(Scale1*40));   
  XDrawLine(d,mtg->w6,mtg->gc6,0,WIN_HEIGHT-(int)(Scale1*30),5,WIN_HEIGHT-(int)(Scale1*30));   
  XDrawLine(d,mtg->w6,mtg->gc6,0,WIN_HEIGHT-(int)(Scale1*20),5,WIN_HEIGHT-(int)(Scale1*20));   
  XDrawLine(d,mtg->w6,mtg->gc6,0,WIN_HEIGHT-(int)(Scale1*10),5,WIN_HEIGHT-(int)(Scale1*10)); 
  XDrawLine(d,mtg->w6,mtg->gc6,(int)(SC3*anaST),WIN_HEIGHT,(int)(SC3*anaST),WIN_HEIGHT-5);
  XDrawLine(d,mtg->w6,mtg->gc6,(int)(SC3*anaST/4),WIN_HEIGHT,(int)(SC3*anaST/4),WIN_HEIGHT-5);
  XDrawLine(d,mtg->w6,mtg->gc6,(int)(SC3*anaST/2),WIN_HEIGHT,(int)(SC3*anaST/2),WIN_HEIGHT-5);
  XDrawLine(d,mtg->w6,mtg->gc6,(int)(SC3*anaST*3/4),WIN_HEIGHT,(int)(SC3*anaST*3/4),WIN_HEIGHT-5);


  /********** window7 **********/
  XDrawLine(d,mtg->w7,mtg->gc7,0,WIN_HEIGHT-(int)(Scale1*50),5,WIN_HEIGHT-(int)(Scale1*50)); 
  XDrawLine(d,mtg->w7,mtg->gc7,0,WIN_HEIGHT-(int)(Scale1*37.5),5,WIN_HEIGHT-(int)(Scale1*37.5)); 
  XDrawLine(d,mtg->w7,mtg->gc7,0,WIN_HEIGHT-(int)(Scale1*25),5,WIN_HEIGHT-(int)(Scale1*25)); 
  XDrawLine(d,mtg->w7,mtg->gc7,0,WIN_HEIGHT-(int)(Scale1*12.5),5,WIN_HEIGHT-(int)(Scale1*12.5)); 
  XDrawLine(d,mtg->w7,mtg->gc7,(int)(SC3*anaST),WIN_HEIGHT,(int)(SC3*anaST),WIN_HEIGHT-5);
  XDrawLine(d,mtg->w7,mtg->gc7,(int)(SC3*anaST/2),WIN_HEIGHT,(int)(SC3*anaST/2),WIN_HEIGHT-5);
  XDrawLine(d,mtg->w7,mtg->gc7,(int)(SC3*anaST*3/4),WIN_HEIGHT,(int)(SC3*anaST*3/4),WIN_HEIGHT-5);
  XDrawLine(d,mtg->w7,mtg->gc7,(int)(SC3*anaST/4),WIN_HEIGHT,(int)(SC3*anaST/4),WIN_HEIGHT-5);
  XSetLineAttributes (d,mtg->gc7,1,LineOnOffDash, CapButt, JoinMiter);
  XDrawLine(d,mtg->w7,mtg->gc7,0,WIN_HEIGHT-(int)(Scale1*25),(int)(SC3*anaST),WIN_HEIGHT-(int)(Scale1*25));
  XDrawLine(d,mtg->w7,mtg->gc7,0,WIN_HEIGHT-(int)(Scale1*50),(int)(SC3*anaST),WIN_HEIGHT-(int)(Scale1*50));
  XSetLineAttributes (d,mtg->gc7,1,LineSolid, CapButt,JoinMiter);


  /********** window8 **********/
  XDrawLine(d,mtg->w8,mtg->gc8,0,WIN_HEIGHT-(int)(Scale1*50),5,WIN_HEIGHT-(int)(Scale1*50));
  XDrawLine(d,mtg->w8,mtg->gc8,0,WIN_HEIGHT-(int)(Scale1*40),5,WIN_HEIGHT-(int)(Scale1*40)); 
  XDrawLine(d,mtg->w8,mtg->gc8,0,WIN_HEIGHT-(int)(Scale1*30),5,WIN_HEIGHT-(int)(Scale1*30)); 
  XDrawLine(d,mtg->w8,mtg->gc8,0,WIN_HEIGHT-(int)(Scale1*20),5,WIN_HEIGHT-(int)(Scale1*20)); 
  XDrawLine(d,mtg->w8,mtg->gc8,0,WIN_HEIGHT-(int)(Scale1*10),5,WIN_HEIGHT-(int)(Scale1*10)); 
  XDrawLine(d,mtg->w8,mtg->gc8,(int)(SC3*anaST),WIN_HEIGHT,(int)(SC3*anaST),WIN_HEIGHT-5);
  XDrawLine(d,mtg->w8,mtg->gc8,(int)(SC3*anaST/2),WIN_HEIGHT,(int)(SC3*anaST/2),WIN_HEIGHT-5);
  XDrawLine(d,mtg->w8,mtg->gc8,(int)(SC3*anaST*3/4),WIN_HEIGHT,(int)(SC3*anaST*3/4),WIN_HEIGHT-5);
  XDrawLine(d,mtg->w8,mtg->gc8,(int)(SC3*anaST/4),WIN_HEIGHT,(int)(SC3*anaST/4),WIN_HEIGHT-5);

  /********** window9 **********/
  XDrawLine(d,mtg->w9,mtg->gc9,WIN_WIDTH/2,0,WIN_WIDTH/2,WIN_HEIGHT);
  XDrawLine(d,mtg->w9,mtg->gc9,0,WIN_HEIGHT/2,WIN_WIDTH,WIN_HEIGHT/2);

  return;
}
