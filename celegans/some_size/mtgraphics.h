/*
 * Author: Akatsuki Kimura <akkimura@nig.ac.jp>
 *         Akira Funahashi <funa@bio.keio.ac.jp>
 * Last modified: Thu, 25 Jul 2013 02:38:53 +0900
 */
#ifndef __mtgraphics__
#define __mtgraphics__

#include <X11/Xlib.h>
#include <X11/Xutil.h>

typedef struct _mtGraphics {
  int Scale1;
  Display *d;
  Colormap cm;
  Window w1,w2,w3,w4,w5,w6,w7,w8,w9;
  GC gc1, gc2, gc3, gc4, gc5, gc6, gc7, gc8, gc9, gc_clr;
  Pixmap pixmap;
  XSetWindowAttributes att;
  XColor c1,c2,c3,c4,c5,c6,cc;
} mtGraphics;

#endif
