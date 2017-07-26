/* 
 * filename: color_setting.c
 * this code is to discriminate the results from distinct parameter sets with
 * different colors: case0=black, case1=violet(pink), case2=blue, case3=green, case4=orenge, and case5=red 
 * Author: Akatsuki Kimura <akkimura@nig.ac.jp>
 *         Akira Funahashi <funa@bio.keio.ac.jp>
 * Last modified: Thu, 25 Jul 2013 02:37:07 +0900
 */
#include "mtsim.h"

void color_setting(mtGraphics *mtg, int p) {
  /* color settings */
  switch(p) {
    case 1:
      XSetForeground(mtg->d, mtg->gc2, mtg->c1.pixel);
      XSetForeground(mtg->d, mtg->gc3, mtg->c1.pixel);
      XSetForeground(mtg->d, mtg->gc4, mtg->c1.pixel);
      XSetForeground(mtg->d, mtg->gc5, mtg->c1.pixel);
      XSetForeground(mtg->d, mtg->gc6, mtg->c1.pixel);
      XSetForeground(mtg->d, mtg->gc7, mtg->c1.pixel);
      XSetForeground(mtg->d, mtg->gc8, mtg->c1.pixel);
      break;
    case 2:
      XSetForeground(mtg->d, mtg->gc2, BlackPixel(mtg->d,0));
      XSetForeground(mtg->d, mtg->gc3, BlackPixel(mtg->d,0));
      XSetForeground(mtg->d, mtg->gc4, BlackPixel(mtg->d,0));
      XSetForeground(mtg->d, mtg->gc5, BlackPixel(mtg->d,0));
      XSetForeground(mtg->d, mtg->gc6, BlackPixel(mtg->d,0));
      XSetForeground(mtg->d, mtg->gc7, BlackPixel(mtg->d,0));
      XSetForeground(mtg->d, mtg->gc8, BlackPixel(mtg->d,0));
      break;
    case 3:
      XSetForeground(mtg->d, mtg->gc2, mtg->c2.pixel);
      XSetForeground(mtg->d, mtg->gc3, mtg->c2.pixel);
      XSetForeground(mtg->d, mtg->gc4, mtg->c2.pixel);
      XSetForeground(mtg->d, mtg->gc5, mtg->c2.pixel);
      XSetForeground(mtg->d, mtg->gc6, mtg->c2.pixel);
      XSetForeground(mtg->d, mtg->gc7, mtg->c2.pixel);
      XSetForeground(mtg->d, mtg->gc8, mtg->c2.pixel);
      break;
    case 4:
      XSetForeground(mtg->d, mtg->gc2, mtg->c3.pixel);
      XSetForeground(mtg->d, mtg->gc3, mtg->c3.pixel);
      XSetForeground(mtg->d, mtg->gc4, mtg->c3.pixel);
      XSetForeground(mtg->d, mtg->gc5, mtg->c3.pixel);
      XSetForeground(mtg->d, mtg->gc6, mtg->c3.pixel);
      XSetForeground(mtg->d, mtg->gc7, mtg->c3.pixel);
      XSetForeground(mtg->d, mtg->gc8, mtg->c3.pixel);
      break;
    case 0:
      XSetForeground(mtg->d, mtg->gc2, mtg->c4.pixel);
      XSetForeground(mtg->d, mtg->gc3, mtg->c4.pixel);
      XSetForeground(mtg->d, mtg->gc4, mtg->c4.pixel);
      XSetForeground(mtg->d, mtg->gc5, mtg->c4.pixel);
      XSetForeground(mtg->d, mtg->gc6, mtg->c4.pixel);
      XSetForeground(mtg->d, mtg->gc7, mtg->c4.pixel);
      XSetForeground(mtg->d, mtg->gc8, mtg->c4.pixel);
      break;
    case 5:
      XSetForeground(mtg->d, mtg->gc2, mtg->c5.pixel);
      XSetForeground(mtg->d, mtg->gc3, mtg->c5.pixel);
      XSetForeground(mtg->d, mtg->gc4, mtg->c5.pixel);
      XSetForeground(mtg->d, mtg->gc5, mtg->c5.pixel);
      XSetForeground(mtg->d, mtg->gc6, mtg->c5.pixel);
      XSetForeground(mtg->d, mtg->gc7, mtg->c5.pixel);
      XSetForeground(mtg->d, mtg->gc8, mtg->c5.pixel);
      break;
  }
}
