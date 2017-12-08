/*
 * filename: usage.c
 * Author: Akira Funahashi <funa@bio.keio.ac.jp>
 * Last modified: Thu, 25 Jul 2013 02:58:12 +0900
 */
#include "mtsim.h"

void usage(char* myname) {
  printf("Usage : %s [option]\n", myname);
  printf(" -h   : Show this message\n");
  printf(" -c   : Check output with fixed random seed (0.0)\n");
  printf(" -v   : Verbose output\n");
  printf(" -m # : specify model (ex. -m 1 )\n");
  printf("        0: MT angle fixed\n");
  printf("        1: MT angle variable\n");
  printf(" -l # : metaspindle length (ex. -l 1 )\n");
  printf("        0: Rad vs metaspindle length\n");
  printf("        1: aspect ratio vs metaspindle length / Rad\n");
}
