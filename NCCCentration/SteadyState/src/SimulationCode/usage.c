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
  /* printf(" -t # : specify simulation time (ex. -t 100 )\n"); */
  /* printf(" -s # : specify simulation step (ex. -s 100 )\n"); */
  /* printf(" -d # : specify simulation delta (ex. -d 0.01 [default:1/4096])\n"); */
  printf(" -s # : specify strain (ex. -s 5 )\n");
  printf("        0: WT\n");
  printf("        1: par-2 w/o LET-99\n");
  printf("        2: par-2 with LET-99\n");
  printf("        3: par-3 w/o LET-99\n");
  printf("        4: par-3 with LET-99\n");
  printf("        5: let-99\n");
  printf("        6: ric-8\n");
  printf("        7: let-99;ric-8\n");
  printf("        8: gpr-1/2\n");
  printf("        9: PosteriorCortexPullingOnly\n");
}
