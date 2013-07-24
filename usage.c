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
  printf(" -t # : specify simulation time (ex. -t 100 )\n");
  printf(" -s # : specify simulation step (ex. -s 100 )\n");
  printf(" -d # : specify simulation delta (ex. -d 0.01 [default:1/4096])\n");
  printf(" -m # : specify model (ex. -m 5 )\n");
  printf("        0: Pushing\n");
  printf("        1: Pulling\n");
  printf("        2: Cortical-Anchoring\n");
  printf("        3: Pulling Supplemental 1\n");
  printf("        4: Pushing Supplemental 2\n");
  printf("        5: Length-Dependent + Cortical-Anchoring\n");
}
