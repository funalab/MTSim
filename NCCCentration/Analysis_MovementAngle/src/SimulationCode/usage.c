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
  printf(" -m # : specify model (ex. -m 1 )\n");
  printf("        0: MT_angle_fixed\n");
  printf("        1: MT_angle_variable\n");
  printf("        2: Cortex polarity(anterior:MT tips fixed, posterior:MT tips movable)\n");
}
