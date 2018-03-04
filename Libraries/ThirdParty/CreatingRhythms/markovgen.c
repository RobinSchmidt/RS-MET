#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  markovgen.c
 *  Copyright (C) 2014 Exstrom Laboratories LLC
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  A copy of the GNU General Public License is available on the internet at:
 *  http://www.gnu.org/copyleft/gpl.html
 *
 *  or you can write to:
 *
 *  The Free Software Foundation, Inc.
 *  675 Mass Ave
 *  Cambridge, MA 02139, USA
 *
 *  Exstrom Laboratories LLC contact:
 *  stefan(AT)exstrom.com
 *
 *  Exstrom Laboratories LLC
 *  Longmont, CO 80503, USA
 *
 */

// Compile: gcc -lm -o markovgen markovgen.c

int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
      printf("usage: %s mfile s n\n", argv[0]);
      printf("  Generates random numbers using a Markov chain\n");
      printf("  mfile = transition matrix file name\n");
      printf("  s = starting state\n");
      printf("  n = how many random numbers to generate\n");
      exit( -1 );
    }

  int i, j, k;
  FILE *fp = fopen(argv[1], "r");
  int ns; // number of states
  fscanf(fp,"%d",&ns);
  double *M = (double *)malloc(ns*ns*sizeof(double));
  double *pM, u, x;
  int s = atoi(argv[2]);
  int n = atoi(argv[3]);

  for(i=0, pM=M; i<ns; ++i, pM+=ns)
    for(j=0; j<ns; ++j) fscanf(fp,"%lf", pM + j);
  fclose(fp);

  for(i=0, pM=M; i<ns; ++i, pM+=ns){
    for(j=0; j<ns; ++j) printf("%lf ", pM[j]);
    printf("\n");}

  srand(time(0));

  for(i=0; i<n; ++i){
    printf("%d ", s);
    pM = M + ns * s;
    u = (double)rand()/RAND_MAX;
    for(j=0, x=0.0; j<ns; ++j){
      x += pM[j];
      if(u < x){s = j; break;}}}

  printf("\n");
  return(0);
}
