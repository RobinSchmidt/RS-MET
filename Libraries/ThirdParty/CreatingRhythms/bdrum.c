#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* Reads a binary rhythm definition file and sends the equivalent abc notation to standard output. */

/*
 *                            COPYRIGHT
 *
 *  bdrum.c
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

// Compile: gcc -lm -o bdrum bdrum.c

/***********************************************************************************/

int main(int argc, char *argv[])
{
  if( argc < 3 )
    {
      printf("\nUsage: %s filename repeats\n", argv[0]);
      printf("    filename - name of rhythm definition file\n");
      printf("    repeats - number of times to repeat the rhythm\n");
      return(-1);
    }

  FILE *fp = fopen(argv[1], "r");
  int i, j, k, speed, nbeat, ndrum;
  fscanf(fp,"%d %d %d", &speed, &nbeat, &ndrum);
  int *drum = (int *)malloc(ndrum*sizeof(int));
  char **pat = (char **)malloc(ndrum*sizeof(char *));
  for(i=0; i<ndrum; ++i){
    pat[i] = (char *)malloc(nbeat*sizeof(char));
    fscanf(fp, "%d %s", drum+i, pat[i]);}
  fclose(fp);
  int repeats = atoi(argv[2]);

  printf("X: 1\n");
  printf("T: %s\n", argv[1]);
  printf("M: %d/4\n", nbeat);
  printf("K: C\n");
  printf("Q: %d\n", speed);

  //  for(i=0; i<ndrum; ++i)
  //    printf("%d %s\n", drum[i], pat[i]);

  char *beats = (char *)malloc(2*nbeat*sizeof(char));
  char *pbeats;
  char note;

  for(i=0; i<ndrum; ++i){
    beats[0]='\0';
    pbeats = beats;
    printf("V:%d clef=perc\n", i+1);
    printf("L: 1/4\n");
    printf("%%%%MIDI channel 10\n");
    printf("%%%%MIDI drummap %c %d\n", 'A'+i, drum[i]);
    for(j=0; j<nbeat;){
      k=1;
      note = pat[i][j] == '0' ? 'z' : 'A'+i;
      while(pat[i][++j]!='1' && j<nbeat) ++k;
      if(k==1)
        pbeats += sprintf(pbeats,"%c", note);
      else
        pbeats += sprintf(pbeats,"%c%d", note, k);}
    for(j=0; j<repeats; ++j) printf("| %s ", beats);
    printf("|\n");}
  return(0);
}
