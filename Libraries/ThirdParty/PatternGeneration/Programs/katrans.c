#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  katrans.c
 *  Copyright (C) 2013 Exstrom Laboratories LLC
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

/* Compile: gcc -lm -o katrans katrans.c */
/* Debug compile: gcc -g -D DEBUG -lm -o katrans katrans.c */

typedef struct {
  int sn;        // state number
  char *out;     // output word
  int *trans;    // list of transitions
} Tstate;

typedef struct {
  int ns;        // number of states
  int nt;        // number of transitions
  Tstate *state; // list of states
} Tautomaton;

Tautomaton aut;

int ReadAutomaton(char *filename)  // read the automaton file
{
  FILE *fp;
  size_t nline=128;
  char *line = (char *)malloc(nline*sizeof(char));
  int c,i,j,sn;
 
  fp = fopen(filename, "r");
  while((c=fgetc(fp))==';'){
    while(fgetc(fp)!='\n'){}
  }
  ungetc(c,fp);

  fscanf(fp,"%d %d",&(aut.ns),&(aut.nt));
  while(fgetc(fp)!='\n'){}
  aut.state = (Tstate *)malloc(aut.ns*sizeof(Tstate));

  for(i=0; i<aut.ns; ++i)
    aut.state[i].trans = (int *)malloc(aut.nt*sizeof(int));

  for(i=0; i<aut.ns; ++i)
    {
      fscanf(fp, "%d", &sn);
      if(sn<0 || sn>=aut.ns) return(0);
      aut.state[sn].sn = sn;

      for(j=0; j<aut.nt; ++j)
        fscanf(fp, "%d", &(aut.state[sn].trans[j]));

      j=getline( &line, &nline, fp)-1;
      line[j]='\0';
      j=0;
      while(line[j]==' ') ++j;
      aut.state[sn].out = strdup(line+j);
    }

  free(line);
  return(1);
}

int main( int argc, char *argv[] )
{
  if( argc < 2 )
    {
      printf( "usage: %s file.kat", argv[0] );
      printf( "  Translates an input stream using a k-automaton\n" );
      printf( "  file.kat  = k-automaton file\n" );
      exit( -1 );
    }

  if( ReadAutomaton(argv[1]) == 0 )
    {
      perror( "Error reading k-automaton file.\n" );
      return(-1);
    }

  int c;
  int is=0;

  while( (c=getchar()) != '\n' )
    {
      is = aut.state[is].trans[c-'0'];
      printf("%s",aut.state[is].out);
    }
  printf("\n");

  return(0);
}
