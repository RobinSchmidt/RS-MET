#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  autogen.c
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

/* Compile: gcc -lm -o autogen autogen.c */
/* Debug compile: gcc -D DEBUG -lm -o autogen autogen.c */

typedef struct {
  char name;     // state name (single character)
  int nt;        // number of transitions
  char **trans;  // list of transitions in the form: aS meaning the symbol a
} Tstate;        // causes a transition to state S

typedef struct {
  int ns;        // number of states
  Tstate *state; // list of states
} Tautomaton;

Tautomaton aut;
char *aword;
int nstep;
char startstate;
int nes;       // number of end states
char *endstate; // list of end states

int getntrans(char *line)
{
  int nt = 0;
  int n = strlen(line);
  int i=0;

  while(i<n)
    {
      while(isspace(line[i])){++i;}
      if(isgraph(line[i])) ++nt;
      while(isgraph(line[i])){++i;}
    }
  return(nt);
}

int ReadAutomaton(char *filename)  // read the automaton file
{
  FILE *fp;
  size_t nline=128;
  char *line = (char *)malloc(nline*sizeof(char));
  int c,ns,nt;
  int i,j;
 
  fp = fopen(filename, "r");
  while((c=fgetc(fp))==';'){
    while(fgetc(fp)!='\n'){}
  }
  ungetc(c,fp);
  fscanf(fp,"%d",&ns);
  while(fgetc(fp)!='\n'){}
  aut.ns = ns;
  aut.state = (Tstate *)malloc(ns*sizeof(Tstate));

  for(i=0; i<ns; ++i)
    {
      fscanf(fp, "%c", &(aut.state[i].name));
      getline( &line, &nline, fp);
      nt = getntrans(line);
      aut.state[i].nt = nt;
      aut.state[i].trans = (char **)malloc(nt*sizeof(char *));
      aut.state[i].trans[0] = strdup(strtok(line," \t\n"));
      for(j=1; j<nt; ++j)
        aut.state[i].trans[j] = strdup(strtok(NULL," \t\n"));
    }

  free(line);
  return(1);
}

void step(int istep)
{
  int i,istate;

  if(istep == nstep)
    {
      for(i=0; i<nes; ++i)
	if(aword[istep] == endstate[i]) // aword[istep] = current state
          { // we are at an end state so the word is good
	    aword[istep] = '\0';
            printf("%s\n", aword);
	    break;
	  }
      return;
    }

  if(istep < nstep)
    {
      for(i=0; i<aut.ns; ++i)
	  if(aword[istep] == aut.state[i].name)
	    { // index of the current state is i
	      istate=i;
	      break;
	    }
      for(i=0; i<aut.state[istate].nt; ++i)
	{ // execute all transitions from the current state
          strcpy(aword+istep,aut.state[istate].trans[i]);
          step(istep+1);
	}
    }
}

int main( int argc, char *argv[] )
{
  if( argc < 5 )
    {
      printf( "usage: %s file.aut n s e1 e2 ...\n", argv[0] );
      printf( "  Generates all words of a given length accepted by an automaton\n" );
      printf( "  file.aut  = automaton file\n" );
      printf( "  n = length of words\n" );
      printf( "  s = start state\n" );
      printf( "  ei = end state i\n" );
      exit( -1 );
    }

  if( ReadAutomaton(argv[1]) == 0 )
    {
      perror( "Error reading automaton file.\n" );
      return(-1);
    }

  int i;

#ifdef DEBUG
  int j;
  printf("Number of states = %d\n",aut.ns);
  for(i=0; i<aut.ns; ++i)
    {
      printf("%c -> ", aut.state[i].name);
      for(j=0; j<aut.state[i].nt; ++j) printf(" %s",aut.state[i].trans[j]);
      printf("\n");
    }
#endif

  nstep = atoi(argv[2]);
  startstate = argv[3][0];
  aword = (char *)malloc((nstep+2)*sizeof(char));
  aword[0]=startstate;
  aword[1]='\0';
  nes = argc - 4; // number of endstates
  endstate = (char *)malloc(nes*sizeof(char));
  for(i=0; i<nes; ++i) endstate[i]=argv[4+i][0];
  step(0);

  return(0);
}
