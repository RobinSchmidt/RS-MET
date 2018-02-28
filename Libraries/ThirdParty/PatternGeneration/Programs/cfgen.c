#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  cfgen.c
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

/* Compile: gcc -lm -o cfgen cfgen.c */
/* Debug compile: gcc -g -D DEBUG -lm -o cfgen cfgen.c */

typedef struct {
  char name;    // variable name (single character)
  int np;       // number of productions
  char **prod;  // list of productions
} Tvar;

typedef struct {
  int nv;     // number of variables
  Tvar *var;  // list of variables and productions
} Tgrammar;

Tgrammar cfg; // the grammar used to generate strings
char *cfword;
int nstep;
char *varlist;

int getnprod(char *line)
{
  int np = 0;
  int n = strlen(line);
  int i=0;

  while(i<n)
    {
      while(isspace(line[i])){++i;}
      if(isgraph(line[i])) ++np;
      while(isgraph(line[i])){++i;}
    }
  return(np);
}

char *prodstring(char *prod)
{
  char *pchar;

  if(strcmp(prod,"epsilon")==0)
    {
      pchar = (char *)malloc(sizeof(char));
      pchar[0] = '\0';
    }
  else
    pchar = strdup(prod);

  return(pchar);
}

int ReadGrammar(char *filename)  // read the grammar file
{
  FILE *fp;
  size_t nline=128;
  char *line = (char *)malloc(nline*sizeof(char));
  int c,nv,np;
  int i,j;
 
  fp = fopen(filename, "r");
  while((c=fgetc(fp))==';'){
    while(fgetc(fp)!='\n'){}
  }
  ungetc(c,fp);
  fscanf(fp,"%d",&nv);
  while(fgetc(fp)!='\n'){}
  cfg.nv = nv;
  cfg.var = (Tvar *)malloc(nv*sizeof(Tvar));

  for(i=0; i<nv; ++i)
    {
      fscanf(fp,"%c",&(cfg.var[i].name));
      getline( &line, &nline, fp);
      np = getnprod(line);
      cfg.var[i].np = np;
      cfg.var[i].prod = (char **)malloc(np*sizeof(char *));
      cfg.var[i].prod[0] = prodstring(strtok(line," \t\n"));
      for(j=1; j<np; ++j)
        cfg.var[i].prod[j] = prodstring(strtok(NULL," \t\n"));
    }

  free(line);
  return(1);
}

void step(int istep)
{
  int i,ivar=-1;
  int nword = strlen(cfword);
  char *pchar;
  char *cfwordtail;

  //  printf("%d %s\n", istep, cfword);

  if(istep == nstep)
    {
      for(i=0; i<nword; ++i)
	if(strchr(varlist,cfword[i])) return;
      printf("%s\n", cfword);
      return;
    }

  if(istep < nstep)
    {
      for(i=0; i<nword; ++i)
	if(pchar=strchr(varlist,cfword[i]))
	  {
	    ivar = pchar - varlist;
            pchar = cfword + i;
            break;
	  }
      if(ivar<0)
        printf("%s\n", cfword);
      else
	{
	  cfwordtail = strdup(pchar);
	  for(i=0; i<cfg.var[ivar].np; ++i)
	    { // production needs to be inserted not copied
	      strcpy(pchar,cfg.var[ivar].prod[i]);
	      strcat(cfword,cfwordtail+1);
	      step(istep+1);
	    }
	  strcpy(pchar,cfwordtail);
	  free(cfwordtail);
	}
    }
}

int main( int argc, char *argv[] )
{
  if(argc < 4)
    {
      printf("usage: %s file.cfg n s\n", argv[0]);
      printf("  Generates words of a context free grammar.\n");
      printf("  file.cfg  = grammar file\n");
      printf("  n = number of derivation steps\n");
      printf("  s = start varialble\n");
      exit(-1);
    }

  if( ReadGrammar(argv[1]) == 0 )
    {
      perror( "Error reading grammar file.\n" );
      return(-1);
    }

  int i;

#ifdef DEBUG
  int j;
  for(i=0; i<cfg.nv; ++i)
    {
      printf("%c -> ", cfg.var[i].name);
      for(j=0; j<cfg.var[i].np; ++j) printf(" %s",cfg.var[i].prod[j]);
      printf("\n");
    }
#endif

  nstep = atoi(argv[2]);
  cfword = (char *)malloc(5*nstep*sizeof(char));
  cfword[0]=argv[3][0];
  cfword[1]='\0';
  varlist = (char *)malloc((cfg.nv+1)*sizeof(char));
  for(i=0; i<cfg.nv; ++i) varlist[i]=cfg.var[i].name;
  varlist[cfg.nv]='\0';

#ifdef DEBUG
  printf("%d %s %s\n", nstep, cfword, varlist);
#endif

  step(0);

  return(0);
}
