#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*
 *                            COPYRIGHT
 *
 *  lsysgen.c
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

// Compile: gcc -lm -o lsysgen lsysgen.c

#define MAXLINE 256
#define MAXSYMB 32
#define MAXPROD 32

typedef struct {
  int prodcount;
  char *prodstr[MAXPROD];
} TSYMB;

char *axiom;
char alphabet[MAXSYMB+1];
TSYMB symbol[MAXSYMB];
int nsymbol = 0;

/***********************************************************************************/

void generate(int igen, char *prodstr)
{
    int i,j;
    char *ptr;
    int n = strlen(prodstr);

    if( igen == 0 )
      printf("%s", prodstr);
    else
      for( i=0; i<n; ++i )
      {
        ptr = strchr(alphabet, prodstr[i]);
        if( ptr == NULL )
          printf("%c", prodstr[i]);
        else
        {
          j = ptr-alphabet;
          generate(igen-1, symbol[j].prodstr[0]);
        }
      }
}

/***********************************************************************************/

void addproduction(char *production)
{
  char *ptr;
  int i,j;

  if(nsymbol == 0)
  {
    alphabet[0] = production[0];
    alphabet[1] = '\0';
    symbol[0].prodstr[0] = strdup(&production[3]);
    symbol[0].prodcount = 1;
    ++nsymbol;
  }
  else
  {
    ptr = strchr(alphabet, production[0]);
    if( ptr == NULL )
    {
      if( nsymbol < MAXSYMB )
      {
        alphabet[nsymbol] = production[0];
        alphabet[nsymbol+1] = '\0';
        symbol[nsymbol].prodstr[0] = strdup(production+3);
        symbol[nsymbol].prodcount = 1;
        ++nsymbol;
      }
    }
    else
    {
      i=ptr-alphabet;
      j=symbol[i].prodcount;
      if(j < MAXPROD)
      {
        symbol[i].prodstr[j] = strdup(production+3);
        symbol[i].prodcount++;
      }
    }
  }
}

/***********************************************************************************/

int Readldf(char *filename)  // read the ldf file
{
    FILE *fp;
    char *line = (char *)malloc(MAXLINE*sizeof(char));
    size_t n,nline=MAXLINE;

    fp = fopen(filename, "r");

    do  //skip comments
    {
      n = getline( &line, &nline, fp);
    }while( line[0] == ';');

    line[n-1]='\0';
    axiom = strdup(line);

    while( !feof(fp) )
    {
      n = getline( &line, &nline, fp);
      line[n-1]='\0';
      addproduction(line);
    }

    fclose(fp);
    free(line);
    return(1);
}

/***********************************************************************************/

int main( int argc, char *argv[] )
{
  if( argc < 2 )
    {
      printf("\nUsage: %s file.ldf n\n", argv[0]);
      printf("  file.ldf = name of file with axiom and production rules.\n");
      printf("  n = number of times to iterate.\n");
      printf("Generates an L-system string given axiom and production rules\n");
      return(-1);
    }

  /**** Read ldf file ****/
  if( Readldf( argv[1] ) == 0 )
    {
      perror( "Error reading ldf file.\n" );
      return(-1);
    }

  int ngen = atoi( argv[2] );
  generate(ngen, axiom );
  printf("\n");
  return(0);
}
