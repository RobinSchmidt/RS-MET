#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  compa.c
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

// Compile: gcc -lm -o compa compa.c

int *parts;
int *aparts; // allowed parts
int nap; // number of allowed parts

int allowed(int p)
{
  int i;
  for(i=0; i<nap; ++i)
    if(p == aparts[i]) return(1);
  return(0);
}

void compose(int n, int p, int m)
{
  if( n==0 ){
    if(allowed(p)){
      for(; n<m; ++n) printf("%d ",parts[n]);
      printf("%d\n", p);}
    return;}

  if(allowed(p)){ parts[m]=p; compose(n-1,1,m+1); }
  compose(n-1,p+1,m);
}

/*******************************************************************/

int main( int argc, char *argv[] )
{
  if( argc < 3 )
    {
      printf("usage: %s n p1 p2 ...\n", argv[0]);
      printf("  Generates compositions of n with allowed parts pi.\n");
      exit( -1 );
    }

  int i, n = atoi(argv[1]);
  parts = (int *)malloc(n*sizeof(int));
  nap = argc - 2;
  aparts = (int *)malloc(nap*sizeof(int));
  for(i=0; i<nap; ++i) aparts[i] = atoi(argv[2+i]);

  compose(n-1,1,0);
  free(parts);
  free(aparts);
  return(0);
}
