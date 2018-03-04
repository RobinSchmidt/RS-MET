#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  partam.c
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

// Compile: gcc -lm -o partam partam.c

int *parts;
int *aparts; // allowed parts
int nap; // number of allowed parts
int mp;

int allowed(int p)
{
  int i;
  for(i=0; i<nap; ++i)
    if(p == aparts[i]) return(1);
  return(0);
}

void partition(int n, int p, int m)
{
  if( n==0 ){
    if(m == mp && allowed(p)){
      for(; n<m; ++n) printf("%d ",parts[n]);
      printf("%d\n", p);}
    return;}

  if( n<0 ) return;

  if(m < mp && allowed(p)){ parts[m]=p; partition(n-p,p,m+1); }
  partition(n-1,p+1,m);
}

/*******************************************************************/

int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
      printf("usage: %s n m p1 p2 ...n\n", argv[0]);
      printf("  Generates all partitions of n with m parts form the set (p1 p2 ...)\n");
      exit( -1 );
    }

  int i, n = atoi(argv[1]);
  parts = (int *)malloc(n*sizeof(int));
  mp = atoi(argv[2])-1;
  nap = argc - 3;
  aparts = (int *)malloc(nap*sizeof(int));
  for(i=0; i<nap; ++i) aparts[i] = atoi(argv[3+i]);

  partition(n-1,1,0);
  free(parts);
  free(aparts);
  return(0);
}
