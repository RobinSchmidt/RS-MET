#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  neckm.c
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

// Compile: gcc -lm -o neckm neckm.c

int n, n1;
int *b;

// k = length of necklace
// l = length of longest prefix that is a lyndon word
// m = number of parts (ones)

void neckbin(int k, int l, int m)
{
  if(k > n){
    if((n%l)==0 && m==n1){
      for(k=1; k<n+1; ++k) printf("%d",b[k]);
      printf("\n");}}
  else{
    b[k]=b[k-l];
    if(b[k]==1){
      neckbin(k+1,l,m+1);
      b[k]=0;
      neckbin(k+1,k,m);}
    else
      neckbin(k+1,l,m);}
}

/*******************************************************************/

int main( int argc, char *argv[] )
{
  if(argc < 2)
    {
      printf("usage: %s n m\n", argv[0]);
      printf("  Generates all binary necklaces of length n with m ones.\n");
      exit(-1);
    }

  n = atoi(argv[1]);
  n1 = atoi(argv[2]);
  b = (int *)malloc((n+2)*sizeof(int));
  b[0] = 1;

  neckbin(1,1,0);
  free(b);
  return(0);
}
