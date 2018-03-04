#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  necka.c
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

// Compile: gcc -lm -o necka necka.c

int n;
int *b;
int *aparts;
int nap;

int allowed(int p)
{
  int i;
  for(i=0; i<nap; ++i)
    if(p == aparts[i]) return(1);
  return(0);
}

// k = length of necklace
// l = length of longest prefix that is a lyndon word
// m = number of parts (ones)
// p = size of the next part

void neckbin(int k, int l, int p)
{
  if(k > n){
    if((n%l)==0 && allowed(p) && p<=n){
      for(k=1; k<n+1; ++k) printf("%d",b[k]);
      printf("\n");}}
  else{
    b[k]=b[k-l];
    if(b[k]==1){
      if(allowed(p) || k==1) neckbin(k+1,l,1);
      b[k]=0;
      neckbin(k+1,k,p+1);}
    else
      neckbin(k+1,l,p+1);}
}

/*******************************************************************/

int main( int argc, char *argv[] )
{
  if(argc < 3){
    printf("usage: %s n p1 p2 ...\n", argv[0]);
    printf("  Generates binary necklaces of length n with allowed parts pi.\n");
    exit(-1);}

  n = atoi(argv[1]);
  b = (int *)malloc((n+2)*sizeof(int));
  b[0] = 1;

  nap = argc - 2;
  aparts = (int *)malloc(nap*sizeof(int));
  int i;
  for(i=0; i<nap; ++i) aparts[i] = atoi(argv[2+i]);

  neckbin(1,1,1);
  free(b);
  free(aparts);
  return(0);
}
