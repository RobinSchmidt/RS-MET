#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  debruijn.c
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

// Compile: gcc -lm -o debruijn debruijn.c

int n;
int *b;
int idbs, ndbs;
char *dbs;

// k = length of necklace
// l = length of longest prefix that is a lyndon word

void neckbin(int k, int l)
{
  if(k > n){
    if((n%l)==0){
      for(k=0; k<l; ++k) dbs[idbs+k] = b[k+1]==0 ? '0' : '1';
      idbs += l;}}
  else{
    b[k]=b[k-l];
    if(b[k]==1){
      neckbin(k+1,l);
      b[k]=0;
      neckbin(k+1,k);}
    else
      neckbin(k+1,l);}
}

/*******************************************************************/

int main( int argc, char *argv[] )
{
  if(argc < 2)
    {
      printf("usage: %s n\n", argv[0]);
      printf("  Generates the largest de Bruijn sequence of order n.\n");
      exit(-1);
    }

  n = atoi(argv[1]);
  ndbs = 1 << n;
  idbs = 0;

  b = (int *)malloc((n+2)*sizeof(int));
  b[0] = 1;
  dbs = (char *)malloc((ndbs+1)*sizeof(char));

  neckbin(1,1);
  dbs[ndbs]='\0';
  printf("%s\n", dbs);

  free(b);
  free(dbs);
  return(0);
}
