#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  pfold.c
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

/* Compile: gcc -lm -o pfold pfold.c */

// oddeven finds a and b such that n = 2^a * (2*b+1)
void oddeven( unsigned int n, unsigned int *a, unsigned int *b )
{
  unsigned int k,l;

  // two's complement of n = -n or ~n + 1
  l = n & -n;  // this is 2^a
  *b = (n / l - 1)/2;
  for(k=0; l>1; ++k) l>>=1;
  *a = k;
  return;
}

int main( int argc, char *argv[] )
{
  if(argc < 4)
    {
      printf("usage: %s n m f\n", argv[0]);
      printf("  Generates fold sequences.\n");
      printf("  n = number of terms, 1,3,7,15,31,63,127,...\n");
      printf("  m = number of bits\n");
      printf("  f = function number 0 -> 2^m-1\n");
      exit(-1);
    }

  unsigned int n = (unsigned int)strtoul(argv[1],NULL,10);
  unsigned int m = (unsigned int)strtoul(argv[2],NULL,10);
  unsigned int f = (unsigned int)strtoul(argv[3],NULL,10);
  unsigned int i, j, k;
  unsigned int b;

  for(i=1; i<=n; ++i)
    {
      oddeven(i, &k, &j);
      k = k % m;
      b =  f & (1 << k) ? 1 : 0;
      if((2*j+1) % 4 > 1) b = 1 - b;
      printf("%u", b);
    }

  printf("\n");
  return(0);
}
