#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  rndint.c
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

// Compile: gcc -lm -o rndint rndint.c

int main( int argc, char *argv[] )
{
  if( argc < 5 )
    {
      printf("usage: %s m s c n\n", argv[0]);
      printf("  Generates random numbers with specified correlation\n");
      printf("  m = range of numbers, 0 to m\n");
      printf("  s = starting number, 0 to m\n");
      printf("  c = degree of correlation\n");
      printf("      0 = total correlation (all numbers = s)\n");
      printf("      m = no correlation (each number is independent)\n");
      printf("  n = how many random numbers to generate\n");
      exit( -1 );
    }

  int i, j, k;
  int m = atoi(argv[1]);
  int s = atoi(argv[2]);
  int c = atoi(argv[3]);
  int n = atoi(argv[4]);

  srand(time(0));

  //  printf("%d ", n);

  for(i=0; i<n; ++i){
    printf("%d ", s);
    if(c > 0){
      for(j=m; j>m-c; --j){
	k = RAND_MAX / j;
        if(rand() < s*k) --s;}
      k = RAND_MAX / 2;
      for(j=0; j<c; ++j)
        if(rand() < k) ++s;}}

  printf("\n");
  return(0);
}
