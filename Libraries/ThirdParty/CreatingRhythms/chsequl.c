#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  chsequl.c
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

/* Compile: gcc -lm -o chsequl chsequl.c */

int main( int argc, char *argv[] )
{
  if(argc < 3)
    {
      printf("usage: %s t p q n\n", argv[0]);
      printf("  Generates the upper or lower Christoffel word for p/q\n");
      printf("  t = type of word\n");
      printf("    u = upper\n");
      printf("    l = lower\n");
      printf("  p = numerator\n");
      printf("  q = denominator\n");
      printf("  n = number of terms to generate, default=p+q\n");
      exit(-1);
    }

  unsigned long a, p = strtoul(argv[2],NULL,10);
  unsigned long b, q = strtoul(argv[3],NULL,10);
  unsigned long n = argc == 5 ? strtoul(argv[4],NULL,10) : p+q;
  unsigned long i = 0;

  do{
    printf(argv[1][0]=='u' ? "1" : "0");
    ++i;
    for(a=p, b=q; a!=b && i<n; ++i)
      if(a>b){
	printf("1");
	b+=q;}
      else{
	printf("0");
	a+=p;}
    if(a==b && i<n){ printf(argv[1][0]=='u' ? "0" : "1"); ++i;}
  }while(i<n);

  printf("\n");
  return(0);
}
