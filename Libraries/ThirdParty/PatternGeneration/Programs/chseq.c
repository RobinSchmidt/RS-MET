#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  chseq.c
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

/* Compile: gcc -lm -o chseq chseq.c */

int main( int argc, char *argv[] )
{
  if(argc < 3)
    {
      printf("usage: %s p q n\n", argv[0]);
      printf("  Generates the christoffel word for p/q\n");
      printf("  p = numerator\n");
      printf("  q = denominator\n");
      printf("  n = number of terms to generate, default=p+q\n");
      exit(-1);
    }

  unsigned long a, p = strtoul(argv[1],NULL,10);
  unsigned long b, q = strtoul(argv[2],NULL,10);
  unsigned long n = argc == 4 ? strtoul(argv[3],NULL,10) : p+q;
  unsigned long i = 0;

  do{
    printf("0");
    ++i;
    for(a=p, b=q; a!=b && i<n;)
      {
	if(a>b)
	  {
	    printf("1");
	    b+=q;
	  }
	else
	  {
	    printf("0");
	    a+=p;
	  }
        ++i;
	if(a==b && i<n)
	  {
	    printf("1");
	    ++i;
	  }
      }
  }while(i<n);

  printf("\n");
  return(0);
}
