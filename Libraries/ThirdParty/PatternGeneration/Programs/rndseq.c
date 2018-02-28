#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  rndseq.c
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

/* Compile: gcc -lm -o rndseq rndseq.c */

#define RAND0 (RAND_MAX - 1)/2

int main( int argc, char *argv[] )
{
  if(argc < 2)
    {
      printf("usage: %s n\n", argv[0]);
      printf("  Generates a random sequence of binary digits\n");
      printf("  n = number of digits to generate\n");
      exit(-1);
    }

  unsigned int i, n = strtoul(argv[1],NULL,10);
  unsigned int seed = (unsigned int)(argc > 2 ? strtoul(argv[2],NULL,10) : time(NULL));
  srand(seed);

  unsigned int i0, i1;

  for(i=0, i0=0, i1=0; i<n; ++i)
    if(rand() <= RAND0)
    {
      putchar('0');
      ++i0;
    }
    else
    {
      putchar('1');
      ++i1;
    }

  putchar('\n');
  return(0);
}
