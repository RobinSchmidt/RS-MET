#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  rsseq.c
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

/* Compile: gcc -lm -o rsseq rsseq.c */

int main( int argc, char *argv[] )
{
  if(argc < 2)
    {
      printf("usage: %s n\n", argv[0]);
      printf("  Generates the Rudin-Shapiro sequence\n");
      printf("  n = number of terms\n");
      exit(-1);
    }

  unsigned int i, j, k, n = (unsigned int)strtoul(argv[1],NULL,10);

  for(i=0; i<n; ++i)
    {
      for(j=i, k=0; j>0; j>>=1) if((j & 3) == 3) ++k;
      putchar(k & 1 ? '0' : '1');
    }

  putchar('\n');
  return(0);
}
