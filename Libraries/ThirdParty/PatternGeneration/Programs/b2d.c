#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  b2d.c
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

/* Compile: gcc -lm -o b2d b2d.c      */
/* Converts binary number to decimal. */
/* usage: b2d n                       */

int main( int argc, char *argv[] )
{
  if(argc < 2) exit(-1);

  int i, l = strlen(argv[1]);
  unsigned long n = 0;

  for(i=0; i<l; ++i) if(argv[1][i]=='1') n += 1 << (l - 1 - i);
  printf("%lu\n", n);
  return(0);
}
