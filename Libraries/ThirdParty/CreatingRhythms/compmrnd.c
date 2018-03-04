#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  compmrnd.c
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

// Compile: gcc -lm -o compmrnd compmrnd.c

int main( int argc, char *argv[] )
{
  if( argc < 3 )
    {
      printf("usage: %s n m\n", argv[0]);
      printf("  Generate random composition of n into m parts\n");
      exit( -1 );
    }

  int j, p;
  int np, n = atoi(argv[1]);
  int mp, m = atoi(argv[2]);

  srand(time(0));

  for(mp=m-1, np=n-1, j=1; mp>0; --np){
    p = mp*(RAND_MAX/np);
    if(rand() < p){
      printf("%d ", j);
      --mp;
      j=1;}
    else
      ++j;}

  printf("%d\n", j+np);
  return(0);
}
