#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  b2int.c
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

// Compile: gcc -lm -o b2int b2int.c

/*
  Reads binary strings from stdin and converts them
  to interval notation.
  Example: 1010010001001000 -> 2 3 4 3 4
*/

int main(int argc, char *argv[])
{
  char line[128];
  int j, k, nbit;

  while(gets(line)){
    nbit = strlen(line);
    for(j=0; j<nbit;){
      k=1;
      while(line[++j]!='1' && j<nbit) ++k;
      printf("%d ", k);}
    printf("\n");}
  return(0);
}
