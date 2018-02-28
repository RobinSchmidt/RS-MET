#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  cfsurd.c
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

/* Compile: gcc -lm -o cfsurd cfsurd.c */

int main( int argc, char *argv[] )
{
  if(argc < 2)
    {
      printf("usage: %s n A B\n", argv[0]);
      printf("  Calculates continued fraction for: (A+sqrt(n))/B.\n");
      printf("  The periodic part is in parenthesis.\n");
      printf("  n, A, B = integers with n > 0\n");
      exit(-1);
    }

  long int n = strtol(argv[1],NULL,10);
  long int A0, A = strtol(argv[2],NULL,10);
  long int B0, B = strtol(argv[3],NULL,10);
  long int a;
  int periodic = 0;
  double x = (A+sqrt(n))/B;
  double y = (A-sqrt(n))/B;

  if(x < 0)
    {
      printf("Error: surd is negative.\n");
      exit(-1);
    }

  if( (n - A*A) % B > 0)
    { // redefine A, B, and n so they remain integers durring calculation
      A*=B;
      B*=B;
      n*=B;
    }

  a = x;
  if( x > 1 && y < 0 && y > -1)
    {
      printf("(%ld", a);
      A0=A;
      B0=B;
      periodic=1;
    }
  else
    printf("%ld", a);

  while( !periodic )
    {
      A=a*B-A;
      B=(n-A*A)/B;
      x=(A+sqrt(n))/B;
      y=(A-sqrt(n))/B;
      a=x;
      if( x > 1 && y < 0 && y > -1)
	{
          printf(" (%ld", a);
          A0=A;
          B0=B;
          periodic=1;
	}
      else
        printf(" %ld", a);
    }

  while( 1 )
    {
      A=a*B-A;
      B=(n-A*A)/B;
      if(A==A0 && B==B0) break;
      x=(A+sqrt(n))/B;
      a=x;
      printf(" %ld", a);
    }

  printf(")\n");
  return(0);
}
