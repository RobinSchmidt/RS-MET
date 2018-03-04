#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *                            COPYRIGHT
 *
 *  permi.c
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

// Compile: gcc -lm -o permi permi.c

int main( int argc, char *argv[] )
{
  if( argc < 2 )
    {
      printf("usage: %s a1 a2 ... an\n", argv[0]);
      printf("  Generates permutations of the integers ai>=0.\n");
      printf("  To generate all permutations they must be ordered a1<a2<...<an.\n");
      printf("  Any other order will only generate permutations that are larger\n");
      printf("  in lexicographic order.\n");
      exit( -1 );
    }

  int i, j, k, m, n = argc - 1;
  int *a = (int *)malloc((n+1)*sizeof(int));
  a[0] = -1;
  for(i=1; i<=n; ++i) a[i] = atoi(argv[i]);
  do{
    for(i=1; i<=n; ++i) printf("%d ",a[i]);
    printf("\n");
    for(i=n-1; i>0 && a[i]>=a[i+1]; --i);
    if(i == 0) break;
    for(j=n; a[i]>=a[j]; --j);
    m = a[j];
    a[j]=a[i];
    a[i]=m;
    for(j=i+1, k=n; j<k; ++j, --k){
      m = a[j];
      a[j]=a[k];
      a[k]=m;}} while(1);

  free(a);
  return(0);
}
