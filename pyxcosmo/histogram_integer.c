/*
 example interface c and python  with ctypes
https://scipy-lectures.github.io/advanced/interfacing_with_c/interfacing_with_c.html

see histo_integer.c

*/
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*  
   Compute the histogram
 */ 
void histogram_integer(long * array1, long * array2, int size1, int size2, int size3, long  * ri){
    long MYMAX2 = size1*size2;
    long *hh = (long *) calloc(MYMAX2, sizeof(long));
    long totalh = 0;
    long *array3 = (long *) calloc(size3, sizeof(long));
    long i = 0;
    long j = 0;
    long maxri = MYMAX2+size3+1;


	
  for (i = 0; i <size3; i++)
    {array3[i]=array1[i] + size1*array2[i];}


  // Histogramme
  for (i = 0; i <size3; i++)
    {
      hh[array3[i]]++;
    }  // end loop on i


  // Compute total h
  totalh = size3;

  // http://www.idlcoyote.com/tips/histogram_tutorial.html
  // store the reverse indices a la IDL
  // not done by numpy.digitize and numpy.searchsort

  long *pointpos = (long*) calloc(MYMAX2, sizeof(long));
  // store the indexes	


  ri[0] = MYMAX2+1;
  pointpos[0] = ri[0];

  for (i=1; i<MYMAX2;i++)
    {
      ri[i] = ri[i-1] + hh[i-1];
      pointpos[i] = ri[i];
    }

  for (i = 0; i <size3; i++)
    {
      ri[pointpos[array3[i]]] = i;
      pointpos[array3[i]]++;
    }

  ri[MYMAX2] = totalh;

  // free memory
  free(hh);
  free(array3);
  free(pointpos);

  // return 1; 
  
}

