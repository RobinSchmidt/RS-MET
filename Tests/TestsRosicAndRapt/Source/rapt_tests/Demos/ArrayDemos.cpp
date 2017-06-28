#include "ArrayDemos.h"

// maybe turn this into a unit test
void convolutionDemo()
{
  // Demonstrates use of teh convolution function. It shows that you can use different types for
  // input signal, impulse response and output signal. everything will be converted automatically
  // in the convolution function.

  static const int Nx = 10;       // length of input signal x
  static const int Nh = 5;        // length of impulse response h
  static const int Ny = Nx+Nh-1;  // length of output signal

  //// allocate memory for input, impulse response and output:
  //int    x[Nx];
  //double h[Nh];
  //float  y[Ny];

  //// fill arrays:
  //rsFillWithValue(x, Nx, 3);
  //rsFillWithValue(h, Nh, 2.0);
  //rsFillWithZeros(y, Ny);    

  //// convolve:
  //rsConvolve(x, Nx, h, Nh, y);

  int dummy = 0;
}