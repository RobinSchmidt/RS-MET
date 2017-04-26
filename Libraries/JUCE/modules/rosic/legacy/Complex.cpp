#include "Complex.h"

Complex::Complex()
{
 re = im = 0.0;
}
Complex::Complex(double Re)
{
 re = Re;
 im = 0.0;
}
Complex::Complex(double Re, double Im)
{
 re = Re;
 im = Im;
}
Complex::~Complex()
{

}

//magnitude and phase:
double Complex::getMagnitude()
{
 return ::sqrt(re*re + im*im);
}
double Complex::getPhase()
{
 if((re==0.0) && (im==0))
  return 0;
 else
  return ::atan2(im, re);
}
//complex conjugate:
Complex Complex::getConj()
{
 Complex result;
 result.re =  re;
 result.im = -im;
 return result;
}
//output to the screen:
void Complex::print()
{
 printf("%s %.5f %s %.5f %s", "z = ", re, " + ", im, "j\n");
}

//---------------------------------------------------------------------------------------
//mathematical functions:
Complex Complex::exp(Complex z)
{
 Complex result;

 double temp = ::exp(z.re);
 double real = temp * ::cos(z.im);
 double imag = temp * ::sin(z.im);

 result.re = real;
 result.im = imag;

 return result;
}

Complex Complex::sqrt(Complex z)
{
 Complex result;

 double mag  = ::sqrt(z.getMagnitude());
 double phs  = 0.5*(z.getPhase());
 double real = mag * ::cos(phs);
 double imag = mag * ::sin(phs);

 result.re = real;
 result.im = imag;

 return result;
}