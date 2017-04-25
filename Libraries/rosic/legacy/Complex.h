#ifndef Complex_h
#define Complex_h

//#include "VstTools.h"

#include "Definitions.h"
#include <math.h>
#include <stdio.h>  // for printf

/**

This is a class for complex numbers.

*/

class Complex  
{
public:
 
	Complex();
 /**< Constructor. Initializes real and imaginary part with zero. */

 Complex(double Re);
 /**< Constructor. Initializes real part with the argument "Re" and imaginary 
      part with zero. */

 Complex(double Re, double Im);
 /**< Constructor. Initializes real and imaginary parts with the parameters. */

	virtual ~Complex(); ///< Destructor.

 //public member variables:
 double re;  ///< Real part. 
 double im;  ///< Imaginary part.

 //overloaded operators:
 /** Defines the negative of a complex number. */
 Complex operator-()
 {
  Complex result;
  result.re = -re;
  result.im = -im;
  return result;
 }

 /** Adds two complex numbers. */
 Complex operator+(const Complex& operand2) const  
 {
  Complex result;
  result.re = re + operand2.re;
  result.im = im + operand2.im;
  return result;
 }

 /** Subtracts two complex numbers. */
 Complex operator-(const Complex& operand2) const  
 {
  Complex result;
  result.re = re - operand2.re;
  result.im = im - operand2.im;
  return result;
 }

 /** Multiplies two complex numbers. */
 Complex operator*(const Complex& operand2) const  
 {
  Complex result;
  result.re = re*operand2.re - im*operand2.im;
  result.im = re*operand2.im + im*operand2.re;
  return result;
 }

 /** Divide two complex numbers. */
 Complex operator/(const Complex& operand2) const  
 {
  Complex result;
  result.re = (re*operand2.re + im*operand2.im) / 
              (operand2.re*operand2.re + operand2.im*operand2.im);
  result.im = (-re*operand2.im + im*operand2.re) /
              (operand2.re*operand2.re + operand2.im*operand2.im);
  return result;
 }

 /** Multiplies a complex number with a real number. */
 Complex operator*(const double& operand2) const  
 {
  Complex result;
  result.re = re*operand2;
  result.im = im*operand2;
  return result;
 }

 /** Divides a complex number by a real number. */
 Complex operator/(const double& operand2) const  
 {
  Complex result;
  result.re = re/operand2;
  result.im = im/operand2;
  return result;
 }

 //public member functions:

 double   getMagnitude();
 /**< Returns the magnitude of "this" complex number. */

 double   getPhase();
 /**< Returns the phase of "this" complex number. */

 Complex getConj();
 /**< Returns the complex conjugate of "this" complex number. */

 void    print();  
 /**< Prints "this" number to standard out. */

 //mathematical functions for complex numbers:
 static Complex exp(Complex z);   
 /**< Calculates the complex exponential of a complex number. */

 static Complex sqrt(Complex z);  
 /**< Calculates the (primitive) square root of a complex number.
      The second square root is obtained by using the negative value. */

 static Complex sin(Complex z);   //complex sine
 static Complex cos(Complex z);   //complex cosine

protected:


};

#endif // Complex_h
