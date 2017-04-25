#ifndef DoubleNumberPair_h
#define DoubleNumberPair_h

#include <intrin.h> 
#include <xmmintrin.h>

/**

This class encapsulates the __mm128d data-type which is used to represent a 
pair (vector) of double precision numbers. It is made to take advantage of 
vector operations while still using stadard syntax for arithmetic operations.
It can be used to represent and process stereo-signals as one entity.

*/

class DoubleNumberPair  
{
public:
 
	DoubleNumberPair();
 /**< Constructor. Initializes left and right value with zero. */

 DoubleNumberPair(double values);
 /**< Constructor. Initializes bot values with "values"*/

 DoubleNumberPair(__m128d newPair);
 /**< Constructor. Initializes the number pair from a vector variable*/

 DoubleNumberPair(double leftValue, double rightValue);
 /**< Constructor. Initializes left and right value seperately. */

 ~DoubleNumberPair(); ///< Destructor.

 void   setLeft(double newLeft);
 double getLeft();
 void   setRight(double newRight);
 double getRight();

 //overloaded operators:
 /** Defines the negative of a pair of numbers. */
 /*
 DoubleNumberPair operator-()
 {
  return -theNumberPair;
 }
 */

 /** Adds two pairs of numbers. */
 DoubleNumberPair operator+(const DoubleNumberPair& operand2) const  
 {
  return DoubleNumberPair(_mm_add_pd(theNumberPair, operand2.theNumberPair));
 }

 /** Subtracts two pairs of numbers. */
 DoubleNumberPair operator-(const DoubleNumberPair& operand2) const  
 {
  return DoubleNumberPair(_mm_sub_pd(theNumberPair, operand2.theNumberPair));
 }

 /** Multiplies two pairs of numbers. */
 DoubleNumberPair operator*(const DoubleNumberPair& operand2) const  
 {
  return DoubleNumberPair(_mm_mul_pd(theNumberPair, operand2.theNumberPair));
 }

 /** Divides two pairs of numbers. */
 DoubleNumberPair operator/(const DoubleNumberPair& operand2) const  
 {
  return DoubleNumberPair(_mm_div_pd(theNumberPair, operand2.theNumberPair));
 }


 /** Adds a number pair and a single number. */
 DoubleNumberPair operator+(const double& operand2) const  
 {
  return DoubleNumberPair(_mm_add_pd(theNumberPair, _mm_set1_pd(operand2)));
 }

 /** subtracts a single number from a number pair. */
 DoubleNumberPair operator-(const double& operand2) const  
 {
  return DoubleNumberPair(_mm_sub_pd(theNumberPair, _mm_set1_pd(operand2)));
 }

 /** Multiplies a number pair with a single number. */
 DoubleNumberPair operator*(const double& operand2) const  
 {
  return DoubleNumberPair(_mm_mul_pd(theNumberPair, _mm_set1_pd(operand2)));
 }

 /** Divides a complex number by a real number. */
 DoubleNumberPair operator/(const double& operand2) const  
 {
  return DoubleNumberPair(_mm_div_pd(theNumberPair, _mm_set1_pd(operand2)));
 }

 __m128d theNumberPair;

};

#endif // DoubleNumberPair_h
