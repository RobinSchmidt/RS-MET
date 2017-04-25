#ifndef RandomNumberGenerator_h
#define RandomNumberGenerator_h

/**

This is a class which generates pseudo-random integers via the linear 
congruential method: x[n] = (a * x[n-1] + c) % m where x is the current 
output, and a,c,m are constants. for good randomness it is recommended to use 
a power of 2 for m, c an odd number, and a % 4 = 1

*/

class RandomNumberGenerator
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

 RandomNumberGenerator();
 ~RandomNumberGenerator();

 //---------------------------------------------------------------------------
 // setup:
 void setFactor(unsigned long newFactor);
 /**< Sets the 'a' in the equation x[n] = (a * x[n-1] + c) % m */

 void setState(unsigned long newState);
 /**< Sets the 'x[n-1]' in the equation x[n] = (a * x[n-1] + c) % m */

 void setAdditiveConstant(unsigned long newConstant);
 /**< Sets the 'c' in the equation x[n] = (a * x[n-1] + c) % m */

 void setModulus(unsigned long newModulus);
 /**< Sets the 'm' in the equation x[n] = (a * x[n-1] + c) % m */

 void setFactorFromString(char* theString);
 /**< Sets the 'a' in the equation x[n] = (a * x[n-1] + c) % m from a string. 
      The string must contain at least 8 characters, otherwise access 
      vioalation will occur.*/

 void setStateFromStringAndNumber(char* theString, unsigned long theNumber);
 /**< Sets the 'x[n-1]' in the equation x[n] = (a * x[n-1] + c) % m from a string. 
      The string must contain at least 8 characters, otherwise access 
      vioalation will occur. The number can be used as product-index*/

 void setAdditiveConstantFromString(char* theString);
 /**< Sets the 'c' in the equation x[n] = (a * x[n-1] + c) % m from a string. 
      The string must contain at least 8 characters, otherwise access 
      vioalation will occur.*/

 //---------------------------------------------------------------------------
 // processing:

 unsigned long getRandomNumber();
 /**< Generates a random number via x[n] = (a * x[n-1] + c) % m  and updates 
      the state variable x.*/

 //===========================================================================

protected:

 unsigned long a, x, c, m;

};

#endif // RandomNumberGenerator_h
