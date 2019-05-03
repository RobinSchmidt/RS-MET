#ifndef rosic_RandomNumberGenerator01_h
#define rosic_RandomNumberGenerator01_h

namespace rosic
{

/** This is a class which generates pseudo-random integers via the linear congruential method:
x[n] = (a * x[n-1] + c) % m where x is the current output, and a,c,m are constants. For good
randomness it is recommended to use a power of 2 for m, c an odd number, and a % 4 = 1.

ToDo: make a class that can generate random numbers with a gaussian distribution using the 
Box-Muller transformation:
http://www.design.caltech.edu/erik/Misc/Gaussian.html */

class RandomNumberGenerator01
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Construction/Destruction */

  /** Constructor. */
  RandomNumberGenerator01();

  /** Destructor. */
  ~RandomNumberGenerator01();

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the 'a' in the equation x[n] = (a * x[n-1] + c) % m. It may modify the passed number in 
  order to satisfy some conditions for good randomness. */
  void setFactor(unsigned long newFactor);

  /** Sets the 'x[n-1]' in the equation x[n] = (a * x[n-1] + c) % m. */
  void setState(unsigned long newState);

  /** Sets the 'c' in the equation x[n] = (a * x[n-1] + c) % m. It may modify the passed number in 
  order to satisfy some conditions for good randomness. */
  void setAdditiveConstant(unsigned long newConstant);

  /** Sets the 'm' in the equation x[n] = (a * x[n-1] + c) % m. */
  void setModulus(unsigned long newModulus);

  /** Sets the 'x[n-1]' in the equation x[n] = (a * x[n-1] + c) % m from a string. The string
  must contain exactly 8 characters (excluding the terminating NULL). */
  void setStateFromString(char* theString);

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the maximum possible number which may come out of this pseudo random number
  generator (equals modulus-1). */
  unsigned long getMaximum() const { return m-1; }

  //-----------------------------------------------------------------------------------------------
  /** Random Number Generation */

  /** Generates a random number via x[n] = (a * x[n-1] + c) % m  and updates
  the state variable x.*/
  unsigned long getRandomNumber();
   // maybe inline this function


protected:

  unsigned long a, x, c, m;

};

} // end namespace rosic

#endif // rosic_RandomNumberGenerator01_h
