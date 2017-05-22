#ifndef rosic_RandomNumberGenerator01_h
#define rosic_RandomNumberGenerator01_h

//#include <string.h>
//#include <math.h>

namespace rosic
{

  /**

  This is a class which generates pseudo-random integers via the linear congruential method: 
  x[n] = (a * x[n-1] + c) % m where x is the current output, and a,c,m are constants. For good 
  randomness it is recommended to use a power of 2 for m, c an odd number, and a % 4 = 1.

  */

  class RandomNumberGenerator01
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    RandomNumberGenerator01();

    /** Destructor. */
    ~RandomNumberGenerator01();

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets the 'a' in the equation x[n] = (a * x[n-1] + c) % m */
    void setFactor(unsigned long newFactor);

    /** Sets the 'x[n-1]' in the equation x[n] = (a * x[n-1] + c) % m */
    void setState(unsigned long newState);

    /** Sets the 'c' in the equation x[n] = (a * x[n-1] + c) % m */
    void setAdditiveConstant(unsigned long newConstant);

    /** Sets the 'm' in the equation x[n] = (a * x[n-1] + c) % m */
    void setModulus(unsigned long newModulus);

    /** Sets the 'x[n-1]' in the equation x[n] = (a * x[n-1] + c) % m from a string. The string 
    must contain exactly 8 characters (excluding the terminating NULL). */
    void setStateFromString(char* theString);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the maximum possible number which may come out of this pseudo random number 
    generator (equals modulus-1). */
    unsigned long getMaximum() const { return m-1; }

    //---------------------------------------------------------------------------------------------
    // random number generation:

    /** Generates a random number via x[n] = (a * x[n-1] + c) % m  and updates 
    the state variable x.*/
    unsigned long getRandomNumber();

    //=============================================================================================

  protected:

    unsigned long a, x, c, m;

  };

} // end namespace rosic

#endif // rosic_RandomNumberGenerator01_h
