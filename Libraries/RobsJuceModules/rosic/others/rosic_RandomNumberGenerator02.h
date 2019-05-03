#ifndef rosic_RandomNumberGenerator02_h
#define rosic_RandomNumberGenerator02_h

namespace rosic
{

  /**

  This is a class which generates pseudo-random integers via ...

  */

  class RandomNumberGenerator02
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    RandomNumberGenerator02();
    /**< Constructor. */

    ~RandomNumberGenerator02();
    /**< Destructor. */

    //---------------------------------------------------------------------------------------------
    // setup:

    void setState(unsigned long newState);
    /**< Sets the old value of 'x' in the assignement x = (x + (x*x | c)) % m. */

    //---------------------------------------------------------------------------------------------
    // random number generation:

    unsigned long getRandomNumber();
    /**< Generates a random number via x = (x + (x*x | c)) % m and updates the state 
    variable x.*/

    //=============================================================================================

  protected:

    unsigned long x, c, m;

  };

} // end namespace rosic

#endif // rosic_RandomNumberGenerator02_h
