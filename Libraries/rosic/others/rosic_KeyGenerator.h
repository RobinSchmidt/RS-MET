#ifndef rosic_KeyGenerator_h
#define rosic_KeyGenerator_h

#include <stdlib.h>

#include "../basics/rosic_HelperFunctions.h"
#include "../math/rosic_PrimeNumbers.h"
#include "../math/rosic_ElementaryFunctionsReal.h"
#include "rosic_ProductIndices.h"
#include "rosic_RandomNumberGenerator01.h"
#include "rosic_RandomNumberGenerator02.h"

namespace rosic
{

  /**

  This class generates key-strings for the commercial product-line.

  \todo make a KeyValidator-class which does only generate certain parts of the key

  */

  class KeyGenerator
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    KeyGenerator();
    /**< Constructor */

    ~KeyGenerator();
    /**< Destructor */

    //---------------------------------------------------------------------------------------------
    // setup:

    void setProductIndex(unsigned long newProductIndex);
    /**< Sets the product index (see enum above) for which the key has to be 
    generated. */

    void setSerialNumber(unsigned long newSerialNumber);
    /**< Sets the serial number for which the key has to be generated. */

    void setEncodedSerialNumber(char* newEncodedSerialNumber);
    /**< Sets the serial number for which the key has to be generated in encoded form. The 
    character array must be of length 8 plus the terminating zero. */

    void setLicenseeNameAndAddress(char* newLicenseeNameAndAddress);
    /**< Sets the name of the licensee for whom the key is to be generated. */

    //---------------------------------------------------------------------------------------------
    // inquiry:

    unsigned long getSerialNumber();
    /**< Returns the serial number as a clear integer number. */

    INLINE char* getEncodedSerialNumber();
    /**< Returns a pointer to a zero terminated string which represents the serial number in 
    encoded form. The encoded form will consist of 8 characters, additionaly there is the 
    teminating zero. The calling function must take care of deleting the array when it is not 
    needed anymore. */

    char* getKeyString(int numChars);
    /**< Generates the key according to the settings of productIndex and serialNumber. The calling 
    function must take care of deleting the array when it is not needed anymore. The key will 
    consist of numChars characters plus the terminating zero. */

    bool testKeyString(char* keyToTest, int numChars);
    /**< Tests, if the the key is valid (fits to the settings of productIndex and serialNumber). 
    If it is valid, it returns true, otherwise false. */

    //---------------------------------------------------------------------------------------------
    // de/encoding

    INLINE void encodeCharacterArray(char* x, unsigned char n);
    /**< Encodes a character array of length n. */

    void decodeCharacterArray(char* y, unsigned char n);
    /**< Decodes a character array of length n. */

    //=============================================================================================

  protected:

    INLINE unsigned long scrambleBits(unsigned long x);
    /**< Scrambles the bits of an unsigned long number. */

    unsigned long descrambleBits(unsigned long y);
    /**< De-scrambles the bits of an unsigned long number. */

    unsigned char fourBitMapping1(unsigned char x);
    /**< Maps a 4-bit word (expected to be contained in the 4 least significant bits) to another
    4 bit word. */

    unsigned char inverseFourBitMapping1(unsigned char y);
    /**< Inverts the mapping of fourBitMapping1. */

    INLINE unsigned char fourBitMapping2(unsigned char x);
    /**< Maps a 4-bit word (expected to be contained in the 4 least significant bits) to another
    4 bit word. */

    unsigned char inverseFourBitMapping2(unsigned char y);
    /**< Inverts the mapping of fourBitMapping1. */

    static const unsigned long xorConstant1 = 0x55555555;
    static const unsigned long xorConstant2 = 0x95468725;
    static const unsigned long initialState = 13;

    unsigned long productIndex;
    unsigned long serialNumber;

    //unsigned char licenseeName[9]; 
    char* licenseeNameAndAddress; 

    // we have an embedded RandomNumberGenerator object, which will be set up 
    // according to productIndex and serialNumber:
    RandomNumberGenerator01 randomNumberGenerator01;
    RandomNumberGenerator02 randomNumberGenerator02;

  };

  //-----------------------------------------------------------------------------------------------
  // these functions should be inlined to make sure that their code does not appear in any product
  // but the KeyGenerator to prevent crackers from extracting that part of code from a product:

  INLINE char* KeyGenerator::getEncodedSerialNumber()
  {
    char* encodedSerial = new char[9];
    unsigned long a = scrambleBits(serialNumber);
    unsigned long aTmp;
    char          currentChar;
    for(int i=0; i<8; i++)
    {
      aTmp             = a & 0x0000000F;              // mask to retain only the 4 rightmost bits
      currentChar      = (char) (aTmp+65);            // map 4-bit pattern to a character
      encodedSerial[i] = currentChar;                 // write character into string
      a                = a >> 4;                      // rightshift by 4 bits
    }
    encodeCharacterArray(encodedSerial, 8);
    encodedSerial[8] = '\0';                          // zero termination
    return encodedSerial;
  }

  INLINE void KeyGenerator::encodeCharacterArray(char* x, unsigned char n)
  {
    int i;
    for(i=0; i<n; i++)
      x[i] -= 65;

    unsigned char tmpState = (initialState*productIndex) % 16;
    x[0] = ( fourBitMapping1(tmpState) ^ fourBitMapping2(x[0]) ) ;
    for(i=1; i<n; i++)
      x[i] = ( fourBitMapping1(x[i-1]) ^ fourBitMapping2(x[i]) ) ;

    for(i=0; i<n; i++)
      x[i] += 65;
  }

  INLINE unsigned long KeyGenerator::scrambleBits(unsigned long x)
  {
    unsigned long y = x;

    for(int i=1; i<27; i++)
    {
      y = y ^ xorConstant1;                    // bitwise XOR
      y = bitReverse(y, (unsigned long) 32);   // permutation
      y = _rotl(y, 7);                         // permutation
      y = y ^ xorConstant2;                    // bitwise XOR
    }

    return y;
  }

  INLINE unsigned char KeyGenerator::fourBitMapping2(unsigned char x)
  {
    switch(x)
    {
    case  0: return 12;
    case  1: return  7;
    case  2: return  0;
    case  3: return  8;
    case  4: return 15;
    case  5: return 10;
    case  6: return  3;
    case  7: return  1;
    case  8: return 14;
    case  9: return  5;
    case 10: return  2;
    case 11: return 13;
    case 12: return  4;
    case 13: return  9;
    case 14: return  6;
    case 15: return 11;
    }
    //DEBUG_BREAK; // x should be in the range 0...15
    return 16;
  }

} // end namespace rosic

#endif // rosic_KeyGenerator_h
