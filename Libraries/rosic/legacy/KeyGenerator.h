#ifndef KeyGenerator_h
#define KeyGenerator_h

#include "ProductIndices.h"
#include "RandomNumberGenerator.h"

/**

This class generates key-strings for the commercial product-line.

*/

class KeyGenerator
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

 KeyGenerator();
 ~KeyGenerator();

 //---------------------------------------------------------------------------
 // setup:

 void setProductIndex(unsigned long newProductIndex);
 /**< Sets the product index (see enum above) for which the key has to be 
      generated. */

 void setSerialNumber(unsigned long newSerialNumber);
 /**< Sets the serial number for which the key has to be generated. */



 //---------------------------------------------------------------------------
 // processing:

 void getKeyString(char* keyToWrite, int numChars);
 /**< Generates the key according to the settings of productIndex and 
      serialNumber. */

 bool testKeyString(char* keyToTest, int numChars);
 /**< Tests, if the the key is valid (fits to the settings of productIndex and 
      serialNumber). If it is valid, it returns true, otherwise false. */

 //===========================================================================

protected:

 unsigned long productIndex;
 unsigned long serialNumber;

 // we have an embedded RandomNumberGenerator object, which will be set up 
 // according to productIndex and serialNumber:
 RandomNumberGenerator randomNumberGenerator;
};

#endif // KeyGenerator_h
