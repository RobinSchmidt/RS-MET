#ifndef KeyValidator_h
#define KeyValidator_h

/**

This class generates key-strings for the commercial product-line.

*/

class KeyValidator
{

public:

 /**< This is an enumeration of the product for which keys can be generated. */
 enum productIndices
 {
  MAGIC_CARPET = 1,
  FDN_REVERB,
  HIGH_ORDER_EQUALIZER
 };

 //---------------------------------------------------------------------------
 // construction/destruction:

 KeyValidator();
 ~KeyValidator();

 //---------------------------------------------------------------------------
 // setup:

 void setProductIndex(unsigned long newProductIndex);
 /**< Sets the product index (see enum above) for which the key has to be 
      generated. */

 void setSerialNumber(unsigned long newSerialNumber);
 /**< Sets the serial number for which the key has to be generated. */



 //---------------------------------------------------------------------------
 // processing:

 void getKeyString(char* keyToWrite, int keyLength);
 /**< Generates the key according to the settings of productIndex and 
      serialNumber*/

 //===========================================================================

protected:

 unsigned long productIndex;
 unsigned long serialNumber;
};

#endif // KeyValidator_h
