#include "KeyValidator.h"

//----------------------------------------------------------------------------
// construction/destruction:

KeyValidator::KeyValidator()
{
 productIndex = 0;
 serialNumber = 0;
}

KeyValidator::~KeyValidator()
{

}

//----------------------------------------------------------------------------
// setup:

void KeyValidator::setProductIndex(unsigned long newProductIndex)
{
 productIndex = newProductIndex;
}

void KeyValidator::setSerialNumber(unsigned long newSerialNumber)
{
 serialNumber = newSerialNumber;
}

void KeyValidator::getKeyString(char *keyToWrite, int keyLength)
{

}





