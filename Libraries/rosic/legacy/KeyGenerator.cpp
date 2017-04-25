#include "KeyGenerator.h"

//----------------------------------------------------------------------------
// construction/destruction:

KeyGenerator::KeyGenerator()
{
 productIndex = 0;
 serialNumber = 0;
}

KeyGenerator::~KeyGenerator()
{

}

//----------------------------------------------------------------------------
// setup:

void KeyGenerator::setProductIndex(unsigned long newProductIndex)
{
 productIndex = newProductIndex;
}

void KeyGenerator::setSerialNumber(unsigned long newSerialNumber)
{
 serialNumber = newSerialNumber;
}

void KeyGenerator::getKeyString(char* keyToWrite, int numChars)
{
 // initialize the (pseudo-) random number generator:
 unsigned long m = 32768;  // 32768=2^15
 unsigned long a = (36739 * productIndex)  % m;
 unsigned long c = (17389 * serialNumber)  % m;
 unsigned long s = serialNumber            % m;
 randomNumberGenerator.setFactor(a);
 randomNumberGenerator.setAdditiveConstant(c);
 randomNumberGenerator.setState(s);

 // fill the char-array with random characters, derived from the generated
 // random numbers:
 unsigned long currentLongInteger;
 unsigned long bitMask = 127;
 char          currentCharacter;
 for(int i=0; i<numChars; i++)
 {
  currentLongInteger = randomNumberGenerator.getRandomNumber();

  currentCharacter   = (char) (currentLongInteger & bitMask);
   // currentCharacter is in the range 0...127

  currentCharacter   = currentCharacter % 26;
   // currentCharacter is now in the range 0...25

  currentCharacter  += 65;
   // currentCharacter is now in the range 65...90, these are the ASCII-codes
   // for uppercase letters A-Z

  keyToWrite[i]      = currentCharacter;
 }
 keyToWrite[numChars] = '\0'; // zero termination
}

bool KeyGenerator::testKeyString(char *keyToTest, int numChars)
{
 // initialize the (pseudo-) random number generator:
 unsigned long m = 32768;  // 32768=2^15
 unsigned long a = (36739 * productIndex)  % m;
 unsigned long c = (17389 * serialNumber)  % m;
 unsigned long s = serialNumber            % m;
 randomNumberGenerator.setFactor(a);
 randomNumberGenerator.setAdditiveConstant(c);
 randomNumberGenerator.setState(s);

 // fill the char-array with random characters, derived from the generated
 // random numbers:
 unsigned long currentLongInteger;
 unsigned long bitMask = 127;
 char          currentCharacter;
 bool          keyIsValid = true;
 for(int i=0; i<numChars; i++)
 {
  currentLongInteger = randomNumberGenerator.getRandomNumber();

  currentCharacter   = (char) (currentLongInteger & bitMask);
   // currentCharacter is in the range 0...127

  currentCharacter   = currentCharacter % 26;
   // currentCharacter is now in the range 0...25

  currentCharacter  += 65;
   // currentCharacter is now in the range 65...90, these are the ASCII-codes
   // for uppercase letters A-Z

  if( keyToTest[i] != currentCharacter )
   keyIsValid = false;
 }

 return keyIsValid;
}







