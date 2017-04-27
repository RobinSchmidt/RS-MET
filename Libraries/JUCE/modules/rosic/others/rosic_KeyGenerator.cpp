#include "rosic_KeyGenerator.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

KeyGenerator::KeyGenerator()
{
  productIndex           = 0;
  serialNumber           = 0;
  licenseeNameAndAddress = NULL;
}

KeyGenerator::~KeyGenerator()
{
  // free memory for the name:
  if( licenseeNameAndAddress != NULL )
  {
    delete[] licenseeNameAndAddress;
    licenseeNameAndAddress = NULL;
  }
}

//-------------------------------------------------------------------------------------------------
// setup:

void KeyGenerator::setProductIndex(unsigned long newProductIndex)
{
  productIndex = newProductIndex;
}

void KeyGenerator::setSerialNumber(unsigned long newSerialNumber)
{
  serialNumber = newSerialNumber;
}

void KeyGenerator::setEncodedSerialNumber(char *newEncodedSerialNumber)
{
  int  i;
  char decodedCharArray[8];
  for(i=0; i<8; i++)
    decodedCharArray[i] = newEncodedSerialNumber[i];
  decodeCharacterArray(decodedCharArray, 8);

  unsigned long accu = 0;
  unsigned long tmp;
  char          currentChar;
  for(i=0; i<8; i++)
  {
    currentChar  = decodedCharArray[7-i];
    tmp          = (unsigned long) currentChar - 65;
    accu         = accu << 4;
    accu        += tmp;
  }
  serialNumber = descrambleBits(accu);
}

void KeyGenerator::setLicenseeNameAndAddress(char *newLicenseeNameAndAddress)
{
  // free old and allocate new memory for the name:
  if( licenseeNameAndAddress != NULL )
  {
    delete[] licenseeNameAndAddress;
    licenseeNameAndAddress = NULL;
  }
  if( newLicenseeNameAndAddress != NULL )
  {
    int newLength = (int) strlen(newLicenseeNameAndAddress);
    licenseeNameAndAddress  = new char[newLength+1];
    for(int c=0; c<=newLength; c++) // the <= is valid here, because we have one more cell allocated
      licenseeNameAndAddress[c] = newLicenseeNameAndAddress[c];
  }
}

//-------------------------------------------------------------------------------------------------
// inquiry:

unsigned long KeyGenerator::getSerialNumber()
{
  return serialNumber;
}

/*
char* KeyGenerator::getEncodedSerialNumber()
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
*/

char* KeyGenerator::getKeyString(int numChars)
{
  char* keyString = new char[numChars+1];

  // initialize the (pseudo-) random number generator:
  unsigned long m = 32768;  // 32768=2^15
  unsigned long a = (36739 * productIndex)  % m;
  unsigned long c = (17389 * serialNumber)  % m;
  unsigned long s = serialNumber            % m;

  randomNumberGenerator01.setFactor(a);
  randomNumberGenerator01.setAdditiveConstant(c);
  randomNumberGenerator01.setStateFromString(licenseeNameAndAddress);
  randomNumberGenerator02.setState(s);

  // fill the char-array with random characters, derived from the generated random numbers:
  unsigned long currentLongInteger;
  unsigned long bitMask = 127;
  char          currentCharacter;
  for(int i=0; i<numChars; i++)
  {
    if( PrimeNumbers::isPrime(i) )
      currentLongInteger = randomNumberGenerator01.getRandomNumber();
    else
      currentLongInteger = randomNumberGenerator02.getRandomNumber();

    if( ((i+1) % 3) == 0 )
    {
      unsigned long tmp = randomNumberGenerator01.getRandomNumber();
      randomNumberGenerator02.setState(tmp);
      tmp = randomNumberGenerator02.getRandomNumber();
      randomNumberGenerator01.setState(tmp);
    }

    currentCharacter   = (char) (currentLongInteger & bitMask);
     // currentCharacter is in the range 0...127

    currentCharacter   = currentCharacter % 26;
     // currentCharacter is now in the range 0...25

    currentCharacter  += 65;
     // currentCharacter is now in the range 65...90, these are the ASCII-codes
     // for uppercase letters A-Z

    keyString[i]      = currentCharacter;
  }

  keyString[numChars] = '\0'; // zero termination
  return keyString;
}

bool KeyGenerator::testKeyString(char *keyToTest, int numChars)
{
  char* correctKey = getKeyString(numChars);
  bool keyIsValid = true;
  for(int i=0; i<numChars; i++)
  {
    if( keyToTest[i] != correctKey[i] )
      keyIsValid = false;
  }
  delete[] correctKey;
  return keyIsValid;
}

//-------------------------------------------------------------------------------------------------
// others:

/*
void KeyGenerator::encodeCharacterArray(char *x, unsigned char n)
{
  int i;
  for(i=0; i<n; i++)
    x[i] -= 65;

  unsigned long tmpState = (initialState*productIndex) % 16;
  x[0] = ( fourBitMapping1(tmpState) ^ fourBitMapping2(x[0]) ) ;
  for(i=1; i<n; i++)
    x[i] = ( fourBitMapping1(x[i-1]) ^ fourBitMapping2(x[i]) ) ;

  for(i=0; i<n; i++)
    x[i] += 65;
}
*/

void KeyGenerator::decodeCharacterArray(char *y, unsigned char n)
{
  int i;
  for(i=0; i<n; i++)
    y[i] -= 65;
  char* yTmp = new char[n];
  for(i=0; i<n; i++)
    yTmp[i] = y[i];

  unsigned char tmpState = (initialState*productIndex) % 16;
  y[0] = inverseFourBitMapping2( yTmp[0] ^ fourBitMapping1(tmpState) );
  for(i=1; i<n; i++)
    y[i] = inverseFourBitMapping2( yTmp[i] ^ fourBitMapping1(yTmp[i-1]) );

  for(i=0; i<n; i++)
    y[i] += 65;
  delete[] yTmp;
}

/*
unsigned long KeyGenerator::scrambleBits(unsigned long x)
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
*/

unsigned long KeyGenerator::descrambleBits(unsigned long y)
{
  unsigned long x = y;

  for(int i=1; i<27; i++)
  {
    x = x ^ xorConstant2;
    x = _rotr(x, 7);
    x = bitReverse(x, (unsigned long) 32);
    x = x ^ xorConstant1;
  }

  return x;
}

unsigned char KeyGenerator::fourBitMapping1(unsigned char x)
{
  switch(x)
  {
  case  0: return  9;
  case  1: return  3;
  case  2: return  8;
  case  3: return 14;
  case  4: return 11;
  case  5: return  7;
  case  6: return  0;
  case  7: return 12;
  case  8: return 15;
  case  9: return  2;
  case 10: return 13;
  case 11: return  4;
  case 12: return  1;
  case 13: return  6;
  case 14: return  5;
  case 15: return 10;
  }
  //DEBUG_BREAK; // x should be in the range 0...15
  return 16;
}

unsigned char KeyGenerator::inverseFourBitMapping1(unsigned char y)
{
  switch(y)
  {
  case  0: return  6;
  case  1: return 12;
  case  2: return  9;
  case  3: return  1;
  case  4: return 11;
  case  5: return 14;
  case  6: return 13;
  case  7: return  5;
  case  8: return  2;
  case  9: return  0;
  case 10: return 15;
  case 11: return  4;
  case 12: return  7;
  case 13: return 10;
  case 14: return  3;
  case 15: return  8;
  }
  //DEBUG_BREAK; // y should be in the range 0...15
  return 16;
}

/*
unsigned char KeyGenerator::fourBitMapping2(unsigned char x)
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
*/

unsigned char KeyGenerator::inverseFourBitMapping2(unsigned char y)
{
  switch(y)
  {
  case  0: return  2;
  case  1: return  7;
  case  2: return 10;
  case  3: return  6;
  case  4: return 12;
  case  5: return  9;
  case  6: return 14;
  case  7: return  1;
  case  8: return  3;
  case  9: return 13;
  case 10: return  5;
  case 11: return 15;
  case 12: return  0;
  case 13: return 11;
  case 14: return  8;
  case 15: return  4;
  }
  //DEBUG_BREAK; // y should be in the range 0...15
  return 16;
}


