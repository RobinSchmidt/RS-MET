//#include "rosic_StringTests.h"
using namespace rotes;

//#include "rosic/rosic.h"
using namespace rosic;

bool rotes::testRosicString()
{
  bool ok = true;
  ok &= testCharacterComparisons();
  ok &= testStringBufferCopying();
  ok &= testStringIntConversions();
  //testStringDoubleConversions();  // fails - why?
  return ok;
}

bool rotes::testCharacterComparisons()
{
  bool ok = true;

  char n = '0';
  char A = 'A';
  char B = 'B';
  char a = 'a';
  char b = 'b';

  ok &= rsString::compareCharacters(n, A) == -1;
  ok &= rsString::compareCharacters(n, a) == -1;
  ok &= rsString::compareCharacters(A, A) ==  0;
  ok &= rsString::compareCharacters(A, a) == -1;
  ok &= rsString::compareCharacters(a, A) == +1;
  ok &= rsString::compareCharacters(a, a) ==  0;
  ok &= rsString::compareCharacters(A, B) == -1;
  ok &= rsString::compareCharacters(A, b) == -1;
  ok &= rsString::compareCharacters(a, B) == -1;
  ok &= rsString::compareCharacters(a, b) == -1;
  ok &= rsString::compareCharacters(B, A) == +1;
  ok &= rsString::compareCharacters(B, a) == +1;
  ok &= rsString::compareCharacters(b, A) == +1;
  ok &= rsString::compareCharacters(b, a) == +1;

  return ok;
}

bool rotes::testStringBufferCopying()
{
  bool ok = true;

  static const int charBufferLength = 20;
  char charBufferOriginal[     charBufferLength];
  char charBufferReconstructed[charBufferLength];
  RAPT::rsArrayTools::fillWithValue(charBufferOriginal,      charBufferLength, 'X');
  RAPT::rsArrayTools::fillWithValue(charBufferReconstructed, charBufferLength, 'X');
  strcpy(charBufferOriginal, "0123456789");

  rsString string;
  string.readFromBuffer(charBufferOriginal);
  string.writeIntoBuffer(charBufferReconstructed, 11);

  bool equal = strcmp(charBufferOriginal, charBufferReconstructed) == 0;
  ok &= equal == true;

  return ok;
}

bool rotes::testStringIntConversions(int numIterations)
{
  bool ok = true;
  int    numberOriginal, numberReconstructed;
  rsString numString;
  for(int i=0; i<numIterations; i++)
  {
    numberOriginal      = (int) RAPT::rsRandomUniform(INT_MIN, INT_MAX);
    numString           = numberOriginal;
    numberReconstructed = numString.asInt();
    ok &= numberReconstructed == numberOriginal;
  }
  return ok;
}

void rotes::testStringDoubleConversions()
{
  testStringDoubleConversionsRandom(10000);
  testStringDoubleConversionsSpecialValues();
  testStringDoubleConversionsDenormals();
  testStringDoubleConversionsLarge();
}

void rotes::testStringDoubleConversionsRandom(int numIterations)
{
  double numberOriginal, numberReconstructed;
  rsString numString;
  for(int i=0; i<numIterations; i++)
  {
    numberOriginal      = RAPT::rsRandomUniform(-1000000.0, 1000000.0);
    numString           = numberOriginal;
    numberReconstructed = numString.asDouble();
    rassert( numberReconstructed == numberOriginal );
  }
}

void rotes::testStringDoubleConversionsSpecialValues()
{
  double numberOriginal      = INF;
  rsString numString           = numberOriginal;
  double numberReconstructed = numString.asDouble();
  rassert( numberReconstructed == numberOriginal );


  numberOriginal      = -INF;
  numString           = numberOriginal;
  numberReconstructed = numString.asDouble();
  rassert( numberReconstructed == numberOriginal );

  numberOriginal      = NAN;
  numString           = numberOriginal;
  numberReconstructed = numString.asDouble();
  if( numberReconstructed != numberOriginal )
  {
    //if( !(_isnan(numberReconstructed) && _isnan(numberOriginal)) )
    if( !(RAPT::rsIsNaN(numberReconstructed) && RAPT::rsIsNaN(numberOriginal)) )
      DEBUG_BREAK; // direct equality check on NaNs seems to always return false
  }
}

void rotes::testStringDoubleConversionsDenormals()
{
  testStringDoubleConversionsGeometricProgression(1.0,  0.5/SQRT2);
  testStringDoubleConversionsGeometricProgression(1.0, -0.5/SQRT2);
}


void rotes::testStringDoubleConversionsLarge()
{
  testStringDoubleConversionsGeometricProgression(1.0,  SQRT2);
  testStringDoubleConversionsGeometricProgression(1.0, -SQRT2);
}

void rotes::testStringDoubleConversionsGeometricProgression(double start, double factor)
{
  double numberOriginal      = start;
  rsString numString           = numberOriginal;
  double numberReconstructed = numString.asDouble();
  int    iteration           = 0;

  double limit;
  double absFactor = fabs(factor);
  if( absFactor < 1.0 )
  {
    limit = 0.0;
    rassert(absFactor < 0.5);
    // factors between 0.5...1.0 will make the iteration stall at a finite denormal number
  }
  else if( absFactor > 1.0 )
    limit = INF;

  while( fabs(numberOriginal) != limit )
  {
    numberOriginal     *= factor;
    numString           = numberOriginal;
    numberReconstructed = numString.asDouble();
    rassert( numberReconstructed == numberOriginal );
    iteration++;
  }
}

rsString rotes::createStringWithAllCharacters()
{
  char cString[256];
  for(int i=0; i<256; i++)
    cString[256-i-1] = i;  // backwards for compliance with C strings
  return rsString(cString);
}

rsString rotes::createStringWithAllPrintableCharacters()
{
  static const int firstPrintableIndex = 32;   // whitespace ' '
  static const int lastPrintableIndex  = 127;  // tilde '~'
  static const int numPrintables       = lastPrintableIndex-firstPrintableIndex;

  char cString[numPrintables+1];
  for(int i=0; i<numPrintables; i++)
    cString[i] = i+firstPrintableIndex;
  cString[numPrintables] = '\0';
  return rsString(cString);
}
