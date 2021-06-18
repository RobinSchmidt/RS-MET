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
  ok &= testStringDoubleConversions();  // fails!
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

bool rotes::testStringDoubleConversions()
{
  bool ok = true;
  ok &= testStringDoubleConversionsRandom(10000);    // fails!
  ok &= testStringDoubleConversionsSpecialValues();
  ok &= testStringDoubleConversionsDenormals();      // fails!
  ok &= testStringDoubleConversionsLarge();          // fails!
  return ok;
}

bool rotes::testStringDoubleConversionsRandom(int numIterations)
{
  bool ok = true;
  double numberOriginal, numberReconstructed;
  rsString numString;
  for(int i=0; i<numIterations; i++)
  {
    numberOriginal      = RAPT::rsRandomUniform(-1000000.0, 1000000.0);
    numString           = numberOriginal;
    numberReconstructed = numString.asDouble();
    ok &= numberReconstructed == numberOriginal;
    //rsAssert(ok);
  }
  return ok;
  // This test fails, because the double -> string conversion (done via sprintf in 
  // rsString::initFromDoubleValue) sometimes produces wrong last decimal digits. Probably some
  // numeric precision issue in sprintf. See:
  // https://stackoverflow.com/questions/62661223/sprintf-formatting-problem-for-doubles-with-high-precision
  // https://www.exploringbinary.com/incorrect-round-trip-conversions-in-visual-c-plus-plus/
  // Maybe use:
  // https://github.com/ulfjack/ryu or
  // https://github.com/jk-jeon/Grisu-Exact
  // https://github.com/jk-jeon/Grisu-Exact11
  // https://github.com/jk-jeon/fp
  // or drag in the code from RSLib for rsBigFloat - i have implemented a parsing algo there, too 
  // and it may be useful anyway
}

bool rotes::testStringDoubleConversionsSpecialValues()
{
  bool ok = true;
  double numberOriginal      = INF;
  rsString numString           = numberOriginal;
  double numberReconstructed = numString.asDouble();
  ok &= numberReconstructed == numberOriginal;

  numberOriginal      = -INF;
  numString           = numberOriginal;
  numberReconstructed = numString.asDouble();
  ok &= numberReconstructed == numberOriginal;

  numberOriginal      = NAN;
  numString           = numberOriginal;
  numberReconstructed = numString.asDouble();
  ok &= RAPT::rsIsNaN(numberReconstructed) && RAPT::rsIsNaN(numberOriginal);
  // Direct equality check on NaNs returns always false

  /*
  // obsolete:
  if( numberReconstructed != numberOriginal )
  {
    //if( !(_isnan(numberReconstructed) && _isnan(numberOriginal)) )
    if( !(RAPT::rsIsNaN(numberReconstructed) && RAPT::rsIsNaN(numberOriginal)) )
      DEBUG_BREAK; // direct equality check on NaNs seems to always return false
  }
  */
  return ok;
}

bool rotes::testStringDoubleConversionsDenormals()
{
  bool ok = true;
  ok &= testStringDoubleConversionsGeometricProgression(1.0,  0.5/SQRT2);
  ok &= testStringDoubleConversionsGeometricProgression(1.0, -0.5/SQRT2);
  return ok;
}


bool rotes::testStringDoubleConversionsLarge()
{
  bool ok = true;
  ok &= testStringDoubleConversionsGeometricProgression(1.0,  SQRT2);
  ok &= testStringDoubleConversionsGeometricProgression(1.0, -SQRT2);
  return ok;
}

bool rotes::testStringDoubleConversionsGeometricProgression(double start, double factor)
{
  bool ok = true;

  double numberOriginal      = start;
  rsString numString           = numberOriginal;
  double numberReconstructed = numString.asDouble();
  int    iteration           = 0;

  double limit;
  double absFactor = fabs(factor);
  if( absFactor < 1.0 )
  {
    limit = 0.0;
    ok &= absFactor < 0.5;
    // factors between 0.5...1.0 will make the iteration stall at a finite denormal number
  }
  else if( absFactor > 1.0 )
    limit = INF;

  while( fabs(numberOriginal) != limit )
  {
    numberOriginal     *= factor;
    numString           = numberOriginal;
    numberReconstructed = numString.asDouble();
    ok &= numberReconstructed == numberOriginal;
    iteration++;
  }
  return ok;
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
