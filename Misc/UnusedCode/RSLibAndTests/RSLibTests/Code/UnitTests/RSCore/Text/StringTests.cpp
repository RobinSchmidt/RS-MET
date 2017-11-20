#include "StringTests.h"

rsString createStringWithAllCharacters()
{
  rsString theString;
  theString.ensureAllocatedSize(256);
  for(int i=0; i<256; i++)
    theString.appendElement(256-i-1); // backwards for compliance with C strings - 0 is last element
  return theString;
}

rsString createStringWithAllPrintableCharacters()
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

bool testString(std::string &reportString)
{
  std::string testName = "rsString";
  bool testResult = true;

  testResult &= testCharacterComparisons(reportString);
  testResult &= testStringComparisons(reportString);
  testResult &= testStringBufferCopying(reportString);
  testResult &= testStringCharacterRemoval(reportString);
  testResult &= testStringIntConversions(reportString);
  testResult &= testStringDoubleConversions(reportString);
  testResult &= testSubStringFunctions(reportString);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testCharacterComparisons(std::string &reportString)
{
  std::string testName = "rsCharacterComparisons";
  bool testResult = true;

  char n = '0';
  char A = 'A';
  char B = 'B';
  char a = 'a';
  char b = 'b';
  testResult &= ( rsString::compareCharacters(n, A) == -1 );
  testResult &= ( rsString::compareCharacters(n, a) == -1 );
  testResult &= ( rsString::compareCharacters(A, A) ==  0 );
  testResult &= ( rsString::compareCharacters(A, a) == -1 );
  testResult &= ( rsString::compareCharacters(a, A) == +1 );
  testResult &= ( rsString::compareCharacters(a, a) ==  0 );
  testResult &= ( rsString::compareCharacters(A, B) == -1 );
  testResult &= ( rsString::compareCharacters(A, b) == -1 );
  testResult &= ( rsString::compareCharacters(a, B) == -1 );
  testResult &= ( rsString::compareCharacters(a, b) == -1 );
  testResult &= ( rsString::compareCharacters(B, A) == +1 );
  testResult &= ( rsString::compareCharacters(B, a) == +1 );
  testResult &= ( rsString::compareCharacters(b, A) == +1 );
  testResult &= ( rsString::compareCharacters(b, a) == +1 );

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testStringComparisons(std::string &reportString)
{
  std::string testName = "rsStringComparisons";
  bool testResult = true;

  rsString a   = rsString("a");
  rsString b   = rsString("b");
  rsString aa  = rsString("aa");
  rsString ab  = rsString("ab");
  rsString ba  = rsString("ba");
  rsString bb  = rsString("bb");

  testResult &= (   a  <  b   );
  testResult &= (   aa <  b   );
  testResult &= (   a  <  aa  );
  testResult &= (   b  <  ba  );
  testResult &= ( !(ab <  ab) );

  testResult &= (   a  <= b  );
  testResult &= (   aa <= b  );
  testResult &= (   a  <= aa );
  testResult &= (   b  <= ba );
  testResult &= (   ab <= ab );

  testResult &= (   b  >  a   );
  testResult &= (   b  >  aa  );
  testResult &= (   aa >  a   );
  testResult &= (   ba >  b   );
  testResult &= ( !(ab >  ab) );

  testResult &= (   b  >= a   );
  testResult &= (   b  >= aa  );
  testResult &= (   aa >= a   );
  testResult &= (   ba >= b   );
  testResult &= (   ab >= ab  );

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testStringBufferCopying(std::string &reportString)
{
  std::string testName = "rsStringBufferCopying";
  bool testResult = true;

  /*
  static const int charBufferLength = 20;
  char charBufferOriginal[     charBufferLength];
  char charBufferReconstructed[charBufferLength];
  rs::fillWithValue(charBufferOriginal,      charBufferLength, 'X');
  rs::fillWithValue(charBufferReconstructed, charBufferLength, 'X');
  strcpy(charBufferOriginal, "0123456789");

  String string;
  string.readFromBuffer(charBufferOriginal);
  string.writeIntoBuffer(charBufferReconstructed, 11);

  bool equal = strcmp(charBufferOriginal, charBufferReconstructed) == 0;
  testResult &= ( equal == true );
  */

  //String *testString = new String("Balh");

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testStringCharacterRemoval(std::string &reportString)
{
  std::string testName = "rsStringCharacterRemoval";
  bool testResult = true;

  rsString testString("abcdefghijklmnopqrstuvwxyz");

  testString.removeAllCharactersBut("abcijklmnopvwxyz");
  testResult &= ( testString == rsString("abcijklmnopvwxyz") );

  testString.removeCharacters("cij");
  testResult &= ( testString == rsString("abklmnopvwxyz") );

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}


bool testStringIntConversions(std::string &reportString)
{
  std::string testName = "rsStringIntConversions";
  bool testResult = true;

  int      numIterations = 10000;
  int      numberOriginal, numberReconstructed;
  rsString numString;
  for(int i=0; i<numIterations; i++)
  {
    numberOriginal      = (int) rsRandomUniform(INT_MIN, INT_MAX);
    numString           = numberOriginal;
    numberReconstructed = numString.asInt();
    testResult &= ( numberReconstructed == numberOriginal );
  }

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testStringDoubleConversions(std::string &reportString)
{
  std::string testName = "rsStringDoubleConversions";
  bool testResult = true;

  testResult &= testStringDoubleConversionsRandom(reportString);
  testResult &= testStringDoubleConversionsSpecialValues(reportString);
  testResult &= testStringDoubleConversionsDenormals(reportString);
  testResult &= testStringDoubleConversionsLarge(reportString);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testStringDoubleConversionsRandom(std::string &reportString)
{
  std::string testName = "rsStringDoubleConversionsRandom";
  bool testResult = true;

  int numIterations = 10000;
  double numberOriginal, numberReconstructed;
  rsString numString;
  for(int i=0; i<numIterations; i++)
  {
    numberOriginal      = rsRandomUniform(-1000000.0, 1000000.0);
    numString           = numberOriginal;
    numberReconstructed = numString.asDouble();
    testResult &= ( numberReconstructed == numberOriginal );
  }

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testStringDoubleConversionsSpecialValues(std::string &reportString)
{
  std::string testName = "rsStringDoubleConversionsSpecialValues";
  bool testResult = true;

  double   numberOriginal      = rsInfDouble;
  rsString numString           = numberOriginal;
  double   numberReconstructed = numString.asDouble();
  testResult &= ( numberReconstructed == numberOriginal );

  numberOriginal      = -rsInfDouble;
  numString           = numberOriginal;
  numberReconstructed = numString.asDouble();
  testResult &= ( numberReconstructed == numberOriginal );

  numberOriginal      = rsQuietNaNDouble; 
  numString           = numberOriginal;
  numberReconstructed = numString.asDouble();
  testResult &= rsIsNaN(numberReconstructed);

  numberOriginal      = rsSignalingNaNDouble; 
  numString           = numberOriginal;
  numberReconstructed = numString.asDouble();
  testResult &= rsIsNaN(numberReconstructed);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testStringDoubleConversionsDenormals(std::string &reportString)
{
  std::string testName = "rsStringDoubleConversionsDenormals";
  bool testResult = true;

  testResult &= testStringDoubleConversionsGeometricProgression(reportString, 1.0,  0.5/SQRT2);
  testResult &= testStringDoubleConversionsGeometricProgression(reportString, 1.0, -0.5/SQRT2);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}


bool testStringDoubleConversionsLarge(std::string &reportString)
{
  std::string testName = "rsStringDoubleConversionsLarge";
  bool testResult = true;

  testResult &= testStringDoubleConversionsGeometricProgression(reportString, 1.0,  SQRT2);
  testResult &= testStringDoubleConversionsGeometricProgression(reportString, 1.0, -SQRT2);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testStringDoubleConversionsGeometricProgression(std::string &reportString, double start, double factor)
{
  bool   testResult            = true;
  double   numberOriginal      = start;
  rsString numString           = numberOriginal;
  double   numberReconstructed = numString.asDouble();
  int      iteration           = 0;

  double limit;
  double absFactor = fabs(factor);
  if( absFactor < 1.0 )
  {
    limit = 0.0;
    testResult &= (absFactor < 0.5);
    // factors between 0.5...1.0 will make the iteration stall at a finite denormal number
  }
  else if( absFactor > 1.0 )
    limit = rsInfDouble;

  while( fabs(numberOriginal) != limit )
  {
    numberOriginal     *= factor;
    numString           = numberOriginal;
    numberReconstructed = numString.asDouble();
    testResult &= ( numberReconstructed == numberOriginal );
    iteration++;
  }

  return testResult;
}

bool testSubStringFunctions(std::string &reportString)
{
  std::string testName = "rsStringSubStringFunctions";
  bool testResult = true;

  rsString testString = rsString("abcdefghabcdefg");
  int index = testString.findFirstOccurrenceOf(rsString("def"));
  testResult &= ( index == 3 );

  rsRange<int> range = testString.findRangeEnclosedBy(rsString("cde"), rsString("bcd"), true);
  testResult &= ( range.getMin() == 2 && range.getMax() == 11 );

  range = testString.findRangeEnclosedBy(rsString("cde"), rsString("bcd"), false);
  testResult &= ( range.getMin() == 5 && range.getMax() == 8 );

  testString = rsString("___a_sdf_fs___sff_dfs__sf____sdf_sdf___");
  testString.removeRepeatedCharacters('_');
  testResult &= ( testString == rsString("_a_sdf_fs_sff_dfs_sf_sdf_sdf_") );

  testString = rsString("___asdadf____");
  testString.removeLeadingCharacters('_');
  testResult &= ( testString == rsString("asdadf____") );

  testString.removeTrailingCharacters('_');
  testResult &= ( testString == rsString("asdadf") );

  testString = rsString("0123456789");
  testString = testString.getSubString(3, 6);
  testResult &= ( testString == rsString("3456") );

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}


