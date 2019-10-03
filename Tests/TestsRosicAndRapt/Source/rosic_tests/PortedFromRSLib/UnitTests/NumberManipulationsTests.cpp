#include "NumberManipulationsTests.h"

bool testNumberManipulations(std::string &reportString)
{
  std::string testName = "rsNumberManipulations";
  bool testResult = true;

  testResult &= testDoubleIntConversions(reportString);
  testResult &= testExponentExtraction(reportString);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testDoubleIntConversions(std::string &reportString)
{
  std::string testName = "rsDoubleIntConversions";
  bool testResult = true;

  int  x;
  x = roundToInt( 1.2); testResult &= (x==1);
  x = roundToInt( 1.8); testResult &= (x==2);
  x = roundToInt( 1.5); testResult &= (x==2);
  x = roundToInt( 2.5); testResult &= (x==3);
  x = roundToInt(-1.2); testResult &= (x==-1);
  x = roundToInt(-1.8); testResult &= (x==-2);
  x = roundToInt(-1.5); testResult &= (x==-1);  // fails - is -2
  x = roundToInt(-2.5); testResult &= (x==-2);  // fails - is -3
  // figure out what is the desired rounding behavior - maybe change the test criterion

  x = floorInt( 1.2); testResult &= (x==1);
  x = floorInt( 1.8); testResult &= (x==1);
  x = floorInt( 1.5); testResult &= (x==1);
  x = floorInt( 2.5); testResult &= (x==2);
  x = floorInt(-1.2); testResult &= (x==-2);
  x = floorInt(-1.8); testResult &= (x==-2);
  x = floorInt(-1.5); testResult &= (x==-2);
  x = floorInt(-2.5); testResult &= (x==-3);

  x = ceilInt( 1.2); testResult &= (x==2);
  x = ceilInt( 1.8); testResult &= (x==2);
  x = ceilInt( 1.5); testResult &= (x==2);
  x = ceilInt( 2.5); testResult &= (x==3);
  x = ceilInt(-1.2); testResult &= (x==-1);
  x = ceilInt(-1.8); testResult &= (x==-1);
  x = ceilInt(-1.5); testResult &= (x==-1);
  x = ceilInt(-2.5); testResult &= (x==-2);

  x = truncateToInt( 1.2); testResult &= (x==1);
  x = truncateToInt( 1.8); testResult &= (x==1);
  x = truncateToInt( 1.5); testResult &= (x==1);
  x = truncateToInt( 2.5); testResult &= (x==2);
  x = truncateToInt(-1.2); testResult &= (x==-1);
  x = truncateToInt(-1.8); testResult &= (x==-1);
  x = truncateToInt(-1.5); testResult &= (x==-1);
  x = truncateToInt(-2.5); testResult &= (x==-2);


  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testExponentExtraction(std::string &reportString)
{
  std::string testName = "rsDoubleIntConversions";
  bool testResult = true;

  int exponent;

  // double:
  exponent = extractExponent(0.49999); testResult &= ( exponent == -2 );
  exponent = extractExponent(0.99999); testResult &= ( exponent == -1 );
  exponent = extractExponent(1.0);     testResult &= ( exponent ==  0 );
  exponent = extractExponent(1.99999); testResult &= ( exponent ==  0 );
  exponent = extractExponent(2.0);     testResult &= ( exponent ==  1 );
  exponent = extractExponent(2.00001); testResult &= ( exponent ==  1 );
  exponent = extractExponent(3.99999); testResult &= ( exponent ==  1 );
  exponent = extractExponent(4.0);     testResult &= ( exponent ==  2 );
  exponent = extractExponent(4.00001); testResult &= ( exponent ==  2 );

  // float:
  exponent = extractExponent(0.49999f); testResult &= ( exponent == -2 );
  exponent = extractExponent(0.99999f); testResult &= ( exponent == -1 );
  exponent = extractExponent(1.0f);     testResult &= ( exponent ==  0 );
  exponent = extractExponent(1.99999f); testResult &= ( exponent ==  0 );
  exponent = extractExponent(2.0f);     testResult &= ( exponent ==  1 );
  exponent = extractExponent(2.00001f); testResult &= ( exponent ==  1 );
  exponent = extractExponent(3.99999f); testResult &= ( exponent ==  1 );
  exponent = extractExponent(4.0f);     testResult &= ( exponent ==  2 );
  exponent = extractExponent(4.00001f); testResult &= ( exponent ==  2 );

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}



