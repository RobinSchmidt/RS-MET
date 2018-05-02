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
  x = rsRoundToInt( 1.2); testResult &= (x==1);
  x = rsRoundToInt( 1.8); testResult &= (x==2);
  x = rsRoundToInt( 1.5); testResult &= (x==2);
  x = rsRoundToInt( 2.5); testResult &= (x==3);
  x = rsRoundToInt(-1.2); testResult &= (x==-1);
  x = rsRoundToInt(-1.8); testResult &= (x==-2);
  x = rsRoundToInt(-1.5); testResult &= (x==-1);
  x = rsRoundToInt(-2.5); testResult &= (x==-2);

  x = rsFloorInt( 1.2); testResult &= (x==1);
  x = rsFloorInt( 1.8); testResult &= (x==1);
  x = rsFloorInt( 1.5); testResult &= (x==1);
  x = rsFloorInt( 2.5); testResult &= (x==2);
  x = rsFloorInt(-1.2); testResult &= (x==-2);
  x = rsFloorInt(-1.8); testResult &= (x==-2);
  x = rsFloorInt(-1.5); testResult &= (x==-2);
  x = rsFloorInt(-2.5); testResult &= (x==-3);

  x = rsCeilInt( 1.2); testResult &= (x==2);
  x = rsCeilInt( 1.8); testResult &= (x==2);
  x = rsCeilInt( 1.5); testResult &= (x==2);
  x = rsCeilInt( 2.5); testResult &= (x==3);
  x = rsCeilInt(-1.2); testResult &= (x==-1);
  x = rsCeilInt(-1.8); testResult &= (x==-1);
  x = rsCeilInt(-1.5); testResult &= (x==-1);
  x = rsCeilInt(-2.5); testResult &= (x==-2);

  x = rsTruncateToInt( 1.2); testResult &= (x==1);
  x = rsTruncateToInt( 1.8); testResult &= (x==1);
  x = rsTruncateToInt( 1.5); testResult &= (x==1);
  x = rsTruncateToInt( 2.5); testResult &= (x==2);
  x = rsTruncateToInt(-1.2); testResult &= (x==-1);
  x = rsTruncateToInt(-1.8); testResult &= (x==-1);
  x = rsTruncateToInt(-1.5); testResult &= (x==-1);
  x = rsTruncateToInt(-2.5); testResult &= (x==-2);


  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testExponentExtraction(std::string &reportString)
{
  std::string testName = "rsDoubleIntConversions";
  bool testResult = true;

  int exponent;

  // double:
  exponent = rsExtractExponentFromDouble(0.49999); testResult &= ( exponent == -2 );
  exponent = rsExtractExponentFromDouble(0.99999); testResult &= ( exponent == -1 );
  exponent = rsExtractExponentFromDouble(1.0);     testResult &= ( exponent ==  0 );
  exponent = rsExtractExponentFromDouble(1.99999); testResult &= ( exponent ==  0 );
  exponent = rsExtractExponentFromDouble(2.0);     testResult &= ( exponent ==  1 );
  exponent = rsExtractExponentFromDouble(2.00001); testResult &= ( exponent ==  1 );
  exponent = rsExtractExponentFromDouble(3.99999); testResult &= ( exponent ==  1 );
  exponent = rsExtractExponentFromDouble(4.0);     testResult &= ( exponent ==  2 );
  exponent = rsExtractExponentFromDouble(4.00001); testResult &= ( exponent ==  2 );

  // float:
  exponent = rsExtractExponentFromFloat(0.49999f); testResult &= ( exponent == -2 );
  exponent = rsExtractExponentFromFloat(0.99999f); testResult &= ( exponent == -1 );
  exponent = rsExtractExponentFromFloat(1.0f);     testResult &= ( exponent ==  0 );
  exponent = rsExtractExponentFromFloat(1.99999f); testResult &= ( exponent ==  0 );
  exponent = rsExtractExponentFromFloat(2.0f);     testResult &= ( exponent ==  1 );
  exponent = rsExtractExponentFromFloat(2.00001f); testResult &= ( exponent ==  1 );
  exponent = rsExtractExponentFromFloat(3.99999f); testResult &= ( exponent ==  1 );
  exponent = rsExtractExponentFromFloat(4.0f);     testResult &= ( exponent ==  2 );
  exponent = rsExtractExponentFromFloat(4.00001f); testResult &= ( exponent ==  2 );

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}



