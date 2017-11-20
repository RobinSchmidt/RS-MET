#include "KeyValueStringTreeTests.h"

bool testKeyValueStringTree(std::string &reportString)
{
  std::string testName = "rsKeyValueStringTree";
  bool testResult = true;

  testResult &= testKeyValueStringTreeToString(reportString);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testKeyValueStringTreeToString(std::string &reportString)
{
  std::string testName = "rsKeyValueStringTreeToString";
  bool testResult = true;

  rsKeyValueStringTree root(rsString("RootNode"));
  root.setStringValue(rsString("SomeStringValue"), rsString("Blah"));
  root.setStringValue(rsString("SomeRealValue"), 1.35);
  root.setIntegerValue(rsString("SomeIntegerValue"), 42);
  root.setIntegerValue(rsString("SomeBooleanValue"), true);

  rsKeyValueStringTree *subNode1 = new rsKeyValueStringTree(rsString("SubNode1"));
  root.hangInSubTree(subNode1);
  subNode1->setStringValue(rsString("SomeStringValue"), rsString("Blub"));
  subNode1->setIntegerValue(rsString("SomeIntegerValue"), 34);

  rsKeyValueStringTree *subNode2 = new rsKeyValueStringTree(rsString("SubNode2"));
  root.hangInSubTree(subNode2);
  subNode2->setStringValue(rsString("SomeStringValue"), rsString("Gobbledegook"));
  subNode2->setRealNumberValue(rsString("SomeRealValue"), PI);

  rsKeyValueStringTree *subSubNode21 = new rsKeyValueStringTree(rsString("SubSubNode21"));
  subNode2->hangInSubTree(subSubNode21);
  subSubNode21->setStringValue(rsString("SomeStringValue"), rsString("Gazonk"));
  subSubNode21->setRealNumberValue(rsString("SomeRealValue"), SQRT2);

  rsKeyValueStringTree *nestedSubNode2 = new rsKeyValueStringTree(rsString("SubNode2"));
  subNode2->hangInSubTree(nestedSubNode2);
  nestedSubNode2->setStringValue(rsString("IAmThe"), rsString("NestedSubNodeWithSameNameAsParent"));

  rsKeyValueStringTree *subNode3 = new rsKeyValueStringTree(rsString("SubNode3"));
  root.hangInSubTree(subNode3);

  rsKeyValueStringTree *subSubNode31 = new rsKeyValueStringTree(rsString("SubSubNode31"));
  subNode3->hangInSubTree(subSubNode31);
  subSubNode31->setStringValue(rsString("SomeStringValue"), rsString("Blooh"));
  subSubNode31->setRealNumberValue(rsString("SomeRealValue"), 1/SQRT2);

  rsString treeString = root.getAsXmlDocument();
  //treeString.printToStandardOutput();

  rsKeyValueStringTree rootReconstructed("RootNodeReconstructed");
  rootReconstructed.setFromXmlDocument(treeString);
  //rootReconstructed.getAsXmlDocument().printToStandardOutput();

  testResult &= ( root == rootReconstructed );

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}
