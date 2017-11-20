#include "TreeTests.h"

bool testTree(std::string &reportString)
{
  std::string testName = "rsTree";
  bool testResult = true;

  testResult &= testTreeInsert(reportString);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testTreeInsert(std::string &reportString)
{
  std::string testName = "rsTreeInsert";
  bool testResult = true;

  // create the nodes:
  rsTree<int> *root = new rsTree<int>(0);
  rsTree<int> *sub1 = new rsTree<int>(1);
  rsTree<int> *sub2 = new rsTree<int>(2);
  rsTree<int> *sub3 = new rsTree<int>(3);
  rsTree<int> *sub4 = new rsTree<int>(4);
  rsTree<int> *sub5 = new rsTree<int>(5);
  rsTree<int> *sub6 = new rsTree<int>(6);

  // connect the nodes to a tree like this (dots for avoiding a multiline comment on backslash):
  //         0
  //        / \        .
  //       1   2
  //      / \          .
  //     3   4
  //        / \        .
  //       5  6

  root->hangInSubTree(sub1);
  root->hangInSubTree(sub2);
  sub1->hangInSubTree(sub3);
  sub1->hangInSubTree(sub4);
  sub4->hangInSubTree(sub5);
  sub4->hangInSubTree(sub6);

  testResult &= (root->getSubTree(0) == sub1);
  testResult &= (root->getSubTree(1) == sub2);
  testResult &= (root->getSubTree(2) == NULL);
  testResult &= (root->getParentNode() == NULL);
  testResult &= (root->getSubTree(0)->getParentNode() == root);

  testResult &= (root->getTotalNumberOfNodes() == 7);
  testResult &= (sub1->getTotalNumberOfNodes() == 5);
  testResult &= (sub2->getTotalNumberOfNodes() == 1);
  testResult &= (sub3->getTotalNumberOfNodes() == 1);
  testResult &= (sub4->getTotalNumberOfNodes() == 3);
  testResult &= (sub5->getTotalNumberOfNodes() == 1);
  testResult &= (sub6->getTotalNumberOfNodes() == 1);

  testResult &= (root->isRoot() == true);
  testResult &= (root->isLeaf() == false);
  testResult &= (sub1->isRoot() == false);
  testResult &= (sub1->isLeaf() == false);
  testResult &= (sub2->isRoot() == false);
  testResult &= (sub2->isLeaf() == true);
  testResult &= (sub3->isRoot() == false);
  testResult &= (sub3->isLeaf() == true);
  testResult &= (sub4->isRoot() == false);
  testResult &= (sub4->isLeaf() == false);
  testResult &= (sub5->isRoot() == false);
  testResult &= (sub5->isLeaf() == true);
  testResult &= (sub6->isRoot() == false);
  testResult &= (sub6->isLeaf() == true);

  delete root; // recursively deletes child-nodes

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}
