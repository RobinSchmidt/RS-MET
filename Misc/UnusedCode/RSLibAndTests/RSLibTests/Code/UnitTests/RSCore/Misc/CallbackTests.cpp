#include "CallbackTests.h"

bool testCallbacks(std::string &reportString)
{
  std::string testName = "Callbacks";
  bool testResult = true;

  testResult &= testCallbackByValueSemantics(reportString);
  testResult &= testCallbackByReferenceSemantics(reportString);

  testResult &= testParameter(reportString);

  //testResult &= testCallback2(reportString);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}


// toy-classes for which member-functions will be called for the tests - n each member-function
// invocation, they'll increment a counter which will be checked for the correct value (after some
// number of invocations have been done) inside the actual unit-tests:

class CalleeClass1
{
public:
  CalleeClass1()
  {
    doubleToDoubleCallCount1 = 0;
    doubleToDoubleCallCount2 = 0;
    floatToVoidCallCount     = 0;
  }
  double doubleToDoubleFunction1(double x)
  {
    doubleToDoubleCallCount1++;
    return x;
  }
  double doubleToDoubleFunction2(double x)
  {
    doubleToDoubleCallCount2++;
    return 2*x;
  }
  void floatToVoidFunction(float x)
  {
    floatToVoidCallCount++;
  }
  int doubleToDoubleCallCount1;
  int doubleToDoubleCallCount2;
  int floatToVoidCallCount;
};
class CalleeClass2
{
public:
  CalleeClass2()
  {
    doubleToDoubleCallCount = 0;
    doubleToVoidCallCount   = 0;
  }
  void doubleToVoidFunction(double x)
  {
    doubleToVoidCallCount++;
  }
  double doubleToDoubleFunction(double x)
  {
    doubleToDoubleCallCount++;
    return x;
  }
  int doubleToVoidCallCount;
  int doubleToDoubleCallCount;
};

// the actual tests:

bool testCallbackByValueSemantics(std::string &reportString)
{
  std::string testName = "CallbackByValueSemantics";
  bool testResult = true;

  // create some callee-objects and callbacks:
  CalleeClass1 callee1, callee2;
  rsCallback1<CalleeClass1, double, double> cb1(&callee1, &CalleeClass1::doubleToDoubleFunction1);
  rsCallback1<CalleeClass1, double, double> cb2(&callee1, &CalleeClass1::doubleToDoubleFunction2);
  rsCallback1<CalleeClass1, double, double> cb3(&callee2, &CalleeClass1::doubleToDoubleFunction1);
  rsCallback1<CalleeClass1, double, double> cb4(&callee2, &CalleeClass1::doubleToDoubleFunction2);
  rsCallback1<CalleeClass1, void,   float>  cb5(&callee1, &CalleeClass1::floatToVoidFunction);

  // invoke the callbacks and see if all member-functions were called the expected number of times:
  int i;
  double doubleResult;
  int numCalls1 = 100; for(i = 1; i <= numCalls1; i++) doubleResult = cb1.call(1.0);
  int numCalls2 = 200; for(i = 1; i <= numCalls2; i++) doubleResult = cb2.call(1.0);
  int numCalls3 = 300; for(i = 1; i <= numCalls3; i++) doubleResult = cb3.call(1.0);
  int numCalls4 = 400; for(i = 1; i <= numCalls4; i++) doubleResult = cb4.call(1.0);
  int numCalls5 = 500; for(i = 1; i <= numCalls5; i++)                cb5.call(1.f);
  testResult &= ( callee1.doubleToDoubleCallCount1 == numCalls1 );
  testResult &= ( callee1.doubleToDoubleCallCount2 == numCalls2 );
  testResult &= ( callee2.doubleToDoubleCallCount1 == numCalls3 );
  testResult &= ( callee2.doubleToDoubleCallCount2 == numCalls4 );
  testResult &= ( callee1.floatToVoidCallCount     == numCalls5 );
  if(doubleResult == 1.0) {} // we just use "doubleResult" once to avoid compiler warning

  // make a copy and see if it calls the same member-function on the same instance:
  rsCallback1<CalleeClass1, double, double> cb1copy = cb1;
  for(i = 1; i <= numCalls1; i++) doubleResult = cb1copy.call(1.0);
  testResult &= ( callee1.doubleToDoubleCallCount1 == 2*numCalls1 );

  // test comparison operators:
  testResult &= ( cb1 == cb1copy );
  testResult &= ( cb1 != cb2     );
  testResult &= ( cb1 != cb3     );
  testResult &= ( cb1 != cb4     );

  // lets store them in a vector and see if that works:
  std::vector<rsCallback1<CalleeClass1, double, double> > callbacks;
  callbacks.push_back(cb1);
  callbacks.push_back(cb2);
  callbacks.push_back(cb3);
  testResult &= ( cb1 == callbacks[0] );
  testResult &= ( cb2 == callbacks[1] );
  testResult &= ( cb3 == callbacks[2] );
  for(i = 1; i <= numCalls2; i++) doubleResult = callbacks[1].call(1.0);
  testResult &= ( callee1.doubleToDoubleCallCount2 == 2*numCalls2 );

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testCallbackByReferenceSemantics(std::string &reportString)
{
  std::string testName = "CallbackByReferenceSemantics";
  bool testResult = true;

  // create some callee-objects and callbacks:
  CalleeClass1 callee1, callee2;
  rsCallbackBase1<double, double> *cb1 = new rsCallback1<CalleeClass1, double, double>
    (&callee1, &CalleeClass1::doubleToDoubleFunction1);
  rsCallbackBase1<double, double> *cb2 = new rsCallback1<CalleeClass1, double, double>
    (&callee1, &CalleeClass1::doubleToDoubleFunction2);
  rsCallbackBase1<double, double> *cb3 = new rsCallback1<CalleeClass1, double, double>
    (&callee2, &CalleeClass1::doubleToDoubleFunction1);
  rsCallbackBase1<double, double> *cb4 = new rsCallback1<CalleeClass1, double, double>
    (&callee2, &CalleeClass1::doubleToDoubleFunction2);
  rsCallbackBase1<void, float> *cb5 =    new rsCallback1<CalleeClass1, void, float>
    (&callee1, &CalleeClass1::floatToVoidFunction);

  // \todo use a createCallback convenience function (see commented out code below)

  // invoke the callbacks and see if all member-functions were called the expected number of times:
  int i;
  double doubleResult;
  int numCalls1 = 100; for(i = 1; i <= numCalls1; i++) doubleResult = cb1->call(1.0);
  int numCalls2 = 200; for(i = 1; i <= numCalls2; i++) doubleResult = cb2->call(1.0);
  int numCalls3 = 300; for(i = 1; i <= numCalls3; i++) doubleResult = cb3->call(1.0);
  int numCalls4 = 400; for(i = 1; i <= numCalls4; i++) doubleResult = cb4->call(1.0);
  int numCalls5 = 500; for(i = 1; i <= numCalls5; i++)                cb5->call(1.f);
  testResult &= ( callee1.doubleToDoubleCallCount1 == numCalls1 );
  testResult &= ( callee1.doubleToDoubleCallCount2 == numCalls2 );
  testResult &= ( callee2.doubleToDoubleCallCount1 == numCalls3 );
  testResult &= ( callee2.doubleToDoubleCallCount2 == numCalls4 );
  testResult &= ( callee1.floatToVoidCallCount     == numCalls5 );
  if(doubleResult == 1.0) {} // we just use "doubleResult" once to avoid compiler warning

  // let's create and use a heterogenic vector of pointers to generic callbacks:
  std::vector<rsCallbackBase1<double, double>* > callbacks;
  CalleeClass2 callee3; // the to-be-called object of another class
  rsCallbackBase1<double, double> *cb6 = new rsCallback1<CalleeClass2, double, double>
    (&callee3, &CalleeClass2::doubleToDoubleFunction);
  callbacks.push_back(cb1);  // callback to object of type CallbackTestClass1
  callbacks.push_back(cb6);  // callback to object of type CallbackTestClass2
  callbacks.push_back(cb2);  // callback to object of type CallbackTestClass1
  int numCalls6 = 50; for(i = 1; i <= numCalls6; i++) doubleResult = callbacks[1]->call(1.0);
  testResult &= ( callee3.doubleToDoubleCallCount == numCalls6 );
  for(i = 1; i <= numCalls2; i++)
    doubleResult = callbacks[2]->call(1.0);
  testResult &= ( callee1.doubleToDoubleCallCount2 == 2*numCalls2 );

  // cleanup and report:
  delete cb1;
  delete cb2;
  delete cb3;
  delete cb4;
  delete cb5;
  delete cb6;
  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testParameter(std::string &reportString)
{
  std::string testName = "Parameter";
  bool testResult = true;

  // create a parameter and add two callbacks to two different objects to it:
  rsNumericParameter param;
  CalleeClass2 callee1, callee2;
  param.addCallback(new rsCallback1<CalleeClass2, void, double>
    (&callee1, &CalleeClass2::doubleToVoidFunction));
  param.addCallback(new rsCallback1<CalleeClass2, void, double>
    (&callee2, &CalleeClass2::doubleToVoidFunction));

  // map range 0...1 to 2...4, check if mapping and clipping works:
  param.setMapper(new rsRealLinearMap(0.0, 2.0, 1.0, 4.0));
  param.setNormalizedValue(0.0);
  testResult &= param.getValue() == 2.0;
  param.setNormalizedValue(0.5);
  testResult &= param.getValue() == 3.0;
  param.setNormalizedValue(1.0);
  testResult &= param.getValue() == 4.0;
  param.setValue(1.5);
  testResult &= param.getValue() == 2.0;
  param.setValue(3.5);
  testResult &= param.getValue() == 3.5;
  param.setValue(4.5);
  testResult &= param.getValue() == 4.0;

  // test exponential mapping from 0.01...100:
  double tol = 1.e-13;
  param.setMapper(new rsRealLinToExpMap(0.0, 0.01, 1.0, 100.0));
  param.setNormalizedValue(0.0);
  testResult &= rsIsCloseTo(param.getValue(), 0.01,  tol);
  param.setNormalizedValue(0.25);
  testResult &= rsIsCloseTo(param.getValue(), 0.1,   tol);
  param.setNormalizedValue(0.5);
  testResult &= rsIsCloseTo(param.getValue(), 1.0,   tol);
  param.setNormalizedValue(0.75);
  testResult &= rsIsCloseTo(param.getValue(), 10.0,  tol);
  param.setNormalizedValue(1.0);
  testResult &= rsIsCloseTo(param.getValue(), 100.0, tol);
  param.setValue(0.1);
  testResult &= rsIsCloseTo(param.getValue(), 0.1,   tol);
  param.setValue(10.0);
  testResult &= rsIsCloseTo(param.getValue(), 10.0,  tol);

  // check, if the number of calls is right for both callee objects:
  testResult &= ( callee1.doubleToVoidCallCount == 13 );
  testResult &= ( callee2.doubleToVoidCallCount == 13 );

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}


/*
// this is a test for an alternative implementation Callbacks2.h which has not been dragged over
// (yet):
bool testCallback2(std::string &reportString)
{
  std::string testName = "CallbackByValueSemantics";
  bool testResult = true;

  CallbackTestClass1 callee01, callee02;
  Callback1<double, double> callback01;

  //callback01 = SpecificCallback1<CallbackTestClass1, double, double>(&callee01, &CallbackTestClass1::doubleToDoubleFunction1);
  callback01 = makeCallback1(&callee01, &CallbackTestClass1::doubleToDoubleFunction1);

  double result = callback01.call(5.0);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}
*/
