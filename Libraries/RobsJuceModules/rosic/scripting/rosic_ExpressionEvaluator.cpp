//#include "rosic_ExpressionEvaluator.h"
//#include "rosic_ExpressionEvaluatorFunctions.h"
//using namespace rosic;

// We need to include the .cpp files of the ExprEval library here
//#define M_PI PI // because it's used in the ExprEval code but not defined - replaced there with the literal number
//#include <math.h>
#include "../_third_party/ExprEval_v3_4/datalist.cpp"
#include "../_third_party/ExprEval_v3_4/except.cpp"
#include "../_third_party/ExprEval_v3_4/expr.cpp"
#include "../_third_party/ExprEval_v3_4/func.cpp"
#include "../_third_party/ExprEval_v3_4/funclist.cpp"
#include "../_third_party/ExprEval_v3_4/node.cpp"
#include "../_third_party/ExprEval_v3_4/parser.cpp"
#include "../_third_party/ExprEval_v3_4/vallist.cpp"

//-------------------------------------------------------------------------------------------------
// construction/destruction:

ExpressionEvaluator::ExpressionEvaluator()
{
  expressionString = NULL;

  // real functions:
  functionList.AddDefaultFunctions();
  functionList.Add(new FunctionFactoryCheby() );
  functionList.Add(new FunctionFactoryGauss() );
  functionList.Add(new FunctionFactoryQuant() );
  functionList.Add(new FunctionFactoryQuantToBits() );
  functionList.Add(new FunctionFactoryPow() );
  functionList.Add(new FunctionFactoryPowBipolar() );
  functionList.Add(new FunctionFactorySaw() );
  functionList.Add(new FunctionFactorySign() );
  functionList.Add(new FunctionFactorySoftClip() );
  functionList.Add(new FunctionFactoryStep() );
  //...

  // complex functions:
  functionList.Add(new FunctionFactoryAbsC());
  functionList.Add(new FunctionFactoryAddC());
  functionList.Add(new FunctionFactoryAngleC());
  functionList.Add(new FunctionFactoryAcosC());
  functionList.Add(new FunctionFactoryAcoshC());
  functionList.Add(new FunctionFactoryAsinC());
  functionList.Add(new FunctionFactoryAsinhC());
  functionList.Add(new FunctionFactoryAtanC());
  functionList.Add(new FunctionFactoryAtanhC());
  functionList.Add(new FunctionFactoryCosC());
  functionList.Add(new FunctionFactoryCoshC());
  functionList.Add(new FunctionFactoryDivC());
  functionList.Add(new FunctionFactoryExpC());
  functionList.Add(new FunctionFactoryLogC());
  functionList.Add(new FunctionFactoryMulC());
  functionList.Add(new FunctionFactoryPowC());
  functionList.Add(new FunctionFactorySinC());
  functionList.Add(new FunctionFactorySinhC());
  functionList.Add(new FunctionFactorySqrtC());
  functionList.Add(new FunctionFactorySubC());
  functionList.Add(new FunctionFactoryTanC());
  functionList.Add(new FunctionFactoryTanhC());

  expression.SetValueList(&valueList);
  expression.SetFunctionList(&functionList);
  setExpressionString("y=x;");

  valueList.AddDefaultValues();
  valueList.Add("pi", PI,    false);
  valueList.Add("e",  EULER, false);
  //initVariableList();

  //...
}

ExpressionEvaluator::~ExpressionEvaluator()
{
  if( expressionString != NULL )
    delete[] expressionString;
}

//-------------------------------------------------------------------------------------------------
// parse/evaluate expressions:

bool ExpressionEvaluator::setExpressionString(const char* newExpressionString)
{
  bool success = false;
  mutex.lock();

  // copy the passed c-string into our member variable:
  if( expressionString != NULL )
  {
    delete[] expressionString;
    expressionString = NULL;
  }
  if( newExpressionString != NULL )
  {
    int length = (int) strlen(newExpressionString);
    expressionString = new char[length+1];
    if( expressionString != NULL )
    {
      for(int c=0; c<=length; c++)
        expressionString[c] = newExpressionString[c];

      // try to parse the expression and return whether this was succesful or not:
      success = parseExpression();
    }
    else
      success = false;
  }
  else
    success = false;

  mutex.unlock();
  return success;
}

bool ExpressionEvaluator::parseExpression()
{
  mutex.lock();
  try
  {
    expression.Parse(expressionString);
  }
  catch( ExprEval::Exception &theException )
  {
    theException.GetValue(); // dummy instruction to suppress 'unreferenced local variable' warning
    mutex.unlock();
    return false;
  }

  // no exception was caught, so the parsing was successful:
  mutex.unlock();
  return true;
}

double ExpressionEvaluator::evaluateExpression()
{
  double result;
  mutex.lock();
  try
  {
    result = expression.Evaluate();
  }
  catch(...)
  {
    // catch evaluation errors like division by zero etc. and return 0.0 as result:
    result = 0.0;
  }
  mutex.unlock();
  return result;
}

//-------------------------------------------------------------------------------------------------
// variable assignments:


void ExpressionEvaluator::initVariableList()
{
  mutex.lock();

  valueList.Clear();
  valueList.AddDefaultValues();  // adds "E" and "PI"
  //valueList.Add("PI", PI,    false);
  //valueList.Add("E",  EULER, false);
  valueList.Add("pi", PI,    false);
  valueList.Add("e",  EULER, false);

  // we need to re-parse after messing with the value list, otherwise the expression object will
  // deal with invalid pointers:
  parseExpression();

  mutex.unlock();
}

void ExpressionEvaluator::initFunctionList()
{
  mutex.lock();

  /*
  functionList.Clear();
  functionList.AddDefaultFunctions();
  functionList.Add(new FunctionFactoryCheby() );
  functionList.Add(new FunctionFactoryGauss() );
  functionList.Add(new FunctionFactoryQuant() );
  functionList.Add(new FunctionFactoryQuantToBits() );
  functionList.Add(new FunctionFactoryPow() );
  functionList.Add(new FunctionFactoryPowBipolar() );
  functionList.Add(new FunctionFactorySaw() );
  functionList.Add(new FunctionFactorySign() );
  functionList.Add(new FunctionFactorySoftClip() );
  functionList.Add(new FunctionFactoryStep() );
  */

  mutex.unlock();
}


void ExpressionEvaluator::assignVariable(const char* name, double value)
{
  mutex.lock();

  // check, if a variable with this name already exists, if yes, overwrite its value, if not add it
  // and assign the value
  double* address = valueList.GetAddress(name);
  if( address != NULL )
  {
    *address = value;
  }
  else
  {
    valueList.Add(name, value, false);
    address = valueList.GetAddress(name);
    *address = value;
  }

  mutex.unlock();
}

double* ExpressionEvaluator::getVariableAddress(const char* name)
{
  mutex.lock();
  double* result = valueList.GetAddress(name);
  mutex.unlock();
  return result;
}

void ExpressionEvaluator::resetVariables()
{
  mutex.lock();
  valueList.Reset();
  mutex.unlock();
}

//-------------------------------------------------------------------------------------------------
// inquiry:
/*
const char* ExpressionEvaluator::getExpressionString()
{
  return expressionString;
}
*/
