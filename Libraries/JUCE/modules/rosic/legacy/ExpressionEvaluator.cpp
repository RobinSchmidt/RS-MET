#include "ExpressionEvaluator.h"
#include "ExpressionEvaluatorFunctions.h"

//----------------------------------------------------------------------------
// construction/destruction:

ExpressionEvaluator::ExpressionEvaluator()
{
 // create the embedded list-"objects" (assume no errors):
 err = exprFuncListCreate(&pFunctionList);
 err = exprValListCreate(&pVariableList);
 err = exprValListCreate(&pConstantList);

 // initialize the constant and function list with the built-in constants
 // and function of the ExprEval library:
 err = exprFuncListInit(pFunctionList);
 err = exprValListInit(pConstantList);

 // create the expression object:
 err = exprCreate(&pExpressionObject, // pointer to a pointer an 
                                      // expression object
                  pFunctionList,      // pointer to a function list
                  pVariableList,      // pointer to a variable list
                  pConstantList,      // pointer to a constant list
                  NULL, 
                  0,                  //
                  NULL);              // pointer to void (userdata)


 //add constants to the constants list:
 exprValListAdd(pConstantList, "pi", 3.1415926535897932384626433832795);

 // add own functions to the function list:
 exprFuncListAdd(pFunctionList, // list, where the function should be added
                 gauss,         // name of the function to be called
                 "gauss",       // name to be typed in by the user
                 1,             // minimum number of arguments
                 3,             // maximum number of arguments
                 0,             // minimum number of reference arguments
                 0);            // maximum number of reference arguments

 // and so on:
 exprFuncListAdd(pFunctionList, belowOrAbove,  "belowOrAbove", 1, 3, 0, 0);
 exprFuncListAdd(pFunctionList, cheby,         "cheby",        1, 2, 0, 0);
 exprFuncListAdd(pFunctionList, foldOver,      "foldOver",     1, 3, 0, 0);
 exprFuncListAdd(pFunctionList, logistic,      "logistic",     1, 2, 0, 0);
 exprFuncListAdd(pFunctionList, minkowski,     "minkowski",    1, 4, 0, 0);
 exprFuncListAdd(pFunctionList, quant,         "quant",        1, 2, 0, 0);
 exprFuncListAdd(pFunctionList, quantToBits,   "quantToBits",  1, 2, 0, 0);
 exprFuncListAdd(pFunctionList, sign,          "sign",         1, 1, 0, 0);
 exprFuncListAdd(pFunctionList, sawWave,       "sawWave",      1, 1, 0, 0);
 exprFuncListAdd(pFunctionList, sqrWave,       "sqrWave",      1, 1, 0, 0);
 exprFuncListAdd(pFunctionList, triWave,       "triWave",      1, 1, 0, 0);


}

ExpressionEvaluator::~ExpressionEvaluator()
{
 exprFree(pExpressionObject);
 exprValListFree(pVariableList);
 exprValListFree(pConstantList);
 exprFuncListFree(pFunctionList);
}

//----------------------------------------------------------------------------
// parse/evaluate expressions:

bool ExpressionEvaluator::parseExpression(char* expression)
{
 err = exprParse(pExpressionObject, expression);

 if( err == EXPR_ERROR_NOERROR )
  return true;
 else
  return false;
}

bool ExpressionEvaluator::evaluateExpression(double* result)
{
 err = exprEval(pExpressionObject, result);

 if( err == EXPR_ERROR_NOERROR )
  return true;
 else
  return false;
}

//----------------------------------------------------------------------------
// variable assignments:

bool ExpressionEvaluator::assignVariable(char* name, double value)
{
 err = exprValListAdd(pVariableList, name, value);

 if( err == EXPR_ERROR_NOERROR )
  return true;
 else
  return false;
}

double* ExpressionEvaluator::getVariableAddress(char* name)
{
  return pVariableList->GetAddress(name);

 //if( err == EXPR_ERROR_NOERROR )
 // return true;
 //else
 // return false;
}



//----------------------------------------------------------------------------
// others:

bool ExpressionEvaluator::clearExpression()
{
 err = exprClear(pExpressionObject);

 if( err == EXPR_ERROR_NOERROR )
  return true;
 else
  return false;
}

int ExpressionEvaluator::getLastError()
{
 return err;
}