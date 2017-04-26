#ifndef ExpressionEvaluator_h
#define ExpressionEvaluator_h

#include "ExprEval\expreval.h"
#include <math.h>

/**

This is a class for parsing and evaluating mathematical expressions.
It encapsulates the ExprEval c-library by Brian Allen Vanderburg II,
thereby making it easier to use and augmenting it by additional functions 
and constants. For more details on the original c-library refer to the 
subdirectory "ExprEval" where the source code together with the original
documentation resides.

Much thanks to Brian Allen Vanderburg II for developing the library and
making it available with such a non-restrictive license.

*/

class ExpressionEvaluator
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

          ExpressionEvaluator();  //< Constructor.
 virtual ~ExpressionEvaluator();  //< Destructor.


 //---------------------------------------------------------------------------
 // parse/evaluate expressions:

 virtual bool parseExpression(char* expression);
 /**< Parses an expression. The parsing is separated from the actual 
      calculation of the result to make it possible to re-assign
      variables to different values and then re-calculate the result without
      needing to parse the expression again. The boolean return value 
      indicates if the parsing was successful. Expressions must end with a
      semicolon. */

 virtual bool evaluateExpression(double* result);
 /**< Evaluates the result of the expression which has been parsed before. */

 //---------------------------------------------------------------------------
 // variable assignments:

 virtual bool assignVariable(char* name, double value);
 /**< Assigns a variable name to a numeric vaule. */

 virtual bool getVariableAdress(char* name, double* adress);
 /**< Returns the memory adress of a variable. Can be used to access the value
      of this variable faster than via calling assignVariable(). */

 //---------------------------------------------------------------------------
 // others:

 virtual bool clearExpression();
 /**< Has to be called, before a new expression string is to be parsed. */

 virtual int getLastError();
 /**< Returns the error code of the most recently performed action. */


 //===========================================================================

protected:

 exprFuncList* pFunctionList;      // pointer to the function list
 exprValList*  pVariableList;      // pointer to the variable list
 exprValList*  pConstantList;      // pointer to the constant list
 exprObj*      pExpressionObject;  // pointer to an expression object

 EXPRTYPE  result;  // variable to store the result,
                    // EXPRTYPE is a typedef for double

 int err;  // variable to hold an error code

};


#endif // ExpressionEvaluator_h