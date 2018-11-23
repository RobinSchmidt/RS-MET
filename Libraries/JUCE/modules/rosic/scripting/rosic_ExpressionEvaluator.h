#ifndef rosic_ExpressionEvaluator_h
#define rosic_ExpressionEvaluator_h

// standard library includes:
//#include <math.h>
//#include <string.h>

//#include "rosic_ExpressionEvaluatorFunctions.h"
//#include "rosic_ExpressionEvaluatorComplexFunctions.h"
//#include "../infrastructure/rosic_MutexLock.h"

namespace rosic
{

  /**

  This is a class for parsing and evaluating mathematical expressions. It encapsulates the ExprEval 
  c++ library by Brian Allen Vanderburg II, thereby making it easier to use and augmenting it by 
  additional functions and constants. 

  Much thanks to Brian Allen Vanderburg II for developing the library and
  making it available with such a non-restrictive license.

  */

  class ExpressionEvaluator
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    ExpressionEvaluator();    //< Constructor.
    ~ExpressionEvaluator();   //< Destructor.

    //---------------------------------------------------------------------------------------------
    // set and evaluate expressions:

    bool setExpressionString(const char* newExpressionString);
    /**< Sets a new expression string. This will automatically trigger a parsing of the expression 
    such that subsequent calls to evaluateExpression will give the result of this expression. */

    double evaluateExpression();
    /**< Evaluates the result of the expression which has been set up before. Make sure that you 
    have assigned all identifiers in the expression with numeric values (via assignVariable()) 
    before you call this. */

    //---------------------------------------------------------------------------------------------
    // variable assignments:

    //void initVariableList();
    /** Clears the variable list and thereafter adds the standard constants. */

    void initFunctionList();
    /** Clears the function list and thereafter adds the standard constants. */

    void assignVariable(const char* name, double value);
    /**< Assigns a variable name to a numeric vaule. This function should be called for each 
    identifier which is containend in the expression (before or after an  expression was set). It
    can be called (for each identifier) as many times as you wish to re-assign identifiers to new
    numeric values without needing to re-parse the expression again (parse once, evaluate many). */

    double* getVariableAddress(const char* name);
    /**< Returns the memory adress of a variable. Can be used to access the value of this variable
    faster than via calling assignVariable(). */

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns a pointer to a c-string containing the expression. */
    const char* getExpressionString() const { return expressionString; }

    //=============================================================================================

  protected:

    bool parseExpression();
    /**< Parses an expression. The parsing is separated from the actual calculation of the result
    to make it possible to re-assign variables to different values and then re-calculate the result
    without needing to parse the expression again. The boolean return value indicates if the 
    parsing was successful. Expressions must end with a semicolon. */

    ExprEval::ValueList    valueList;     // list of values (variables and constants
    ExprEval::FunctionList functionList;  // list of functions
    ExprEval::Expression   expression;    // the expression object

    char* expressionString;     // a c-string containing the expression

    MutexLock mutex;

  };

} // end namespace rosic

#endif // rosic_ExpressionEvaluator_h