#ifndef ExpressionEvaluatorFunctions_h
#define ExpressionEvaluatorFunctions_h

#include "ExprEval\expreval.h"

#include <MoreMath.h>
using namespace MoreMath;

#include <errno.h>
extern int errno;  // is set in functions from math.h and can be checked here

/**

This file defines custom functions (and mathematical constants) for the
ExpressionEvaluator class by using the macros included in the ExprEval
library.

*/



//----------------------------------------------------------------------------
// belowOrAbove function:
int belowOrAbove(struct _exprObj *o, 
                 struct _exprNode *n, 
                 int count,            // number of arguments
                 EXPRTYPE **refitems, 
                 int refcount, 
                 EXPRTYPE *val)        // memory location to store the result
{
	int err;
	EXPRTYPE arg1, arg2, arg3;

	// accept at most 3 arguments:
	if( (count<1) || (count>3) )
		return EXPR_ERROR_BADNUMBERARGUMENTS;

	// evaluate argument 1:
	err = exprEvalNode(o, &(n[0]), &arg1);
	if(err != EXPR_ERROR_NOERROR)
		return err;

	// evaluate argument 2:
 if( count >= 2 )
 {
	 err = exprEvalNode(o, &(n[1]), &arg2);
	 if(err != EXPR_ERROR_NOERROR)
	 	return err;
 }
 else
  arg2 = -1.0;  // default value for the low value

	// evaluate argument 3:
 if( count >= 3 )
 {
	 err = exprEvalNode(o, &(n[2]), &arg3);
	 if(err != EXPR_ERROR_NOERROR)
	 	return err;
 }
 else
  arg3 = 1.0;  // default value for the high value

	// Clear any math routine errors:
	errno = 0;

	// Do the math:
 *val = MoreMath::belowOrAbove(arg1, arg2, arg3); 

	// Check for math errors:
	if(errno == EDOM || errno == ERANGE)
		{
		if(exprGetSoftErrors(o))
			{
			*val = 0.0;
			return EXPR_ERROR_NOERROR;
			}
		else
			{
			return EXPR_ERROR_OUTOFRANGE;
			}
		}

	return EXPR_ERROR_NOERROR;
}

//----------------------------------------------------------------------------
// chebychev function:
int cheby(struct _exprObj *o, 
             struct _exprNode *n, 
             int count,            // number of arguments
             EXPRTYPE **refitems, 
             int refcount, 
             EXPRTYPE *val)        // memory location to store the result
{
	int err;
	EXPRTYPE arg1, arg2;

	// accept at most 2 arguments:
	if( (count<1) || (count>2) )
		return EXPR_ERROR_BADNUMBERARGUMENTS;

	// evaluate argument 1:
	err = exprEvalNode(o, &(n[0]), &arg1);
	if(err != EXPR_ERROR_NOERROR)
		return err;

	// evaluate argument 2:
 if( count >= 2 )
 {
	 err = exprEvalNode(o, &(n[1]), &arg2);
	 if(err != EXPR_ERROR_NOERROR)
	 	return err;
 }
 else
  arg2 = 1.0;  // default value for the order of the polynomial

	// Clear any math routine errors:
	errno = 0;

	// Do the math:
 *val = MoreMath::cheby(arg1, arg2); 

	// Check for math errors:
	if(errno == EDOM || errno == ERANGE)
		{
		if(exprGetSoftErrors(o))
			{
			*val = 0.0;
			return EXPR_ERROR_NOERROR;
			}
		else
			{
			return EXPR_ERROR_OUTOFRANGE;
			}
		}

	return EXPR_ERROR_NOERROR;
}

//----------------------------------------------------------------------------
// foldOver function:
int foldOver(struct _exprObj *o, 
             struct _exprNode *n, 
             int count,            // number of arguments
             EXPRTYPE **refitems, 
             int refcount, 
             EXPRTYPE *val)        // memory location to store the result
{
	int err;
	EXPRTYPE arg1, arg2, arg3;

	// accept at most 3 arguments:
	if( (count<1) || (count>3) )
		return EXPR_ERROR_BADNUMBERARGUMENTS;

	// evaluate argument 1:
	err = exprEvalNode(o, &(n[0]), &arg1);
	if(err != EXPR_ERROR_NOERROR)
		return err;

	// evaluate argument 2:
 if( count >= 2 )
 {
	 err = exprEvalNode(o, &(n[1]), &arg2);
	 if(err != EXPR_ERROR_NOERROR)
	 	return err;
 }
 else
  arg2 = -1.0;  // default value for the min value

	// evaluate argument 3:
 if( count >= 3 )
 {
	 err = exprEvalNode(o, &(n[2]), &arg3);
	 if(err != EXPR_ERROR_NOERROR)
	 	return err;
 }
 else
  arg3 = 1.0;  // default value for the max value

	// Clear any math routine errors:
	errno = 0;

	// Do the math:
 *val = MoreMath::foldOver(arg1, arg2, arg3); 

	// Check for math errors:
	if(errno == EDOM || errno == ERANGE)
		{
		if(exprGetSoftErrors(o))
			{
			*val = 0.0;
			return EXPR_ERROR_NOERROR;
			}
		else
			{
			return EXPR_ERROR_OUTOFRANGE;
			}
		}

	return EXPR_ERROR_NOERROR;
}


//----------------------------------------------------------------------------
// gaussian function:
int gauss(struct _exprObj *o, 
          struct _exprNode *n, 
          int count,            // number of arguments
          EXPRTYPE **refitems, 
          int refcount, 
          EXPRTYPE *val)        // memory location to store the result
{
	int err;
	EXPRTYPE arg1, arg2, arg3;

	// accept at most 3 arguments:
	if( (count<1) || (count>3) )
		return EXPR_ERROR_BADNUMBERARGUMENTS;

	// evaluate argument 1:
	err = exprEvalNode(o, &(n[0]), &arg1);
	if(err != EXPR_ERROR_NOERROR)
		return err;

	// evaluate argument 2:
 if( count >= 2 )
 {
	 err = exprEvalNode(o, &(n[1]), &arg2);
	 if(err != EXPR_ERROR_NOERROR)
	 	return err;
 }
 else
  arg2 = 0.0;  // default value for the mean mu

	// evaluate argument 3:
 if( count >= 3 )
 {
	 err = exprEvalNode(o, &(n[2]), &arg3);
	 if(err != EXPR_ERROR_NOERROR)
	 	return err;
 }
 else
  arg3 = 1.0;  // default value for the standard deviation sigma

	// make sure, that sigma is not 0.0 (division by zero):
	if(arg3 == 0.0)
	{
		if(exprGetSoftErrors(o)) // Soft errors on?
		{			
			*val = 0.0;  // Yes
			return EXPR_ERROR_NOERROR;
		}
		else // No
		{
			return EXPR_ERROR_OUTOFRANGE; // Bad argument
		}
	}

	// Clear any math routine errors:
	errno = 0;

	// Do the math:
 *val = MoreMath::gauss(arg1, arg2, arg3); 

	// Check for math errors:
	if(errno == EDOM || errno == ERANGE)
		{
		if(exprGetSoftErrors(o))
			{
			*val = 0.0;
			return EXPR_ERROR_NOERROR;
			}
		else
			{
			return EXPR_ERROR_OUTOFRANGE;
			}
		}

	return EXPR_ERROR_NOERROR;
}


//----------------------------------------------------------------------------
// logistic function:
int logistic(struct _exprObj *o, 
             struct _exprNode *n, 
             int count,            // number of arguments
             EXPRTYPE **refitems, 
             int refcount, 
             EXPRTYPE *val)        // memory location to store the result
{
	int err;
	EXPRTYPE arg1, arg2;

	// accept at most 2 arguments:
	if( (count<1) || (count>2) )
		return EXPR_ERROR_BADNUMBERARGUMENTS;

	// evaluate argument 1:
	err = exprEvalNode(o, &(n[0]), &arg1);
	if(err != EXPR_ERROR_NOERROR)
		return err;

	// evaluate argument 2:
 if( count >= 2 )
 {
	 err = exprEvalNode(o, &(n[1]), &arg2);
	 if(err != EXPR_ERROR_NOERROR)
	 	return err;
 }
 else
  arg2 = 1.0;  // default value for the slope beta

	// Clear any math routine errors:
	errno = 0;

	// Do the math:
 *val = MoreMath::logistic(arg1, arg2); 

	// Check for math errors:
	if(errno == EDOM || errno == ERANGE)
		{
		if(exprGetSoftErrors(o))
			{
			*val = 0.0;
			return EXPR_ERROR_NOERROR;
			}
		else
			{
			return EXPR_ERROR_OUTOFRANGE;
			}
		}

	return EXPR_ERROR_NOERROR;
}

//----------------------------------------------------------------------------
// minkowsi distribution function:
int minkowski(struct _exprObj *o, 
              struct _exprNode *n, 
              int count,            // number of arguments
              EXPRTYPE **refitems, 
              int refcount, 
              EXPRTYPE *val)        // memory location to store the result
{
	int err;
	EXPRTYPE arg1, arg2, arg3, arg4;

	// accept at most 4 arguments:
	if( (count<1) || (count>4) )
		return EXPR_ERROR_BADNUMBERARGUMENTS;

	// evaluate argument 1:
	err = exprEvalNode(o, &(n[0]), &arg1);
	if(err != EXPR_ERROR_NOERROR)
		return err;

	// evaluate argument 2:
 if( count >= 2 )
 {
	 err = exprEvalNode(o, &(n[1]), &arg2);
	 if(err != EXPR_ERROR_NOERROR)
	 	return err;
 }
 else
  arg2 = 0.0;  // default value for the mean mu

	// evaluate argument 3:
 if( count >= 3 )
 {
	 err = exprEvalNode(o, &(n[2]), &arg3);
	 if(err != EXPR_ERROR_NOERROR)
	 	return err;
 }
 else
  arg3 = 1.0;  // default value for the standard deviation sigma

	// make sure, that sigma is not 0.0 (division by zero):
	if(arg3 == 0.0)
	{
		if(exprGetSoftErrors(o)) // Soft errors on?
		{			
			*val = 0.0;  // Yes
			return EXPR_ERROR_NOERROR;
		}
		else // No
		{
			return EXPR_ERROR_OUTOFRANGE; // Bad argument
		}
	}

	// evaluate argument 4:
 if( count >= 4 )
 {
	 err = exprEvalNode(o, &(n[3]), &arg4);
	 if(err != EXPR_ERROR_NOERROR)
	 	return err;
 }
 else
  arg4 = 2.0;  // default value for the exponent

	// Clear any math routine errors:
	errno = 0;

	// Do the math:
 *val = MoreMath::minkowski(arg1, arg2, arg3, arg4); 

	// Check for math errors:
	if(errno == EDOM || errno == ERANGE)
		{
		if(exprGetSoftErrors(o))
			{
			*val = 0.0;
			return EXPR_ERROR_NOERROR;
			}
		else
			{
			return EXPR_ERROR_OUTOFRANGE;
			}
		}

	return EXPR_ERROR_NOERROR;
}

//----------------------------------------------------------------------------
// quantization function:
int quant(struct _exprObj *o, 
          struct _exprNode *n, 
          int count,            // number of arguments
          EXPRTYPE **refitems, 
          int refcount, 
          EXPRTYPE *val)        // memory location to store the result
{
	int err;
	EXPRTYPE arg1, arg2;

	// accept at most 2 arguments:
	if( (count<1) || (count>2) )
		return EXPR_ERROR_BADNUMBERARGUMENTS;

	// evaluate argument 1:
	err = exprEvalNode(o, &(n[0]), &arg1);
	if(err != EXPR_ERROR_NOERROR)
		return err;

	// evaluate argument 2:
 if( count >= 2 )
 {
	 err = exprEvalNode(o, &(n[1]), &arg2);
	 if(err != EXPR_ERROR_NOERROR)
	 	return err;
 }
 else
  arg2 = 0.1;  // default value for the quantization interval

	// Clear any math routine errors:
	errno = 0;

	// Do the math:
 *val = MoreMath::quant(arg1, arg2); 

	// Check for math errors:
	if(errno == EDOM || errno == ERANGE)
		{
		if(exprGetSoftErrors(o))
			{
			*val = 0.0;
			return EXPR_ERROR_NOERROR;
			}
		else
			{
			return EXPR_ERROR_OUTOFRANGE;
			}
		}

	return EXPR_ERROR_NOERROR;
}

//----------------------------------------------------------------------------
// quantization to bits function:
int quantToBits(struct _exprObj *o, 
                struct _exprNode *n, 
                int count,            // number of arguments
                EXPRTYPE **refitems, 
                int refcount, 
                EXPRTYPE *val)        // memory location to store the result
{
	int err;
	EXPRTYPE arg1, arg2;

	// accept at most 2 arguments:
	if( (count<1) || (count>2) )
		return EXPR_ERROR_BADNUMBERARGUMENTS;

	// evaluate argument 1:
	err = exprEvalNode(o, &(n[0]), &arg1);
	if(err != EXPR_ERROR_NOERROR)
		return err;

	// evaluate argument 2:
 if( count >= 2 )
 {
	 err = exprEvalNode(o, &(n[1]), &arg2);
	 if(err != EXPR_ERROR_NOERROR)
	 	return err;
 }
 else
  arg2 = 5;  // default value for the bit depth

	// Clear any math routine errors:
	errno = 0;

	// Do the math:
 *val = MoreMath::quantToBits(arg1, arg2); 

	// Check for math errors:
	if(errno == EDOM || errno == ERANGE)
		{
		if(exprGetSoftErrors(o))
			{
			*val = 0.0;
			return EXPR_ERROR_NOERROR;
			}
		else
			{
			return EXPR_ERROR_OUTOFRANGE;
			}
		}

	return EXPR_ERROR_NOERROR;
}

//----------------------------------------------------------------------------
// signum function:
int sign(struct _exprObj *o, 
         struct _exprNode *n, 
         int count,            // number of arguments
         EXPRTYPE **refitems, 
         int refcount, 
         EXPRTYPE *val)        // memory location to store the result
{
	int err;
	EXPRTYPE arg1;

	// accept only 1 argument:
	if( (count<1) || (count>1) )
		return EXPR_ERROR_BADNUMBERARGUMENTS;

	// evaluate argument 1:
	err = exprEvalNode(o, &(n[0]), &arg1);
	if(err != EXPR_ERROR_NOERROR)
		return err;

	// Clear any math routine errors:
	errno = 0;

	// Do the math:
 *val = MoreMath::sign(arg1); 

	// Check for math errors:
	if(errno == EDOM || errno == ERANGE)
		{
		if(exprGetSoftErrors(o))
			{
			*val = 0.0;
			return EXPR_ERROR_NOERROR;
			}
		else
			{
			return EXPR_ERROR_OUTOFRANGE;
			}
		}

	return EXPR_ERROR_NOERROR;
}

//----------------------------------------------------------------------------
// a 2*pi periodic sawtooth wave function:
int sawWave(struct _exprObj *o, 
            struct _exprNode *n, 
            int count,            // number of arguments
            EXPRTYPE **refitems, 
            int refcount, 
            EXPRTYPE *val)        // memory location to store the result
{
	int err;
	EXPRTYPE arg1;

	// accept only 1 argument:
	if( (count<1) || (count>1) )
		return EXPR_ERROR_BADNUMBERARGUMENTS;

	// evaluate argument 1:
	err = exprEvalNode(o, &(n[0]), &arg1);
	if(err != EXPR_ERROR_NOERROR)
		return err;

	// Clear any math routine errors:
	errno = 0;

	// Do the math:
 *val = MoreMath::sawWave(arg1); 

	// Check for math errors:
	if(errno == EDOM || errno == ERANGE)
		{
		if(exprGetSoftErrors(o))
			{
			*val = 0.0;
			return EXPR_ERROR_NOERROR;
			}
		else
			{
			return EXPR_ERROR_OUTOFRANGE;
			}
		}

	return EXPR_ERROR_NOERROR;
}

//----------------------------------------------------------------------------
// a 2*pi periodic square wave function:
int sqrWave(struct _exprObj *o, 
            struct _exprNode *n, 
            int count,            // number of arguments
            EXPRTYPE **refitems, 
            int refcount, 
            EXPRTYPE *val)        // memory location to store the result
{
	int err;
	EXPRTYPE arg1;

	// accept only 1 argument:
	if( (count<1) || (count>1) )
		return EXPR_ERROR_BADNUMBERARGUMENTS;

	// evaluate argument 1:
	err = exprEvalNode(o, &(n[0]), &arg1);
	if(err != EXPR_ERROR_NOERROR)
		return err;

	// Clear any math routine errors:
	errno = 0;

	// Do the math:
 *val = MoreMath::sqrWave(arg1); 

	// Check for math errors:
	if(errno == EDOM || errno == ERANGE)
		{
		if(exprGetSoftErrors(o))
			{
			*val = 0.0;
			return EXPR_ERROR_NOERROR;
			}
		else
			{
			return EXPR_ERROR_OUTOFRANGE;
			}
		}

	return EXPR_ERROR_NOERROR;
}

//----------------------------------------------------------------------------
// a 2*pi periodic triangle wave function:
int triWave(struct _exprObj *o, 
            struct _exprNode *n, 
            int count,            // number of arguments
            EXPRTYPE **refitems, 
            int refcount, 
            EXPRTYPE *val)        // memory location to store the result
{
	int err;
	EXPRTYPE arg1;

	// accept only 1 argument:
	if( (count<1) || (count>1) )
		return EXPR_ERROR_BADNUMBERARGUMENTS;

	// evaluate argument 1:
	err = exprEvalNode(o, &(n[0]), &arg1);
	if(err != EXPR_ERROR_NOERROR)
		return err;

	// Clear any math routine errors:
	errno = 0;

	// Do the math:
 *val = MoreMath::triWave(arg1); 

	// Check for math errors:
	if(errno == EDOM || errno == ERANGE)
		{
		if(exprGetSoftErrors(o))
			{
			*val = 0.0;
			return EXPR_ERROR_NOERROR;
			}
		else
			{
			return EXPR_ERROR_OUTOFRANGE;
			}
		}

	return EXPR_ERROR_NOERROR;
}



#endif  // ExpressionEvaluatorFunctions_h