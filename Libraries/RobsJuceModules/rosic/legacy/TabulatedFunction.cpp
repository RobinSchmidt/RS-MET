#include "TabulatedFunction.h"

//----------------------------------------------------------------------------
// construction/destruction:

TabulatedFunction::TabulatedFunction()
{
 // init the tableSize and allocate memory for the table:
 tableSize = 1024;
 funcTbl   = NULL;
 funcTbl   = new doubleA[tableSize];

 bool dummy;

 lowerLimit = -2.0;
 upperLimit =  2.0;

 setRange(lowerLimit, upperLimit);
 dummy = setFunctionString("tanh(x);");
}

TabulatedFunction::~TabulatedFunction()
{
 if(funcTbl)
  delete[] funcTbl;
}

//----------------------------------------------------------------------------
// parameter settings:

bool TabulatedFunction::setTableSize(int newTableSize)
{
 if( newTableSize == tableSize )
  return true;
 else
 {
  // free memory with the old table:
  if(funcTbl)
   delete[] funcTbl;
  funcTbl = NULL;

  // allocate memory for the new table:
  funcTbl = new doubleA[newTableSize];

  // if this memory allocation was successful, set the tableSize to
  // newTableSize, recalculate the table and report success:
  if( funcTbl != NULL )
  {
   tableSize = newTableSize;
   calculateTable();
   return true;
  }
  else
  {
   // try to allocate memory with the old tableSize and report failure:
   funcTbl = new doubleA[tableSize];
   calculateTable();
   return false;
  }

 }

}

void TabulatedFunction::setRange(double newLowerLimit, double newUpperLimit)
{
 if((newUpperLimit-newLowerLimit) != 0.0) //avoid division by zero
 {
  lowerLimit = newLowerLimit;
  upperLimit = newUpperLimit;
 }

 mappingSlopeRec = (tableSize-1) / (upperLimit-lowerLimit);
  //reciprocal of the slope factor a in the mapping x = i*a + b
  //between the input-argument x and the table-index i
}

bool TabulatedFunction::assignVariable(char *name, double value)
{
 bool success = evaluator.assignVariable(name, value);

 if(success)
  calculateTable();

 return success;
}

bool TabulatedFunction::setFunctionString(char *newFunctionString)
{
 // clear the expression in the embedded ExpressionEvaluator object (has 
 // to done each time, a new expression is going to be parsed):
 evaluator.clearExpression();

 // parses the expression (does not yet evaluate it numerically - this
 // will be done inside the calculateTable()-function for different values
 // of x):
 evalSuccess = evaluator.parseExpression(newFunctionString);

 if(!evalSuccess)
  return false;
 else
 {
  calculateTable();
  return true;
 }

}

void TabulatedFunction::calculateTable()
{
 doubleA x;
 doubleA a, b; // x = a*i + b where a and b are determined by 
              // the upper and lower limit (mapping between the 
              // table-index i and the x-variable)

 a = (upperLimit-lowerLimit)/(tableSize-1); // the slope of the mapping
 b = lowerLimit;                            // the constant factor

 for(int i=0; i<tableSize; i++)
 {
  x = a*i + b;  //map table-index to an x-value

  evaluator.assignVariable("x", x); 
   // occurrences of string "x" in newFunctionString will be replaced with
   // the current numeric value of variable x

  evaluator.evaluateExpression( &(funcTbl[i]) );
   // evaluates the expression and stores the result at position i in 
   // the table
 }
}

