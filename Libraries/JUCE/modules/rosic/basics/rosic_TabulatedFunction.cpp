#include "rosic_TabulatedFunction.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

TabulatedFunction::TabulatedFunction()
{
  // init the tableSize and allocate memory for the table:
  tableSize  = 1024;
  funcTbl    = NULL;
  funcTbl    = new doubleA[tableSize];
  lowerLimit = -2.0;
  upperLimit =  2.0;
  setRange(lowerLimit, upperLimit);
  setFunctionString("tanh(x);", true);
}

TabulatedFunction::~TabulatedFunction()
{
  if(funcTbl)
    delete[] funcTbl;
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

bool TabulatedFunction::setTableSize(int newTableSize)
{
  mutex.lock();
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
  mutex.unlock();
}

void TabulatedFunction::setRange(double newLowerLimit, double newUpperLimit)
{
  mutex.lock();
  if((newUpperLimit-newLowerLimit) != 0.0) //avoid division by zero
  {
    lowerLimit = newLowerLimit;
    upperLimit = newUpperLimit;
  }
  mappingSlopeRec = (tableSize-1) / (upperLimit-lowerLimit);
  //reciprocal of the slope factor a in the mapping x = i*a + b
  //between the input-argument x and the table-index i
  mutex.unlock();
}

void TabulatedFunction::assignVariable(const char *name, double value, bool reCalculateTable)
{
  evaluator.assignVariable(name, value);
  if( reCalculateTable == true )
    calculateTable();
}

bool TabulatedFunction::setFunctionString(const char *newFunctionString, bool reCalculateTable)
{
  // parses the expression (does not yet evaluate it numerically - this
  // will be done inside the calculateTable()-function for different values
  // of x):
  evalSuccess = evaluator.setExpressionString(newFunctionString);

  if(!evalSuccess)
    return false;
  else
  {
    if( reCalculateTable == true )
      calculateTable();
    return true;
  }
}

void TabulatedFunction::calculateTable()
{
  mutex.lock();

  doubleA x, y;
  doubleA a, b; // x = a*i + b where a and b are determined by
  // the upper and lower limit (mapping between the
  // table-index i and the x-variable)

  a = (upperLimit-lowerLimit)/(tableSize-1); // the slope of the mapping
  b = lowerLimit;                            // the constant factor

  for(int i=0; i<tableSize; i++)
  {
    x = a*i + b;  // map table-index to an x-value

    evaluator.assignVariable("x", x);
    // occurrences of string "x" in newFunctionString will be replaced with
    // the current numeric value of variable x

    y          = evaluator.evaluateExpression();
    funcTbl[i] = y;
  }

  mutex.unlock();
}

void TabulatedFunction::clipTableValues(double minTableValue, double maxTableValue)
{
  mutex.lock();
  for(int i=0; i<tableSize; i++)
    funcTbl[i] = clip(funcTbl[i], minTableValue, maxTableValue);
  mutex.unlock();
}

