#ifndef TabulatedFunction_h
#define TabulatedFunction_h

#include "AudioModule.h"
#include "ExpressionEvaluator.h"

/**

This class tabulates a mathematical function (within a certain range
for input arguments x) and can return approximate y-values for any input 
arguments within this range by means of linear or cubic interpolation. This
can avoid expensive mathematical operations. However, tests showed that 
it doesn't achieve performance gains for standard-functions such as 
sin(), cos(), tanh() and the like - it is recommended to always test, 
if the actual mathematical operation is really more expensive than the 
use of a table.

*/

class TabulatedFunction : public AudioModule
{
public:

 //---------------------------------------------------------------------------
 // construction/destruction:

	         TabulatedFunction();  ///< Constructor.
	virtual ~TabulatedFunction();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 virtual bool setTableSize(int newTableSize);
 /**< Sets the size (number of values) of the table. The boolean return 
      value indicates if the memory allocation was successful. */

 virtual void setRange(double newLowerLimit, double newUpperLimit);
 /**< Set the supposed range of input arguments - if an argument out of this
      range is passed, the first or last value of the table will be 
      returned. */

 virtual bool assignVariable(char *name, double value);
 /**< Assigns a variable name to a numeric vaule. */

 virtual bool setFunctionString(char *newFunctionString);      
 /**< Accepts a string like "sin(x)", "tanh(x)", ...
      and fills the table with this fucntion. The boolean result indicates,
      if the operation was successful (if the string is valid).  */
                                                                      
 //virtual void setTable(double *newTable);  
 /**< The table passed in the argument will be copied into the internal
      function table (use it for special functions which can not easily
      be expressed with terms like above */  

 virtual void calculateTable();
 /**< Calculates the table from the function string. */ 
 
 /*   This function should be
      called each time after a new function string is passed or the value of
      some variable has changed. It is not called automatically from 
      setFunction() or assignVariable().  ...at the moment it is called
      automatically ....*/


 //---------------------------------------------------------------------------
 // inquiry:

 virtual int getTableSize() { return tableSize; }
 /**< Returns the size (the number of samples) of the table */

 (double*) getTableAdress() { return funcTbl; }
 /**< Returns the adress (pointer to double) where the table begins. */

 //---------------------------------------------------------------------------
 // table readout:

 virtual INLINE double getValueLinear(double x); 
 /**< Calculates the result by means of linear interpolation. */

 //virtual INLINE double getValueCubic(double x); 
 /**< Calculates the result by means of cubic spline interpolation. */

 //---------------------------------------------------------------------------
 // embedded objects:
 ExpressionEvaluator evaluator;

 //===========================================================================

protected:

 int  tableSize;
 bool evalSuccess;

 //parameter variables:
 doubleA lowerLimit, upperLimit; // upper and lower limit for input arguments
 doubleA mappingSlopeRec;

 doubleA tableLengthMinus1Rec;   // reciprocal of the tableSize-1 

 //the actual table:
 doubleA* funcTbl;               // the actual lookup table

};

//----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):

INLINE double TabulatedFunction::getValueLinear(double x)
{
 static doubleA tblPos;     // postion inside the table where the value is read
 static intA tblPosInt;  // integer part of the position
 static doubleA tblPosFrac; // fractional part of the position

 //determine location in the table:
 tblPos = (x-lowerLimit) * mappingSlopeRec;

 //split position into integer and fractional part:
 tblPosInt  = (int) floor(tblPos);
 tblPosFrac = tblPos - tblPosInt;

 if ( tblPosInt >= (tableSize-1) )    
  return funcTbl[(tableSize-1)]; //saturate for x >= upperLimit
 else if( tblPosInt < 0 )
  return funcTbl[0];
 else
 return (   funcTbl[tblPosInt] 
          + tblPosFrac*(funcTbl[tblPosInt+1]
          - funcTbl[tblPosInt]) );
  // read out value from the table by means of linear interpolation
}


#endif // TabulatedFunction_h
