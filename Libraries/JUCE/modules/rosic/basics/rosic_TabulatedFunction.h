#ifndef rosic_TabulatedFunction_h
#define rosic_TabulatedFunction_h

//// rosic-indcludes:
//#include "../scripting/rosic_ExpressionEvaluator.h"
//#include "../infrastructure/rosic_MutexLock.h"

namespace rosic
{

  /**

  This class tabulates a mathematical function (within a certain range for input arguments x) and
  can return approximate y-values for any input arguments within this range by means of linear or
  cubic interpolation. This can avoid expensive mathematical operations. However, tests showed that
  it doesn't achieve performance gains for standard-functions such as sin(), cos(), tanh() and the
  like - it is recommended to always test, if the actual mathematical operation is really more
  expensive than the use of a table.

  \todo: don't re-allocate memory in setTableSize. instead use a fixed maxsize that is passed to
  the constructor -> avoids mutex-lock stuff.

  */

  class TabulatedFunction
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    TabulatedFunction();

    /** Destructor. */
    ~TabulatedFunction();

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the size (number of values) of the table. The boolean return value indicates if the
    memory allocation was successful. */
    bool setTableSize(int newTableSize);

    /** Set the supposed range of input arguments - if an argument out of this range is passed, the
    first or last value of the table will be returned. */
    void setRange(double newLowerLimit, double newUpperLimit);

    /** Assigns a variable name to a numeric vaule and optionally re-calculates the table. If
    false is passed as second argument, you should call calculateTable yorself manually at some
    later stage. */
    void assignVariable(const char *name, double value, bool reCalculateTable);

    /** Accepts a string like "sin(x)", "tanh(x)", and fills the table with this fucntion. The
    boolean result indicates, if the operation was successful (if the string is valid). If
    false is passed as second argument, you should call calculateTable yorself manually at some
    later stage. */
    bool setFunctionString(const char *newFunctionString, bool reCalculateTable);

    /** The table passed in the argument will be copied into the internal
    function table (use it for special functions which can not easily
    be expressed with terms like above */
    //void setTable(double *newTable);

    /** Calculates the table from the function string. This function should be called each time
    after a new function string is passed or the value of some variable has changed. It is not
    called automatically from setFunction() or assignVariable().  ...at the moment it is called
    automatically ....*/
    void calculateTable();

    /** Clips all values which are currently stored in the table to the range between the passed
    minTableValue and maxTableValue. Useful to apply after calculateTable() for functions which
    create very large negative and/or positive values. */
    void clipTableValues(double minTableValue, double maxTableValue);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the size (the number of samples) of the table */
    int getTableSize() const { return tableSize; }

    /** Returns the adress (pointer to double) where the table begins. */
    double* getTableAdress() const { return funcTbl; }

    /** Returns the lower limit for the input variable. */
    double getLowerLimitX() const { return lowerLimit; };

    /** Returns the upper limit for the input variable. */
    double getUpperLimitX() const { return upperLimit; } ;

    /** Returns a pointer to a c-string containing the expression. */
    const char* getExpressionString() const { return evaluator.getExpressionString(); }

    //---------------------------------------------------------------------------------------------
    // table readout:

    /** Aquires the mutex lock for the table. This is supposed to be used, if you want to do
    faster loockup inside some tight loop - you can then acquire and release the lock manually
    outside the tight loop and use getValueLinearNoLock instead of getValueLinear inside the
    tight loop. */
    INLINE void acquireLock();

    /** Releases the mutex lock for the table. */
    INLINE void releaseLock();

    /** Calculates the result by means of linear interpolation. */
    INLINE double getValueLinear(double x);

    /** Calculates the result by means of linear interpolation without aquiring the mutex lock.
    You should hold the lock in you calling function when using this function.
    @see acquireLock(), releaseLock() */
    INLINE double getValueLinearNoLock(double x);

    /** Calculates the result by means of cubic spline interpolation. */
    //INLINE double getValueCubic(double x);

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

    //doubleA tableLengthMinus1Rec;   // reciprocal of the tableSize-1

    //the actual table:
    doubleA* funcTbl;               // the actual lookup table

    MutexLock mutex; // mutex to avoid retrieval of values during tabel-generation

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void TabulatedFunction::acquireLock()
  {
    mutex.lock();
  }

  INLINE void TabulatedFunction::releaseLock()
  {
    mutex.unlock();
  }

  INLINE double TabulatedFunction::getValueLinear(double x)
  {
    double result;
    mutex.lock();
    result = getValueLinearNoLock(x);
    mutex.unlock();
    return result;
  }

  INLINE double TabulatedFunction::getValueLinearNoLock(double x)
  {
    static doubleA tblPos;     // postion inside the table where the value is read
    static intA tblPosInt;  // integer part of the position
    static doubleA tblPosFrac; // fractional part of the position

    //determine location in the table:
    tblPos = (x-lowerLimit) * mappingSlopeRec;

    //split position into integer and fractional part:
    //tblPosInt  = (int) floor(tblPos);
    tblPosInt  = floorInt(tblPos);
    tblPosFrac = tblPos - tblPosInt;

   if ( tblPosInt >= (tableSize-1) )
      return funcTbl[(tableSize-1)]; // saturate for x >= upperLimit
    else if( tblPosInt < 0 )
      return funcTbl[0];
    else
      // read out value from the table by means of linear interpolation:
      return ( funcTbl[tblPosInt] + tblPosFrac*(funcTbl[tblPosInt+1] - funcTbl[tblPosInt]) );


    /*
    double result;
    mutex.lock();
    if ( tblPosInt >= (tableSize-1) )
      result = funcTbl[(tableSize-1)]; //saturate for x >= upperLimit
    else if( tblPosInt < 0 )
      result = funcTbl[0];
    else
    {
      // read out value from the table by means of linear interpolation:
      result = ( funcTbl[tblPosInt] + tblPosFrac*(funcTbl[tblPosInt+1] - funcTbl[tblPosInt]) );
    }
    mutex.unlock();
    return result;
    */
  }

} // end namespace rosic

#endif // rosic_TabulatedFunction_h
