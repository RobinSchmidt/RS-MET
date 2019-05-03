#ifndef rosic_PiecewiseFunction_h
#define rosic_PiecewiseFunction_h

//// includes from the STL:
//#include <iterator>
//#include <vector>
//using std::vector;
//
//// rosic-indcludes:
//#include "../math/rosic_ElementaryFunctionsReal.h"
//#include "../math/rosic_Interpolation.h"
//#include "../math/rosic_LinearAlgebra.h"
//#include "../infrastructure/rosic_MutexLock.h"

namespace rosic
{

  /** A class to combine all the required data for one function breakpoint. */

  class FunctionBreakpoint
  {

    friend class PiecewiseFunction;

  public:

    /** Contructor. */
    FunctionBreakpoint(double x = 0.0, double y = 0.0)
    {
      this->x = x;  
      this->y = y;
    }

  protected:

    // member variables:
    double x, y;

  };

  /**

  This is a class that defines a piecewise function by interpolating between a set of user defined 
  breakpoints via various interpolation methods.

  References: 
   -(1) Bronstein et al: Taschenbuch der Mathematik

  \todo:
  -refine the segment-search - binary search and taking the segment from the previous call as 
   initial guess
  -allow for different boundary conditions for the cubic spline (letting the user choose values for 
   the first derivative)

  */

  class PiecewiseFunction
  {

  public:

    enum interpolationMethods
    {
      CONSTANT = 0,
      LINEAR,
      CUBIC,

      NUM_INTERPOLATION_METHODS
    };

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    PiecewiseFunction();

    /** Destructor. */
    ~PiecewiseFunction();

    /** Copies the data (the content of all member variables) from the passed source into this
    instance of PiecewiseFunction. */
    void copyDataFrom(PiecewiseFunction& source);

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Selects one of the interpolation methods. @see interpolationMethods */
    void setInterpolationMethod(int newInterpolationMethod);

    /** Sets a minimum for the y-values */
    void setMinimumAllowedY(double newMinY);

    /** Sets a maximum for the y-values */
    void setMaximumAllowedY(double newMaxY);

    /** Inserts a new breakpoint into the vector. The new breakpoint must satisfy some
    preconditions in order to be successfully inserted (it is not allowed to be too close to an
    already existing breakpoint, for example). The integer return-value informs, at which index if
    the new breakpoint was inserted. It will return -1, when the breakpoint could not be
    inserted. */
    int  insertBreakpoint(double newX, double newY);

    /** Removes a breakpoint from the vector. The return-value informs, if there was actually a
    breakpoint removed (if you try to remove a non-existing breakpoint, or a breakpoint which
    cannot be removed, it will return false */
    bool removeBreakpoint(int index);

    /** Modifies the data of an existing breakpoint. If you try to modify a non-existent breakpoint
    or try to modify the data in such a way which is not allowed, it will return false. */
    bool modifyBreakpoint(int index, double newX, double newY);

    /** Changes the x-coordinate of one breakpoint and reports about the success of that 
    operation. */
    bool setBreakpointX(int index, double newX);

    /** Changes the y-coordinate of one breakpoint and reports about the success of that 
    operation. */
    bool setBreakpointY(int index, double newY);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the index of the last breakpoint (the end-point). */
    int lastBreakpointIndex() const;

    /** Returns the number of data->breakpoints. */
    int getNumBreakpoints() const;

    /** Returns the x-coordinate of the first breakpoint. */
    double getMinX();

    /** Returns the x-coordinate of the last breakpoint. */
    double getMaxX();

    /** Returns the minimum y-coordinate of the function. */
    double getMinY();

    /** Returns the maximum y-coordinate of the function. */
    double getMaxY();

    /** Returns the x-value of one breakpoint. If the index is out of range, it will return NaN. */
    double getBreakpointX(int index) const;

    /** Returns the minimum value to which the x-coordinate of a breakpoint can be set without
    violating the constraints to not come too close to its neighbours. If the index is out of
    range, it will return Nan. */
    double getBreakpointMinX(int index) const;

    /** Returns the maximum value to which the x-coordinate of a breakpoint can be set without
    violating the constraints to not come too close to its neighbours. If the index is out of
    range, it will return -1.0. */
    double getBreakpointMaxX(int index) const;

    /** Returns the y-coordinate of one breakpoint. If the index is out of range, it will return 
    NaN. */
    double getBreakpointY(int index) const;

    //---------------------------------------------------------------------------------------------
    // function value retrieval:

    /** Calculates the y-value corresponding to the given x-value according to the defined 
    breakpoints and the chosen interpolation method. */
    INLINE double getValue(double x);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Sets the object to an initial state with only two breakpoints at (-1,-1) and (1,1) */
    void initialize();

    //=============================================================================================

  protected:

    /** Returns the index of the next breakpoint which has a x-coordinate strictly larger then the
    x-coordinate of the startIndex - mostly this will be just one number higher than startIndex, but
    when there are more breakpoints at one time instant, it will be higher than that. If the result
    is out of range, it will return Nan. */
    int getNextNonCoincidingIndex(int startIndex);

    /** Return the last index at the same x-coordinate as startIndex - if there is only one 
    breakpoint at this x-coordinate, this will return the startIndex iteself. */
    int getLastCoincidingIndex(int startIndex);

    /** Calculates the coefficients for the cubic interpolation spline. */
    void calculateCubicCoefficients();

    /** Make a copy-constructor unavailable because this class needs deep copying (because of the
    pointers in the MutexLocks). If you need to create a new object based on an existing one,
    simply create a new object with the default constructor and then use copyDataFrom(). */
    PiecewiseFunction(const PiecewiseFunction& piecewiseFunctionToCopy);

    //---------------------------------------------------------------------------------------------

    /** The array of the breakpoints. */
    std::vector<FunctionBreakpoint> breakpoints;

    /** A MutexLock for accessing the array of breakpoints. */
    MutexLock mutex;

    double minY, maxY;
    double minBreakpointDistance;
    int    interpolationMethod;

    double *a, *b, *c, *d;
    int    numSegments;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE double PiecewiseFunction::getValue(double x)
  {
    //double xL, xR, yL, yR; // x- and y-coordinates of the breakpoints left and right to x
    double y;

    mutex.lock();    // access the array of breakpoints only inside mutex locks...
    if( x < breakpoints[0].x )
    {
      y = breakpoints[0].y;
      mutex.unlock();
      return y;
    }
    else if( x > breakpoints[lastBreakpointIndex()].x )
    {
      y = breakpoints[lastBreakpointIndex()].y;
      mutex.unlock();
      return y;
    }
    else
    {
      // find the right segment (\todo: use binary search):
      int segmentEndIndex = 1;
      while( segmentEndIndex <= lastBreakpointIndex() && breakpoints[segmentEndIndex].x <= x )
        segmentEndIndex++;
      int segmentStartIndex = segmentEndIndex-1;

      switch( interpolationMethod )
      {
      case CONSTANT: 
        {
          y = breakpoints[segmentStartIndex].y;
        }
        break;
      case LINEAR:   
        {
          // x- and y-coordinates of the breakpoints left and right to x:
          double xL = breakpoints[segmentStartIndex].x;
          double yL = breakpoints[segmentStartIndex].y;
          double xR = breakpoints[segmentEndIndex].x;
          double yR = breakpoints[segmentEndIndex].y;
          y         = yL + (x-xL)*(yR-yL)/(xR-xL);   
        }
        break;
      case CUBIC:   
        {
          int i     = segmentStartIndex;
          if( i < 0 || i >= numSegments )
          {
            DEBUG_BREAK; // this means that the numSegments and breakpoints.size() are not 
                         // consistent ...some threading/mutexing problem?
            y = 0.0;
          }
          else
          {
            double x1 = x - breakpoints[i].x;                // (x-x_i)
            double x2 = x1*x1;                               // (x-x_i)^2
            double x3 = x2*x1;                               // (x-x_i)^3
            y         = a[i] + b[i]*x1 + c[i]*x2 + d[i]*x3;  // Eq.19.235
          }
        }
        break;

      default: y = breakpoints[segmentStartIndex].y; // piecewise constant by default
      }
    }
    mutex.unlock();

    return y;  
  }

}  // end namespace rosic

#endif // PiecewiseFunction_h
