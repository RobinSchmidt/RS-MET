#ifndef RAPT_STATISTICS_H_INCLUDED
#define RAPT_STATISTICS_H_INCLUDED

/** A collection of functions for statistical data analysis. */

class rsStatistics
{

public:

  /** Given two input arrays x and y of length N where x represents abscissa values and y the 
  corresponding ordinate values, the function computes straight line equation coefficients for
  a line y = a*x + b that best fits the data in a least squares sense. */
  template<class T>
  static void linearRegression(int N, const T* x, const T* y, T& a, T& b);

  /** Like linearRegression, just that we assume the b-value (the y-intercept, offset) to be
  equal to zero. This is appropriate, when a proportionality between x and y is assumed. Also,
  we return the a-coefficient for the line y = a*x in the function's return value here. */
  template<class T>
  static T proportionalRegression(int N, const T* x, const T* y);

  // maybe move these into rsCurveFitter and rename to fitLine, fitLineThroughZero

};

//=================================================================================================

/** A class to remove a linear trend and absoulte offset from data and possibly later re-apply 
it. This is useful for dealing with data that has an intrinsic linear trend, such as unwrapped 
phases over time, but this trend is not wanted in processing algorithms. So the idea is: 
analyze trend (and offset) -> remove it -> do processing -> re-apply trend (and offset). */

template<class T>
class rsDeTrender
{

public:

  /** Analyzes the trend and offset, i.e. the coefficients a,b in t(x) = a*x + b, from input
  data via linear regression. */
  void analyzeTrendAndOffset(int N, const T* x, const T* y)
  {
    rsStatistics::linearRegression(N, x, y, a, b);
  }

  /** Analyzes the trend and offset in the data and removes it, while storing the coefficients
  that have been used in order to re-apply the trend later. If for some reason, the a,b, coeffs are
  still valid from an earlier analysis, the analysia may be skipped by passing false for "analyze".
  The y-axis input is taken from yIn and the output is written to yOut, but these arrays may be the 
  same. */
  void removeTrendAndOffset(int N, const T* x, const T* yIn, T* yOut, bool analyze = true)
  {
    if(analyze)
      analyzeTrendAndOffset(N, x, yIn);
    for(int i = 0; i < N; i++)
      yOut[i] = yIn[i] - a*x[i] - b;
  }

  /** (Re)-applies the trend to the data, according to our stored a,b coefficients from prior 
  analysis. */
  void applyTrendAndOffset(int N, const T* x, const T* yIn, T* yOut)
  {
    for(int i = 0; i < N; i++)
      yOut[i] = yIn[i] + a*x[i] + b;
  }

protected:

  T a = 0, b = 0;

};


#endif
