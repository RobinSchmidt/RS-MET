#ifndef RAPT_STATISTICS_H_INCLUDED
#define RAPT_STATISTICS_H_INCLUDED

/** A collection of functions for statistical data analysis. */

class Statistics
{

public:

  /** Given two input arrays x and y of length N where x represents abscissa values and y the 
  corresponding ordinate values, the function computes straight line equation coefficients for
  a line y = a*x + b that best fits the data in a least squares sense. */
  template<class T>
  static void linearRegression(int N, T* x, T* y, T& a, T& b);

  /** Like linearRegression, just that we assume the b-value (the y-intercept, offset) to be
  equal to zero. This is appropriate, when a proportionality between x and y is assumed. Also,
  we return the a-coefficient in the function's return value here. */
  template<class T>
  static T proportionalRegression(int N, T* x, T* y);

};

#endif
