#ifndef RAPT_STATISTICS_H_INCLUDED
#define RAPT_STATISTICS_H_INCLUDED

/** A collection of functions for statistical data analysis. */

class Statistics
{

public:

  /** Given two input arrays x and y of length N where x represents abscissa values and y the 
  corresponding ordinate values, the function computes straight line equations coefficient for
  a line y = a*x + b that best fits the data in a least squares sense. */
  template<class T>
  static void linearRegression(int N, T* x, T* y, T& a, T& b);

};

#endif
