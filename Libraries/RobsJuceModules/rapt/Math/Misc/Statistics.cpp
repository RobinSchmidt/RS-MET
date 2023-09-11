template<class T>
void rsStatistics::linearRegression(int N, const T* x, const T* y, T& a, T& b)
{
  T xm = rsArrayTools::mean(x, N);
  T ym = rsArrayTools::mean(y, N);
  T xx = rsArrayTools::sumOfSquares(x, N);
  T xy = rsArrayTools::sumOfProducts(x, y, N);
  a = (xy - N*xm*ym) / (xx - N*xm*xm);
  b = ym - a*xm;
}
// todo: use pointers for output variables

template<class T>
T rsStatistics::proportionalRegression(int N, const T* x, const T* y)
{
  T xx = rsArrayTools::sumOfSquares(x, N);
  T xy = rsArrayTools::sumOfProducts(x, y, N);
  return xy / xx;
}

/*

ToDo: 
-Implement multiple linear regression (multiple x-arrays, one y-array)
 linearRegression(int numInputs, int numDataPoints, const T** x, const T* y, T* a, T* b)
 x is numInputs x numDataPoints matrix
-Based on that, implement univariate polynomial regression - it uses, 1,x,x^2,x^3,... as the 
 x-arrays in multiple linear regression - but that should go to CurveFitting.h/cpp

Ideas:
-Fit a line y = a*x + b to model y as function of x, then fit another line x = c*y + d to model y 
 as function of y. Convert the result of x = f(y) into y = f(x) = (x-d)/c = x/c - d/c. By somehow 
 "crossfading" or morphing between the two solutions, the user can decide, which error to minimize
 more strongly - the vertical or the horizontal error. I guess, with an equal mix, we get something
 that minimizes the avarage Euclidean distance between the line and the data points. This technique
 could be related to total least squares: https://en.wikipedia.org/wiki/Total_least_squares  
 ...figure out the details. Maybe implement a proper total least squares fitting, too. This 
 crossfading approach could perhaps be genralized to other functions like general quadratic forms 
 (conics) or linfrac approximants. Such a fitting scheme can be better than regular least squares 
 if we have noise in the x- and y-values. regular least squares assumes exact x-values and noise 
 only in the y-values.
 Q: How do we best morph between two lines a1*x + b1  and  a2*x + b2? We could linearly interpolate
 a1..a2 and b1..b2 separately but maybe that's not the best thing to do? But perhaps it is. If we
 view both lines defined by (a1,b1), (a2,b2) as points in the a-b-plane, the most natural thing to 
 go from (a1,b1) to (a2,b2) would be in a straight line, I guess.



*/