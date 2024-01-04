template<class T>
void rsMinSqrDifFixSum(T* v, int N, T* s, T* w)
{
  rsAssert(N >= 2, "Input array too short in rsMinSqrDifFixSum");

  typedef std::vector<T> Vec;
  int Nv = N;        // number of values in v
  int Ns = N - 1;    // number of sums in s
  int Nm = Nv + Ns;  // number of linear equations, matrix size

  // handle special case for just one given sum:
  if(Ns == 1) {
    v[0] = v[1] = 0.5*s[0];
    return;
  }

  // establish the diagonals for the coefficient matrix:
  Vec d0(Nm), d1(Nm-1), d2(Nm-2);
  int i;
  d0[0] = 2;
  for(i = 1; i < Nm; i += 2) {
    d0[i]   = 0;
    d0[i+1] = 4; }
  d0[Nm-1] = 2;
  for(i = 0; i < Nm-1; i++)
    d1[i] = 1;
  for(i = 0; i < Nm-3; i += 2) {
    d2[i]   = -2;
    d2[i+1] =  0; }
  d2[i] = -2;

  // apply error weights, if desired:
  if( w != nullptr ) {
    for(i = 0; i < Ns; i++)            // modify outer sub- and superdiagonal
      d2[2*i] *= w[i];
    d0[0]    *= w[0];                  // modify main diagonal
    d0[Nm-1] *= w[Ns-1];
    for(i = 1; i < Ns; i++)
      d0[2*i] *= 0.5*(w[i-1]+w[i]);
  }

  // establish right-hand side vector:
  Vec b(Nm);
  int j = 0;
  for(i = 0; i <= Nm-2; i += 2) {
    b[i]   = 0;
    b[i+1] = s[j];
    j++;
  }
  b[Nm-1] = 0;

  // use temporaries, because things get messed up in the solver:
  Vec bt = b, l2 = d2, l1 = d1, d = d0, u1 = d1, u2 = d2;
  Vec x(Nm); 
  rsSolvePentaDiagonalSystem(&l2[0], &l1[0], &d[0], &u1[0], &u2[0], &bt[0], &x[0], Nm);

  // extract output array v (in x, the outputs are interleaved with the Lagrange multipliers):
  for(i = 0; i < Nv; i++)
    v[i] = x[2*i];
}


template<class T>
T rsMinimizer1D<T>::bisection(const std::function<T(T)>& f, T xL, T xR)
{
  rsError("Not yet implemented correctly!");
  // This algorithm is still under construction and does not yet work properly!

  // Idea:
  // -At any stage of the algorithm we have 3 points xL, xM, xR which are the left, middle and 
  //  right point of the current interval inside which we wnat to find a local minimum.
  // -If f(xL) < f(xM) and f(xM) < f(xR) we assume the function is monotonically increasing and 
  //  therefore the minimum is at the left boundary. Likewise, if f(xL) > f(xM) && f(xM) > f(xR) 
  //  then f is monotonically decreasing and therefore the minimum is at the right boundary.
  // -If none of the above two cases occur, then the minimum is either in the left half interval
  //  xL..xM or the right half interval xM..xR.
  // -At any stage we select either the left half interval xL..xM or the right half interval xM..xR 
  //  for further search refinement

  static const int maxIts = 100;

  T xM = 0.5 * (xL + xR); // Midpoint
  T dx = xR - xL;         // Interval length
  T fL = f(xL);           // left function value
  T fM = f(xM);           // middle function value
  T fR = f(xR);           // right function value


  T s  = std::numeric_limits<T>::epsilon();  // scaler for the midpoint
  // ToDo: verify, if that is a good value. It's just a first quick and dirty guess

  // As convergence criterion, we use the relative interval length:
  int its = 0;
  while(rsAbs(dx) > s * rsAbs(xM) && its < maxIts)
  {
    if(fL < fM && fM < fR) return xL; // minimum found at left boundary
    if(fL > fM && fM > fR) return xR; // minimum found at right boundary

    // At this point, we know that the minimum is somewhere between xL and xR. We figure out on 
    // which side of xM it is and update xL or xR accordingly along with the corresponding function
    // values fL, fR:
    if(fL < fM)
    {
      // Minimum is left to xM:
      xR = xM;
      fR = fM;
    }
    else
    {
      // Minimum is right to xM:
      xL = xM;
      fL = fM;
    }

    xM = 0.5 * (xL + xR);   // New midpoint
    dx = xR - xL;           // New interval length
    fM = f(xM);             // New function value at midpoint
    its++;
  }

  return xM;



   
  //return rsMin(bisection(xL, xM), bisection(xM, xR));
    // This recursive implementation illustrates the idea and should produce the right result but 
    // has exponential complexity so it's untenable. The actual code should implement the same idea


  //return rsMin(bisection(xL, xM), bisection(xM, xR));

  // 


  //return rsMin(xL, xM, xR);
}


//=================================================================================================
/*                                          Ideas

Implement this:
https://www.youtube.com/watch?v=j29rVHCpRUY  Minimax Approximation and the Exchange Algorithm
-For the initial guess, maybe start with equidistant or - better - Chebychev nodes or somehow use a
 least-squares solution and find the points of max error of it.
-Finding the max-error position(s) may itself be a root-finding problem (a case for 
 Newton-iteration?)
-Maybe exchanging multiple reference points for actual max-error points per step at once could lead
 to faster overall convergence? Experiment a bit with that.


This:
https://www.youtube.com/watch?v=cLtN6jAG_3Q  A Review of Top 16 Optimizers for Training Neural Networks
explains a lot of variations of the gradient descent method.




*/