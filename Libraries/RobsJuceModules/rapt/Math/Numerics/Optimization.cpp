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
T rsMinimizer1D<T>::goldenSectionMin(const std::function<T(T)>& f, T a, T b)
{
  //rsError("Not yet implemented correctly!");
  //return T(0);
  // This algorithm is still under construction and does not yet work properly!

  // https://en.wikipedia.org/wiki/Golden-section_search

  static const int maxIts = 1000;             // Maximum number of iterations.
  T k   = (sqrt(5.) - 1.) * 0.5;              // The golden ratio - maybe rename k to phi.
  T xL  = b - k * (b - a);                    // Left inner evaluation point (I think - verify!)
  T xR  = a + k * (b - a);                    // Right inner evaluattion point (I think - verify!)
  T fL  = f(xL);                              // Function value at xL
  T fR  = f(xR);                              // Function value at xR
  T eps = std::numeric_limits<T>::epsilon();  // Not sure, if that value is always appropriate.
  int its = 0;
  while(b - a > eps && its < maxIts)
  {
    if(fL < fR) {
      b  = xR;
      xR = xL;
      fR = fL;
      xL = b - k * (b - a);
      fL = f(xL); }
    else {
      a  = xL;
      xL = xR;
      fL = fR;
      xR = a + k * (b - a);
      fR = f(xR); }
    its++;
  }
  return (a + b) / 2.;

  // The implementation here is based on the code posted here
  // https://stackoverflow.com/questions/21144309/method-of-the-golden-ratio
  // but modified to:
  // (1) Accept a std::function (by const-ref) instead of a function pointer
  // (2) Limit the number of iterations by having the 2nd condition "its < maxIts"
  // (3) Using only one evaluation of f per iteration by caching fL and fR
  // (4) Find the minimum instead of maximum by replacing ">" by "<" in "if(fL < fR)"
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