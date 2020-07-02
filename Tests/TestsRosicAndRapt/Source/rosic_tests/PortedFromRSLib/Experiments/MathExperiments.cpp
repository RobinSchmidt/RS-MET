#include "MathExperiments.h"

typedef std::complex<double> rsComplexDbl; // maybe get rid of that



void bandMatrix()
{
  // under construction - we want an algorithm for solving band-diagonal linear systems - ideally 
  // without requiring extra storage 

  static const int N = 8;    // size of the matrix: NxN
  int nu = 2;                // number of superdiagonals (upper, right diagonals)
  int nl = 3;                // number of subdiagonals (lower, left digonals);
  double A[N][N];            // the matrix
  double b[N];               // right hand side
  //double x[N];               // solution vector

  // fill the matrix with zeros:
  int i, j, k;
  for(i = 0; i < N; i++)
    for(j = 0; j < N; j++)
      A[i][j] = 0;

  // fill the diagonals:
  for(i = 0; i < N; i++) {
    int js = RAPT::rsClip(i - nu, 0, N-1);  // start value for colmun index j
    int je = RAPT::rsClip(i + nu, 0, N-1);  // end value for column index j
    for(j = js; j <= je; j++)               // inner loop has (at most) nu + nl + 1 iterations
      A[i][j] = 10*(i+1) + j+1;
  }

  // fill the right hand side vector b:
  for(i = 0; i < N; i++)
    b[i] = 3*(i+1);  // whatever...some values

  // Gaussian elimination where the subtraction of rows is restricted to only those rows which are
  // not zero already:
  int p[N];       // index of pivot row
  double pivAbs;  // absolute value of pivot element
  for(i = 0; i < N; i++) {
    // search for pivot row:
    p[i] = i;     // preliminary
    pivAbs = 0;
    for(j = i; j < i+nl; j++) {
      if(fabs(A[j][i]) > pivAbs) {
        pivAbs = fabs(A[j][i]);
        p[i] = j;
      }
    }

    // subtract pivot row from rows below it (this may still be wrong):
    int je = RAPT::rsClip(i + nu, 0, N-1);
    for(j = i+1; j <= je; j++)
    {
      if(j == p[i])
        continue;
      double multiplier = A[j][i] / A[p[i]][i];
      b[j] -= multiplier * b[p[i]];
      for(k = i; k < N; k++)   // can this loop also be restricted?
        A[j][k] -= multiplier * A[p[i]][k];
    }

    // ....
  }

  // ...this is still under construction

  // todo: Even though the above Gaussian elemination has only O(N * (nl+nu+1)) time complexity, it
  // still uses O(N^2) memory to store the full matrix. We need a storage scheme that only needs
  // O(N * (nl+nu+1)) as well. maybe make a class rsBandMatrix where we my conveniently access 
  // elements via A(i,j) - the () operator "translates" the indices i,j to actual memory locations
  // of the respective storage format (and correctly returns zero, if the combination i,j refers to
  // a position where the matrix element indeed is zero)
}


void pentaDiagnonalMatrix()
{
  // The algorithm is based on the paper "On Solving Pentadiagonal Linear Systems via 
  // Transformations" by KARAWIA, see here: https://arxiv.org/pdf/1409.4802.pdf

  // this is actually obsolete now - there are apparent errors in the paper and i've given up on
  // hunting for less apparent ones or bugs in my code and instead derived an algorithm for 
  // pentadiagonal systems myself - it's implemented in solvePentaDiagonalSystem

  static const int N = 10; // we use the 10x10 matrix in section 4 of the paper
  double d[10] = { 1,2,3,-4,5,6,7,-1,1,8 };       // main diagonal
  double a[9]  = { 2,2,1,5,-7,3,-1,4,5 };         // 1st superdiagonal
  double b[8]  = { 1,5,-2,1,5,2,4,-3 };           // 2nd superdiagonal
  double c[10] = { 0,3,2,1,2,1,2,1,-2,4 };        // 1st subdiagonal
  double e[10] = { 0,0,1,3,1,5,2,2,2,-1 };        // 2nd subdiagonal
  double y[10] = { 8,33,8,24,29,98,99,17,57,108}; // right hand side
  //double y[10] = { 8,33,8,24,29,82,71,17,57,108}; // right hand side
  // in the paper, they filled up the subdiagonals with leading zeros - probably for more 
  // convenient indexing? should we fill the superdiagonals with trailing zeros, too - or are they
  // not referenced?
  // in the paper, the solution vector y is given as rhs in the matrix equation differently than as 
  // vector - apparently, the authors made a mistake



  // the first algorithm "PTRANS-I":
  double mu[10], alpha[10], beta[10], gamma[10], z[10], x[10];
  
  // Steps 1 and 2 just check the matrix for singularity - we assume here that the matrix
  // is nonsingular - actually, we assume even more: that the system is soluble without
  // explicit pivoting (the algo implicitly chooses always the element at [i][i] as pivot)
  
  // Step 3:
  mu[0]    = d[0]; 
  alpha[0] = a[0]/mu[0];
  beta[0]  = b[0]/mu[0];
  z[0]     = y[0]/mu[0];
  
  // Step 4:
  gamma[1] = c[1];
  mu[1]    = d[1] - alpha[0]*gamma[1];
  alpha[1] = (a[1]-beta[0]*gamma[1]) / mu[1];
  beta[1]  = b[1] / mu[1];
  z[1]     = (y[1]-z[0]*gamma[1]) / mu[1];
  
  // Step 5:
  int i;
  for(i = 2; i <= N-3; i++) {
    gamma[i] =  c[i] - alpha[i-2]*e[i];
    mu[i]    =  d[i] - beta[i-2]*e[i] - alpha[i-1]*gamma[i];
    alpha[i] = (a[i] - beta[i-1]*gamma[i]) / mu[i];
    beta[i]  =  b[i] / mu[i];
    z[i]     = (y[i] - z[i-2]*e[i] - z[i-1]*gamma[i]) / mu[i];
  }
  gamma[N-2] =  c[N-2] - alpha[N-4]*e[N-2];
  mu[N-2]    =  d[N-2] - beta[N-4]*e[N-2] - alpha[N-3]*gamma[N-2];
  alpha[N-2] = (a[N-2] - beta[N-3]*gamma[N-2]) / mu[N-2];
  gamma[N-1] =  c[N-1] - alpha[N-3]*e[N-1];
  mu[N-1]    =  d[N-1] - beta[N-3]*e[N-1] - alpha[N-2]*gamma[N-1];
  z[N-2]     = (y[N-2] - z[N-3]*e[N-2] - z[N-3]*gamma[N-2]) / mu[N-2];
  z[N-1]     = (y[N-1] - z[N-2]*e[N-1] - z[N-2]*gamma[N-1]) / mu[N-1];
  
  // Step 6:
  x[N-1] = z[N-1];
  x[N-2] = z[N-2] - alpha[N-2]*x[N-1];
  for(i = N-3; i >= 0; i--)
    x[i] = z[i] - alpha[i]*x[i+1] - beta[i]*x[i+2];

  // the solution vector is supposed to be x = { 1,2,3,4,5,6,7,8,9,10 }  
  // ...nope - that's still wrong...maybe make a unit test
  // the mu-values are correct, though
  

  // the second algorithm "PTRANS-II":
  double sigma[10], phi[10], rho[10], psi[10], w[10];
  
  // Step 3:
  psi[N-1]   = d[N-1];
  sigma[N-1] = c[N-1]/psi[N-1];
  phi[N-1]   = e[N-1]/psi[N-1];
  w[N-1]     = y[N-1]/psi[N-1];
  
  // Step 4:
  rho[N-2]   =  a[N-2];
  psi[N-2]   =  d[N-2] - sigma[N-1]*rho[N-2];
  sigma[N-2] = (c[N-2] - phi[N-1]*rho[N-2]) / psi[N-2];
  phi[N-2]   =  e[N-2] / psi[N-2];
  w[N-2]     = (y[N-2] - w[N-1]*rho[N-2]) / psi[N-2];
  
  // Step 5:
  for(i = N-3; i >= 2; i--) {
    rho[i]   =  a[i] - sigma[i+2]*b[i];
    psi[i]   =  d[i] - phi[i+2]*b[i] - sigma[i+1]*rho[i];
    sigma[i] = (c[i] - phi[i+1]*rho[i]) / psi[i];
    phi[i]   =  e[i] / psi[i];
    w[i]     = (y[i] - w[i+2]*b[i] - w[i+1]*rho[i]) / psi[i];
  }
  rho[1]   =  a[1] - sigma[3]*b[1];
  psi[1]   =  d[1] - phi[3]*b[1] - sigma[2]*rho[1];
  sigma[1] = (c[1] - phi[3]*rho[1]) / psi[1];
  rho[0]   =  a[0] - sigma[2]*b[0];
  psi[0]   =  d[0] - phi[2]*b[0] - sigma[1]*rho[0];
  w[1]     = (y[1] - w[3]*b[1] - w[2]*rho[1]) / psi[1];
  w[0]     = (y[0] - w[2]*b[0] - w[1]*rho[0]) / psi[0];
  
  // Step 6:
  x[0] = w[0];
  x[1] = w[1] - sigma[1]*x[0];
  for(i = 2; i <= N-1; i++)
    x[i] = w[i] - sigma[i]*x[i-1] - phi[i]*x[i-2];

  // the 2nd algo gives also a wrong result - but a different one - double-check everything!
  // psi[0] is wrong

  int dummy = 0;
}


void pentaDiagnonalMatrix2()
{
  typedef std::vector<double> Vec;

  int N = 9;
  Vec d = { 2,0,4,0,4,0,4,0,2 };
  Vec u = { 1,1,1,1,1,1,1,1   };
  Vec v = {-2,0,-2,0,-2,0,-2  };
  Vec l = u, m = v;
  Vec b = { 1,2,3,4,5,6,7,8,9 };

  //int N = 7;
  //Vec d = { 4,4,4,4,4,4,4 };
  //Vec u = { 1,1,1,1,1,1 };
  //Vec v = { 1,1,1,1,1 };
  //Vec l = u, m = v;
  //Vec b = { 1,2,3,4,5,6,7 };

  // create copies of the original vectors to be used and messed up by the algo (V and M are 
  // actually kept intact, but anyway):
  Vec D=d, U=u, V=v, L=l, M=m, B=b;

  // solve the system by Gaussian elimination (without pivoting) and backsubstitution:
  Vec x = solvePentaDiagonalSystem(M, L, D, U, V, B);

  // verify the result by multiplying the original pentadiagonal matrix with the found 
  // solution vector:
  Vec c = pentaDiagMatVecMul(m, l, d, u, v, x); // should be equal to b (up to roundoff)

  // maybe implement an algorithm that does pivot search and compare the precision - with the 9x9
  // matrix, from the freq-estimation problem we seem to already loose 3 decimal digits of 
  // precision which seems not very good considering that in practice we'll have matrices with
  // hundreds or even thousands of entries
  // maybe convert into a regular (non-diagonal) NxN matrix and use the regular Gaussian 
  // elminination solver for comparison - to figure out, if it's worth to implement pivoting for
  // freq-estimation

  int dummy = 0;
}

// for production code later use plain c-arrays..
// N: length of v, v: computed output values, s: desired sums, w: weights for the squared 
// differences (optional, if not given, all weights are unity)
//template<class T>
//void rsMinSqrDifFixSum(int N, T* v, const T* s, const T* w = nullptr)
//{
//
//}
// output array contains sums of adjacent elements of input array
std::vector<double> adjacentSums(const std::vector<double>& v)
{
  int N = (int)v.size();
  std::vector<double> s(N-1);
  for(int i = 0; i < N-1; i++)
    s[i] = v[i] + v[i+1];
  return s;
}
// output array contains differences of adjacent elements of input array
std::vector<double> adjacentDifs(const std::vector<double>& v)
{
  int N = (int)v.size();
  std::vector<double> d(N-1);
  for(int i = 0; i < N-1; i++)
    d[i] = v[i+1] - v[i];
  return d;
}
// (optionally weighted) sum of squared differences
double sqrdDifSum(
  const std::vector<double>& v, 
  const std::vector<double>& w = std::vector<double>() )
{
  std::vector<double> d = adjacentDifs(v);
  double r = 0;
  if(w.empty()) {
    for(size_t i = 0; i < d.size(); i++)
      r += d[i]*d[i];
    return r;
  }
  RAPT::rsAssert(d.size() == w.size());
  for(size_t i = 0; i < d.size(); i++)
    r += w[i]*d[i]*d[i];
  return r;
}
std::vector<double> valuesWithFixedSums(const std::vector<double>& s, double first)
{
  int N = (int)s.size() + 1;
  std::vector<double> v(N);
  v[0] = first;
  for(int i = 1; i < N; i++)
    v[i] = s[i-1] - v[i-1];
  return v;
}
void minSqrdDifsForFixSums()
{
  // Test minimization of the sum of the squares of the differences between adjacent array elements
  // when their sums are given fixed values.

  std::vector<double> v, s, w; // values, sums and weights
  std::vector<double> t, d;    // adjacent sums and differences for verification
  double c;                    // cost as computed by cost function

  s = { 12, 24, 36 };  
  v = rsMinSqrDifFixSum(s);       // v = 4, 8, 16, 20 
  d = adjacentDifs(v);            // d = 4, 8, 4

  s = { 12, 24, 36, 48 };  
  v = rsMinSqrDifFixSum(s);       // v = 3, 9, 15, 21, 27
  d = adjacentDifs(v);            // d = 6, 6,  6,  6
  // with weights:
  w = { 2, 3, 4, 5 };
  v = rsMinSqrDifFixSum(s, w);   // v = 2.57, 9.43, 14.57, 21.43, 26.57
  d = adjacentDifs(v);           // d = 6.86, 5.14,  6.86,  5.14


  s = { 12, 24, 36, 48, 60 };
  w = {  2,  3,  4,  5,  6 };
  v = rsMinSqrDifFixSum(s, w);   // 3.6, 8.4, 15.6, 20.4, 27.6, 32.4
  d = adjacentDifs(v);           // 4.8, 7.2,  4.8,  7.2,  4.8
  t = adjacentSums(v);
  c = sqrdDifSum(v);
  // v = 3.6, 8.4, 15.6, 20.4, 27.6, 32.4 -> d = 4.8, 7.2, 4.8, 7.2, 4.8 -> ?
  // the result seems wrong - the increasing weights should make the differences between adjacent
  // elements go down for later values...instead, they alternate - wtf? ..hmm...or could it be that
  // the weights are indeed without effect, if the given sums are just linearly increasing?
  // the same problem with unit weights:
  v = rsMinSqrDifFixSum(s);      // 3.6, 8.4, 15.6, 20.4, 27.6, 32.4
  d = adjacentDifs(v);           // 4.8, 7.2,  4.8,  7.2,  4.8
  t = adjacentSums(v);
  c = sqrdDifSum(v);
  // gives the same result indeed...i think i really need to move that over to the experiments 
  // section and create a sage notebook that obtains some analytic solutions


  s = { 20 };        
  v = rsMinSqrDifFixSum(s);       // v = 10, 10
  d = adjacentDifs(v);            // d =  0

  s = { 20, 30 };    
  v = rsMinSqrDifFixSum(s);       // v = 7.5, 12.5, 17.5
  d = adjacentDifs(v);            // d = 5,    5

  s = { 20, 30, 40 };  
  v = rsMinSqrDifFixSum(s);       // 8.33, 11.66, 18.33, 21.66
  d = adjacentDifs(v);

  s = { 20, 30, 40, 50 }; 
  v = rsMinSqrDifFixSum(s);       // v = 7.5, 12.5, 17.5, 22.5, 27.5; d = 5,5,5,5
  d = adjacentDifs(v);

  s = { 20, 30, 40, 50, 60 }; 
  v = rsMinSqrDifFixSum(s);       // v = 8, 12, 18, 22, 28, 32, d = 4,6,4,6,4
  d = adjacentDifs(v);
  // with zero weights at the ends:
  w = { 0, 1, 1, 1, 0 }; 
  v = rsMinSqrDifFixSum(s, w);    // NaN!!!! maybe div-by zero due to lack of pivot-search
  d = adjacentDifs(v);
  w = { 0.5, 1, 1, 1, 0.5 }; 
  v = rsMinSqrDifFixSum(s, w);   // 7.5, 12.5, 17.5, 22.5, 27.5, 32.5
  d = adjacentDifs(v);           // 5,    5,    5,    5,    5,    5


  s = { 20, 30, 40, 50, 60, 70 }; 
  v = rsMinSqrDifFixSum(s);       // 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5
  d = adjacentDifs(v);

  s = { 20, 30, 40, 50, 60, 70, 80 }; 
  v = rsMinSqrDifFixSum(s);       // 7.86, 12.14, 17.86, 22.14, 27.86, 32.14, 37.86, 42.14
  d = adjacentDifs(v);            // 4.29,  5.71, ..

  s = { 40, 60, 40, 60 };
  v = rsMinSqrDifFixSum(s);     // 5, 35, 25, 15, 45
  d = adjacentDifs(v); 

  s = { 40, 60, 40, 60, 40 };
  v = rsMinSqrDifFixSum(s);     // 0, 40, 20, 20, 40, 0
  d = adjacentDifs(v); 

  s = { 1010, 990, 1010, 990 };
  v = rsMinSqrDifFixSum(s);     // 520, 490, 500, 510, 480
  d = adjacentDifs(v);

  s = { 1010, 990, 1010, 990, 1010 };
  v = rsMinSqrDifFixSum(s);    // 525, 485, 505, 505, 485, 525
  d = adjacentDifs(v);

  s = { 1010, 990, 1010, 990, 1010, 990 };
  v = rsMinSqrDifFixSum(s);    // 530, 480, 510, 500, 490, 520, 470
  d = adjacentDifs(v);


  s = { 1010, 990, 1010, 990, 1010, 990, 1010 };
  v = rsMinSqrDifFixSum(s);    // 535, 475, 515, 495, 495, 515, 475, 535
  d = adjacentDifs(v);
  c = sqrdDifSum(v);
  rsPlotVector(v);
  //rsPlotVector(d);

  //v = valuesWithFixedSums(s, 534);  // lowering a value < 535 for the first element increases
  v = valuesWithFixedSums(s, 530);  // sum-of squares as well as max-abs-difference
  d = adjacentDifs(v);
  c = sqrdDifSum(v);
  rsPlotVector(v);
  //rsPlotVector(d);

  //v = valuesWithFixedSums(s, 536);   // increaing it has the same effec, but then the max-abs-dif
  v = valuesWithFixedSums(s, 540);  // occurs at the other end of the value array
  d = adjacentDifs(v);
  c = sqrdDifSum(v);
  rsPlotVector(v);
  //rsPlotVector(d);


  // if the sums alternate, there's tendency for the outer values to become more spread out - 
  // that's undesirable - does that really minimize the sum of the squared differences? that
  // seems not really plausible...maybe try to reduce and increase the first value by 1 and 
  // compute all other values according to that new, fixed first value - see, if the total sum of
  // squared differences really increases in both cases - yes, it does - so, with respect to the 
  // first value, we have indeed a minimum - and if the first value is locked in place, then all
  // others follow due to the constraints

  // maybe we should give higher weights to the differences at the ends to counteract this 
  // tendency? ..but why should a difference at the ends get more weight than inner differences?

  // or maybe we should treat it as a one-dimensional optimization problem: use as cost function 
  // the sum of absolute differences...or maybe better: the maximum absolute difference and 
  // minimize that cost function with respect to the first value - maybe the sum of squared
  // differences is a bad cost function...maybe use not the squares of the differences but the
  // 4th power ...or some general even power - higher exponents should give more impact to single
  // high values of the difference and turn the error function towards minimax - oh - but that 
  // doesn't result in a linear system of equations - terms like (v2-v1)^4 become
  // 4*(v2-v1)^3 when taking the partial derivative...hmm...i think, i should do it via nonlinear
  // minimax optimization

  s = { 1010, 990, 1010, 990, 1010, 990, 1010, 990, 1010, 990, 1010, 990, 1010, 990, 1010 };
  //w = {    8,   7,    6,   5,    4,   3,    2,   1,    2,   3,    4,   5,    6,   7,    8 };
  w = {   64,  49,   36,  25,   16,   9,    4,   1,    4,   9,   16,  25,   36,  49,   64 };
  v = rsMinSqrDifFixSum(s, w);
  //v = rsMinSqrDifFixSum(s);
  d = adjacentDifs(v);
  c = sqrdDifSum(v);
  rsPlotVector(v);
  //rsPlotVector(d*d);
  // it seems to make no difference, if we use weights or not...wtf? oh - it seems like this is 
  // true only, if the weight function is symmetrical - maybe in that case, lowering the error at 
  // one end would increase the error at the other in just the right way that the overall cost 
  // remains constant

  int dummy = 0;

  // Observations:

  // -When the sums are limearly increasing (like: s = 20, 30, 40, 50, ... ): 
  //  -when N (length of v) is odd (and length of s is even), we get a constant value for 
  //   differences between adjacent values
  //  -when N is even, the differences alternate
  //  -this alternation seems counterintuitive but gives indeed the smallest sum-of-squares of the
  //   differences - but it actually seems less desirable than a constant value
  //  -it can be fixed, by giving the first and last difference a weight of 1/2
  //  -...maybe we should automatically apply an (additional) weight of 1/2 to the end-values
  //   when N is even?

  // -When s = 12, 24, 36, 48, 60; w = 2, 3, 4, 5, 6, we get the same results with and without
  //  weights - that's somewhat unexpected, too...figure out...

  // -may it have something to do with the total sum of weights and/or the total sum of 
  //  neighbour-sums and/or its relationship to the total sum of values?
  // -maybe derive an "anti-alternation" filter from that: just obtain the neighbour sums from an
  //  original array of values and create the output with the algo from that obtained neigbour-sum
  //  signal


  // when the length of the sum vector s is odd (i.e. length of v is even), the computed 
  // solution looks suboptimal - calculate N=4 case by hand - see, if the matrix looks 
  // different - use v1,v2,v3,v4,w1,w2,w3 it seems, the longer the vector, the less strong the suboptimality...maybe it sort of averages
  // out for longer vectors...figure out

  // no - for s = { 20, 30, 40 } the result v = { 8.33, 11.66, 18.33, 21.66 } gives indeed a lower
  // sum-of-squared-differences than for example v = { 7.5, 12.5, 17.5, 22.5 } - so the math works
  // out correctly. maybe my minimization criterion is not well suited for even N - maybe i should
  // give the two outermost squared differences less weight (maybe 1/2)

  // probably the best and most flexible way to deal with this is to give the user the option to 
  // apply a weight to each squared difference - we will need that anyway for the freq-estimation.
  // is there a meaningful way to apply weights to the constraints as well? or rather some sort of
  // slack-variables maybe instead of requiring the hard constraint 
  // v1 + v2 = s1 -> v1 + v2 - s1 = 0, require k * (v1 + v2 - s1)^2 = min ...but i think, for the 
  // frequency-estimation problem, hard constraints are actually appropriate - but maybe for other
  // problems...wait...that makes the whole thing nonlinear - so screw it

  // todo: make production version of the function that includes weights for the squared 
  // differences, maybe make a version that uses the simple penta-solver and one that uses 
  // rsBandDiagonalSolver  - maybe call the simple one "Fast"

  // to figure out, if the results are really correct, verify a couple of cases for small N 
  // symbolically in sage

  // maybe instead of minimizing the squared differences, it may also make sense to minimize the 
  // squared (numeric) 2nd derivative or a mix of the two or the squared difference of the squared
  // difference - maybe it's worth to experiment with such more elaborate cost functions - using
  // the 2nd derivative would (i think) allow zero cost for linearly increasing (or decreasing) 
  // sums. the same framework/solver can be used, just the matrix would look different
}

double sqrdCurvSum(
  const std::vector<double>& v)
{
  double c = 0;
  for(size_t i = 1; i < v.size()-1; i++) {
    double d = 2*v[i] - v[i-1] - v[i+1];
    c += d*d;
  }
  return c;
}
void minSqrdCurvForFixSums()
{
  // this is similar to minSqrdDifsForFixSums but uses a cost function that minimizes the squared
  // deviations from a straight line, i.e. some sort of curvature or numeric approximation of the 
  // second derivative

  std::vector<double> v, s, w; // values, sums and weights
  double c;                    // cost as computed by cost function
  std::vector<double> u;

  //s = { 1010, 990, 1010, 990, 1010 };
  //v = rsMinSqrCrvFixSum(s);
  //rsPlotVector(v);

  //s = { 1010, 990, 1010, 990, 1010, 990, 1010, 990, 1010, 990, 1010, 990, 1010 };
  s = { 1010, 990, 1010, 990, 1010, 990, 1010, 990, 1010, 990, 1010, 990 };
  v = rsMinSqrCrvFixSum(s);
  u = rsMinSqrDifFixSum(s);

  //rsPlotVectors(u, v);
  // minimizing both cost functions gives the exact same result...at least in this case
  //c = sqrdCurvSum(v);

  // try a sequence of random sums - see if in this case also both error functions produce the same
  // result:
  RAPT::rsArrayTools::fillWithRandomValues(&s[0], (int) s.size(), 10.0, 20.0, 0);
  v = rsMinSqrCrvFixSum(s);
  u = rsMinSqrDifFixSum(s);
  c = sqrdCurvSum(v);
  rsPlotVectors(u, v, u-v);
  // they are not exactly the same, but very similar and their difference shows a clear alternating 
  // pattern
  int dummy = 0;

  // maybe we should include an error term for the difference between first and last datapoint - 
  // that may work at least in case of odd N (then the function is odd and the ends are far apart
  // which would get a penalty).

  // todo: compute actual values of the new error function - verify that we are actually at a 
  // minimum by perturbing the first value (which determines the values of all others) and 
  // compute the error for a small upward and downward perturbation from the found minimum
}


// uses the recursive formula found on page 1170 in 
// Mathematical Methods for Physics and Engineering (Riley, Hobson, Bence), 3rd Edition
template<class T>
void rsBinomialDistribution2(T* P, int n, T p)
{
  P[0] = pow(p, n);
  T r  = p/(1-p);
  for(int k = 0; k < n-1; k++)
    P[k+1] = r * P[k] * (n-k) / (k+1);  // Eq. 30.95
}
// would perhaps be nice to use with a rational number type for T


void binomialDistribution()
{
  int    n = 20;                        // number of coin tosses
  double p = 0.5;                       // probability that the result of a single toss is "heads"
  vector<double> P(n+1), Q(n+1);        // probability of seeing P[k] heads in n tosses
  rsBinomialDistribution( &P[0], n, p); // compute probabilities
  rsBinomialDistribution2(&Q[0], n, p); // dito, but uses another algorithm/formula

  GNUPlotter plt;
  plt.addDataArrays(n+1, &P[0]);
  plt.addDataArrays(n+1, &Q[0]);
  plt.plot();
}

// check, if these are already in RSLib and can be deleted here (i think so, but not sure about the
// one using Newton iteration - but maybe this is not needed anyway):
//void sineAmplitudeAndPhaseViaNewton(double y0, double y1, double w, double *a, double *p)
//{
//  double x  = 0.0;
//  double dx = 1.0;
//  double f, fp;
//  double sx, cx, sxw, cxw;
//
//  while( fabs(dx) > fabs(x*EPS) )
//  {
//    rsSinCos(x,   &sx,  &cx);
//    rsSinCos(x+w, &sxw, &cxw);
//    f  = y1*sx - y0*sxw;       // f(x)  = y1*sin(x) - y0*sin(x+w)
//    fp = y1*cx - y0*cxw;       // f'(x) = y1*cos(x) - y0*cos(x+w)
//    dx = -f/fp;
//    x += dx;
//  }
//
//  while( x < 0.0 )
//    x += 2*PI;
//
//  if( fabs(y1) > fabs(y0) )
//    *a = y1 / sxw;
//  else
//    *a = y0 / sx;
//
//  if( *a < 0.0 )
//  {
//    *a = -(*a);
//    x -= PI;
//  }
//
//  *p = x;
//}
//
//void rsSineAmplitudeAndPhase(double y0, double y1, double w, double *a, double *p)
//{
//  double s, c;
//  rsSinCos(w, &s, &c);
//  *p = atan2(y0*s, y1-y0*c);
//  s  = sin(*p);
//  if( fabs(s) > EPS )
//    *a = y0 / s;
//  else
//    *a = y1 / sin(*p + w);
//}
//
//double rsSineFrequency(double y0, double y1, double y2)
//{
//  rsAssert( fabs(y1) > EPS * (fabs(y0)+fabs(y2)), "y1 (numerically) zero is not allowed");
//  return acos(0.5*(y0+y2)/y1); 
//}

void sineParameters()
{
  // From a know amplitude A, normalized radian frequency w and initial phase w, we compute 3 
  // successive sample values of the sequence y[n] = A * sin(w*n + p) and try to retrieve the 
  // sinusoids parameters A, w, p from these sample values. We have a system of 3 nonlinear 
  // equations:
  // y[n]   := y0 = A * sin(p)
  // y[n+1] := y1 = A * sin(p + w) 
  // y[n+2] := y2 = A * sin(p + 2*w) 

  // sine parameters:
  double A = 1.2;
  double p = 1.3;
  double w = 1.7;

  // compute 3 successive values of the sinusoid:
  double y0 = A * sin(p);
  double y1 = A * sin(p+w);
  double y2 = A * sin(p+2*w);

  // retrieve parameters from values:
  double A2, p2, w2;  // retrieved parameters

  // There's a recursion for the sine y[n] = a1*y[n-1] - y[n-2] where a1 = 2*cos(w) and the states
  // y[n-1], y[n-2] are initialized as y[n-1] = A * sin(p - w), y[n-2] = A * sin(p - 2*w) which in 
  // our notation here translates to y2 = a1*y1 - y0. This leads to a1 = (y0+y2)/y1 and 
  // w = acos(a1/2):

  /*
  double a1;          // recursion coefficient
  a1 = (y0+y2)/y1;    // a1 = 2*cos(w)
  w2 = acos(0.5*a1);  // seems to work
  */

  w2 = rsSineFrequency(y0, y1, y2);
  rsSineAmplitudeAndPhase(y0, y1, w2, &A2, &p2);
  int dummy = 0;

  /*
  // We also have:
  // y1 = a1 * A * sin(p - w)         =: y[n-1]
  // y0 =      A * sin(p - 2*w)       =: y[n-2]
  // we may solve the 2nd equation for A = y0/sin(p-2*w) and put it back into the 1st which gives:
  // y1 = a1 * y0/sin(p-2*w) * sin(p-w)
  // y1 * sin(p-2*w) = a1 * y0 * sin(p-w)
  // we may use Newton iteration to solve f(x) = 0 = y1 * sin(x-2*w) - a1 * y0 * sin(x-w)
  // this is all still wrong:
  double x  = 0.0;  // test: later init with 0
  double dx = 1.0;
  double f, fp;
  while( fabs(dx) > EPS )
  {
    f  = y1 * sin(x-2*w) - a1*y0 * sin(x-w);
    fp = y1 * cos(x-2*w) - a1*y0 * cos(x-w);  // can be optimized
    dx = -f/fp;
    x += dx;
    int dummy = 0;
  }
  //p2 = x;
  p2 = x - 2*w;
  while( p2 < 0.0 )
    p2 += 2*PI;
  A2 = y2/(sin(p2-2*w));
  */
}

// Computes the value of a bandlimited unit step (with unit cutoff frequency) at the given 
// time-instant t.
double rsBandLimitedStep(double t)
{
  return rsSineIntegral(PI*t)/PI + 0.5;  // maybe *ONE_OVER_PI instead of /PI
}
// Creates an approximation to a bandlimited step in which the wiggles are confined into some 
// finite time window. The window is centered around t = 0 and extends halfLength into both 
// directions. The shape of the window is a cos^2.
double rsWindowedBandLimitedStep(double t, double halfLength)
{
  double tmp = 0.0; 
  if( fabs(t) < halfLength )
  {
    tmp  = cos(0.5*PI*t/halfLength);
    tmp *= tmp;
  }
  return (1-tmp)*rsStep(t) + tmp*rsBandLimitedStep(t);
}
void bandLimitedStep()
{
  double tMin = -20.0;
  double tMax = +20.0;
  double T    =  15.0;  // half of window length
  static const int N = 5000;
  double t[N];
  double y[N];
  RAPT::rsArrayTools::fillWithRangeLinear(t, N, tMin, tMax);
  for(int n = 0; n < N; n++)
    y[n] = rsWindowedBandLimitedStep(t[n], T);
  plotData(N, t, y);
}


// this commented code can be removed - it already integrated in streamlined form into RSLib
/*
template<class T>
T rsInterpolateCubicHermite(T x1, T x2, T x3, T x4, T y1, T y2, T y3, T y4, T x)
{
  // normalize coordinates:
  T s = T(1)/(x3-x2);
  x1 = s*(x1-x2);
  x3 = s*(x3-x2);  // == 1
  x4 = s*(x4-x2);
  x  = s*(x -x2);
  x2 = s*(x2-x2);  // == 0

  // compute distances between successive normalized x-values:
  T d1, d2, d3;
  d1 = x2-x1;
  d2 = x3-x2;     // == 1 -> remove
  d3 = x4-x3;

  // compute desired slopes at x2 and x3 as weighted averages of the difference quotients left and
  // right to the respective point (x2 or x3):
  T s2, s3;
  s2 = (d2*(y2-y1)/d1 + d1*(y3-y2)/d2) / (d1+d2);
  s3 = (d3*(y3-y2)/d2 + d2*(y4-y3)/d3) / (d2+d3);

  // fit cubic and evaluate:
  T a[4];
  fitCubicWithDerivativeFixedX(y2, y3, s2, s3, &a[3], &a[2], &a[1], &a[0]);
  return evaluatePolynomialAt(x, a, 3);
}
*/
/*
template<class T>
T rsInterpolateCubicHermite(T x1, T x2, T x3, T x4, T y1, T y2, T y3, T y4, T x)
{
  // normalize coordinates:
  T s = T(1)/(x3-x2);
  x1 = s*(x1-x2);
  x4 = s*(x4-x2);

  // compute distances between successive normalized x-values:
  T d1, d3;
  d1 = -x1;
  d3 = x4-1;

  // compute desired slopes at x2 and x3 as weighted averages of the difference quotients left and
  // right to the respective point (x2 or x3):
  T s2, s3;
  s2 = ((y2-y1)/d1 + (y3-y2)*d1) / (d1+1);
  s3 = ((y3-y2)*d3 + (y4-y3)/d3) / (d3+1);

  // fit cubic and evaluate:
  T a[4];
  fitCubicWithDerivativeFixedX(y2, y3, s2, s3, &a[3], &a[2], &a[1], &a[0]);
  return evaluatePolynomialAt(s*(x-x2), a, 3);
  //return evaluatePolynomialAt(x, a, 3);
}
*/
/*
template<class T>
T rsInterpolateCubicHermite(T x1, T x2, T x3, T x4, T y1, T y2, T y3, T y4, T x)
{
  // compute distances between successive normalized x-values:
  T s = T(1)/(x3-x2);
  T d1, d3;
  d1 = s*(x2-x1);
  d3 = s*(x4-x2)-T(1);

  // compute desired slopes at x2 and x3 as weighted averages of the difference quotients left and
  // right to the respective point (x2 or x3):
  T s2, s3;
  s2 = ((y2-y1)/d1 + (y3-y2)*d1) / (d1+T(1));
  s3 = ((y3-y2)*d3 + (y4-y3)/d3) / (d3+T(1));

  // fit cubic and evaluate:
  T a[4];
  fitCubicWithDerivativeFixedX(y2, y3, s2, s3, &a[3], &a[2], &a[1], &a[0]);
  return evaluatePolynomialAt(s*(x-x2), a, 3);
}
*/

void naturalCubicSpline()
{
  // Tests the natural cubic spline interpolation using the test function 1/(1+x^2).

  std::vector<double> x = { -5, -3, -1, 0, 1, 2, 3, 4, 5 };
  int N = (int)x.size();
  std::vector<double> y(N);
  for(int i = 0; i < N; i++)
    y[i] = 1 / (1 + x[i]*x[i]); // f(x) = 1 / (1 + x^2)

  int Ni = 501;
  std::vector<double> xi(Ni), yi(Ni);
  RAPT::rsArrayTools::fillWithRangeLinear(&xi[0], Ni, -6.0, 6.0);

  //RAPT::rsArrayTools::fillWithZeros(&yi[0], Ni);  // preliminary
  RAPT::rsNaturalCubicSpline(&x[0], &y[0], N, &xi[0], &yi[0], Ni);

  GNUPlotter plt;
  plt.addDataArrays(Ni, &xi[0], &yi[0]);
  plt.addGraph("index 0 using 1:2 with lines lw 2 lc rgb \"#808080\" notitle");
  plt.addDataArrays(N, &x[0], &y[0]);
  plt.addGraph("index 1 using 1:2 with points pt 7 ps 1.2 lc rgb \"#000000\" notitle");
  plt.plot();
}

void naturalCubicSpline2()
{
  // Tests the natural cubic spline interpolation using random values with random x-distances 
  // between them. 

  int N = 30; // number of input data points
  std::vector<double> x(N), y(N);
  double dxMin =  1.0;
  double dxMax =  4.0;
  double yMin  = -1.0;
  double yMax  = +1.0;
  double xScl  =  1.0;   // scale factor for x-values


  createRandomDataXY(&x[0], &y[0], N, xScl*dxMin, xScl*dxMax, yMin, yMax);

  int Ni = 1001;
  std::vector<double> xi(Ni), yi(Ni);
  RAPT::rsArrayTools::fillWithRangeLinear(&xi[0], Ni, x[0], x[N-1]);

  RAPT::rsNaturalCubicSpline(&x[0], &y[0], N, &xi[0], &yi[0], Ni, 1.0); // Natural
  //RAPT::rsInterpolateSpline(&x[0], &y[0], N, &xi[0], &yi[0], Ni, 1);   // Hermite

  // maybe factor out (it's the same as above):
  GNUPlotter plt;
  plt.addDataArrays(Ni, &xi[0], &yi[0]);
  plt.addGraph("index 0 using 1:2 with lines lw 2 lc rgb \"#808080\" notitle");
  plt.addDataArrays(N, &x[0], &y[0]);
  plt.addGraph("index 1 using 1:2 with points pt 7 ps 1.2 lc rgb \"#000000\" notitle");
  plt.plot();

  // ToDo:
  // -maybe plot a cubic Hermite spline interpolant for comparison
  //  -Hermite is only 1st order smooth and prescribes *values* for the derivatives at the 
  //   junctions (which we derive from a numerical derivative) whereas the natural spline just
  //   demands that those values must match but doesn't care what the values should be - that 
  //   "gives back" a degree of freedom per datapoint which can be used to match the 2nd derivative
  //   in the same way
  //  -Hermite splines are local: changing the value of a single datapoint affects only the 
  //   segments that go directly into and out of it and its immediate neighbours
}

void cubicInterpolationNonEquidistant() // turn into unit-test
{
  bool testResult = true;

  double x1, x2, x3, x4;     // abscissa (x-axis) values
  double y1, y2, y3, y4;     // ordinate (y-axis) values
  x1 = 2.5; 
  x2 = 3.0; 
  x3 = 3.3; 
  x4 = 4.0;
  y1 = 5.0;
  y2 = 4.0;
  y3 = 7.0;
  y4 = 9.0;

  // compute distances between successive x-values:
  double d1, d2, d3;
  d1 = x2-x1;
  d2 = x3-x2;
  d3 = x4-x3;

  // compute desired slopes at x2 and x3 as weighted averages of the difference quotients left and
  // right to the respective point (x2 or x3):
  double s2, s3;
  s2 = (d2*(y2-y1)/d1 + d1*(y3-y2)/d2) / (d1+d2);
  s3 = (d3*(y3-y2)/d2 + d2*(y4-y3)/d3) / (d2+d3);

  // find polynomial coefficients:
  double a[4];
  fitCubicWithDerivative(x2, x3, y2, y3, s2, s3, &a[3], &a[2], &a[1], &a[0]);

  // interpolate:
  double x = 3.2; // abscissa value to interpolate at
  double y;       // interpolated value at x, computed via a- or b-coeffs
  y = RAPT::rsPolynomial<double>::evaluate(x,  a, 3);


  // do the same in normalized coordinates, such that x2=0, x3=1 - this is done by computing primed
  // x-values xp = (x-x2)/(x3-x2) = (x-x2)/d2
  double x1p, x2p, x3p, x4p, xp;
  xp  = (x -x2)/d2;
  x1p = (x1-x2)/d2;
  x2p = (x2-x2)/d2;  // == 0.0
  x3p = (x3-x2)/d2;  // == 1.0
  x4p = (x4-x2)/d2; 

  double d1p, d2p, d3p;
  d1p = x2p-x1p;
  d2p = x3p-x2p;    // == 1.0
  d3p = x4p-x3p;

  double s2p, s3p;
  s2p = (d2p*(y2-y1)/d1p + d1p*(y3-y2)/d2p) / (d1p+d2p);
  s3p = (d3p*(y3-y2)/d2p + d2p*(y4-y3)/d3p) / (d2p+d3p);

  double b[4];
  fitCubicWithDerivativeFixedX(y2, y3, s2p, s3p, &b[3], &b[2], &b[1], &b[0]);

  double yp = RAPT::rsPolynomial<double>::evaluate(xp, b, 3);

  double error = fabs(y-yp);
  testResult &= fabs(y-yp) < 1.e-10;

  // now, use the function from RAPT:
  double yq = rsInterpolateCubicHermite(x1, x2, x3, x4, y1, y2, y3, y4, x);
  error = fabs(y-yq);
  testResult &= fabs(y-yq) < 1.e-10;
}

void splineInterpolationNonEquidistant()
{
  static const int N = 7;    // number of non-interpolated values
  static const int M = 5000; // number of interpolated values
  double xiMin =  0.45;      // leftmost abscissa value for inter/extrapolated values
  double xiMax =  3.25;      // rightmost abscissa value for inter/extrapolated values
  double xi[M];              // x-values for interpolated data
  double yiSp0[M];           // y-values for spline-interpolated data, smoothness = 0 (linear)
  double yiSp1[M];           // y-values for spline-interpolated data, smoothness = 1 (cubic)
  double yiSp2[M];           // y-values for spline-interpolated data, smoothness = 2 (quintic)
  double yiSp3[M];           // y-values for spline-interpolated data, smoothness = 3 (septic)

  // create input data:
  double xn[N] = {0.52343, 0.93478, 1.12673, 1.73764, 2.13456, 2.46568, 3.23456};
  double yn[N] = {1.0,     2.0,     4.0,     3.0,     7.0,     2.0,     1.0};

  // just a test
  //rsInterpolateSpline(xn, yn, N, xi, yiSp3, M, 10);

  // interpolate and plot interpolated data:
  RAPT::rsArrayTools::fillWithRangeLinear(xi, M, xiMin, xiMax);
  rsInterpolateSpline(xn, yn, N, xi, yiSp0, M, 0);
  rsInterpolateSpline(xn, yn, N, xi, yiSp1, M, 1);
  rsInterpolateSpline(xn, yn, N, xi, yiSp2, M, 2);
  rsInterpolateSpline(xn, yn, N, xi, yiSp3, M, 3);
  plotData(M, xi, yiSp0, yiSp1, yiSp2, yiSp3);

  // Observations:
  // Higher order splines are smoother at the data-points where successive splines join (this is
  // by construction and desirable), but they tend to oscillate more between the data points which 
  // might be undesirable.

  // Ideas: try to use a transformed input variable z = x^c (or z = sign(x) + |x|^c)
}

void rationalInterpolation()
{
  // Consider the problem of finding a rational function:
  //
  //         p(x)     a0 + a1*x + a2*x^2 + ... + aM*x^M
  // r(x) = ------ = -----------------------------------
  //         q(x)     b0 + b1*x + b2*x^2 + ... + bN*x^N
  //
  // that passes through the points (x0 = 0, y0), (x1 = 1, y1) and also satisfies constraints on a 
  // number of derivatives at these points. So we have the data:
  // y(0), y'(0), y''(0), ..., yK(0) and y(1), y'(1), y''(1), ..., yL(1)
  // where yK(0) is the K-th derivative at x=0 and yL(1) is the L-th derivative at x=1. We may also
  // use the notation y0(0), y1(0) instead of y(0), y'(0) such that the index indicates the number
  // of primes. This leads to a generalization of the Hermite interpolation which is the special
  // case of N=0 with b0=1. It should lead to a linear system of equations in which the matrix and
  // right hand side are to be constructed from coefficient arrays of the derivative polynomials 
  // which have to be constructed using the quotient rule. 
  // Maybe, in the application of the quotient rule, the convolutions that have to be performed 
  // there can be truncated after each step at max(M,N) or max(M+1,N+1) or something to avoid 
  // exponential growth of the sequence lengths with the order of the maximum derivative.
  // Maybe it's also possible to skip over the prescription of some derivatives - for example, we
  // may want to prescribe r(1), r''(1) but not r'(1). A special case is when all derivatives
  // are given at x=0, such that L=0 and they are all equal to zero. In this case only the 
  // denominator polynomial of each derivative has to vanish, leading to a possibly simpler 
  // algorithm that uses a Taylor expansion of r(x) around x=0. The successive derivatives may
  // be found by using rsDeConvolve (i think), see the derivation of the Butterworth responses
  // in Paarmann, page 114. The long division in Eq. 3.3 is actually a deconvolution of the 
  // sequences y = {1}, h = {b0, b1, ..., bN}, i.e. y = x * h is the convolution of an unknown 
  // (infinite) sequence x and the polynomial coefficients in the denominator, whcih we may 
  // truncate appropriately.

  // Ahh - damn - no: the system of equations becomes nonlinear. Even in the simple case:
  //
  //         p(x)     a0 + a1*x + a2*x^2
  // r(x) = ------ = --------------------
  //         q(x)     1  + b1*x
  //
  // we get
  //
  //           (a1+2*a2*x)*(1+b1*x) - b1*(a0+a1*x+a2*x^2)
  // r'(x)  = --------------------------------------------
  //                             (1+b1)^2
  //
  // imposing our constraints: r(0) = y0, r'(0) = y0', r(1) = y1, r'(1) = y1'
  // leads to
  // 1: y0  = a0
  // 2: y0' = (a1-b1*a0) / (1+b1)^2
  // 3: y1  = (a0+a1+a2) / (1+b1)
  // 4: y1' = ((a1+2*a2)*(1+b1)-b1*(a0+a1+a2)) / (1+b1)^2
  // solving 2 for a1 and 3 for a2:
  // a1 = y0'(1+b1)^2+b1*a0, a2=y1(1+b1)-a0-a1
  // which can be plugged into 4, leading to a cubic equation for b1. after b1 has been found,
  // a1 and a2 can be found by backsubstitution




  // ....

  int dummy = 0;
}

// move to RAPT:
void rsHermiteCoeffsS1R1(double *y0, double *y1, double *a)
{
  // smoothness = 1, matching integrals = 1
  double k0, k1, k2;
  a[0] = y0[0];
  a[1] = y0[1];
  k0   = y1[0]-a[0]-a[1];
  k1   = y1[1]-a[1];
  k2   = (y1[0]-y0[0]-a[1])/2;
  a[2] = -(-60*k2-3*k1+24*k0)/2;
  a[3] = -60*k2-4*k1+28*k0;
  a[4] = -(-60*k2-5*k1+30*k0)/2;
}
void rsHermiteCoeffsS1R2(double *y0, double *y1, double *a)
{
  // smoothness = 1, matching integrals = 2
  double k0, k1, k2, k3;
  a[0] = y0[0];
  a[1] = y0[1];
  k0   = y1[0]-a[0]-a[1];
  k1   = y1[1]-a[1];
  k2   = (y1[0]-y0[0]-a[1])/2;
  k3   = (2*y0[0]+y1[0]-a[1])/6 - a[0]/2;
  a[2] = 420*k3-180*k2-2*k1+30*k0;
  a[3] = -1680*k3+780*k2+10*k1-140*k0;
  a[4] = 2100*k3-1020*k2-15*k1+195*k0;
  a[5] = -840*k3+420*k2+7*k1-84*k0;
}
void splineInterpolationAreaNormalized()
{
  // We interpolate between (x=0, y=y0), (x=1, y=y1), with derivatives y'(0)=yp0, y'(1)=yp1
  // with a quartic polynomial. In addition to match function values and derivatives at x=0 and 
  // x=1, we impose as 5th condition that the definite integral of the interpolant from 0 to 1 
  // should be equal to the integral that a linear interpolant would give, which is given by:
  // y0 + (y1-y0)/2. The resulting interpolant can be thought of average-preserving between any two
  // data points.
  //
  // We have:
  //
  // p(x)  = a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4
  // p'(x) = a1 + 2*a2*x + 3*a3*x^2 + 4*a4*x^4
  // P(x)  = a0*x + a1*x^2/2 + a2*x^3/3 + a3*x^4/4 + *a4*x^5/5  P(x) is antiderivative of p(x)
  //
  // applying the conditions gives the linear system of equations:
  //
  // y0  = p(0)   = a0
  // yp0 = p'(0)  = a1
  // y1  = p(1)   = a0 + a1 + a2 + a3 + a4
  // yp1 = p'(1)  = a1 + 2*a2 + 3*a3 + 4*a4
  // y0+(y1-y0)/2 = a0 + a1/2 + a2/3 + a3/4 + a4/5
  //
  // a0, a1 are immediately known as a0 = y0, a1 = yp0 which leaves us with the system:
  //
  // y1-a0-a1     = k0 = a2 + a3 + a4
  // yp1-a1       = k1 = 2*a2 + 3*a3 + 4*a4
  // (y1-y0-a1)/2 = k2 = a2/3 + a3/4 + a4/5
  //
  // which gives:
  // 
  // a2 = -(-60*k2-3*k1+24*k0)/2
  // a3 = -60*k2-4*k1+28*k0
  // a4 = -(-60*k2-5*k1+30*k0)/2
  //
  // To generalize the approach, we also fit higher order integrals to the integral of the linear 
  // interpolant. The linear interpolant is given by:
  //
  // yl(x) = b0 + b1*x, where b0 = y0, b1 = y1-y0
  //
  // and its k-th order definite integral from 0 to 1 is given by:
  //
  //       (k+1)*b0 + b1      k*y0 + y1
  // Ik = --------------- = ------------
  //           (k+1)!          (k+1)!
  //
  // so: 
  // I1 = (2*b0+b1) / 2! = b0 + b1/2 = y0 + (y1-y0)/2 = (y0+y1)/2
  // I2 = (3*b0+b1) / 3! = (2*y0+y1)/6
  // I3 = (4*b0+b1) / 4! = (3*y0+y1)/24
  //
  // To match additional integrals, we go for higher order polynomials, so we have to modify the 
  // existing equations and also get additional equations. In the case of matching 2 integrals
  // using a 5th order polynomial:
  //
  // y0  = a0                                                            (as before)
  // yp0 = a1                                                            (as before)
  // y1  = a0 + a1 + a2 + a3 + a4 + a5                                   (additional a5 term)
  // yp1 = a1 + 2*a2 + 3*a3 + 4*a4 + 5*a5                                (additional 5*a5 term)
  // I1  = a0 + a1/2 + a2/3 + a3/4 + a4/5 + a5/6                         (additional a5/6 term)
  // I2  = a0/2 + a1/(2*3) + a2/(3*4) + a3/(4*5) + a4/(5*6) + a5/(6*7)   (new equation)
  // 
  // In general, the new equations have the form:
  //
  //       k*y0 + y1               
  // Ik = ------------ = sum_{n=0}^N ( a[n] / rsProduct(n+1, n+k) )
  //        (k+1)!
  //
  // but those terms of the right hand side which are already known (a[0] = y(0), a[1] = y'(0),...)
  // have to be brought over to the left hand side when establishing the equation system. Maybe, as
  // a generalization, we could pass the integral values Ik as parameters along with the derivative
  // values and only optionally let them be calculated like this (which would correspond to 
  // evaluating the integral numerically from the function values using the trapezoidal rule). So
  // the interpolation routine would take an array of x-values, an array of y-values, a number M
  // of arrays of derivatives and a number L of arrays of integrals. The polynomial order would be
  // N = 2*M+L+1


  // set y-values and derivatives:
  double y0[2] = { 5, -1 };  // p(0)=5, p'(0)=-1
  double y1[2] = { 9, -4 };  // p(1)=9, p'(1)=-4 

  // compute linear interpolation coefficients for the line y(x) = b0 + b1*x:
  double b0, b1;
  b0 = y0[0];
  b1 = y1[0] - y0[0];

  // compute integrals of linear interpolant:
  double I1, I2;
  I1 = (2*b0+b1) / rsFactorial(2);
  I2 = (3*b0+b1) / rsFactorial(3);

  // compute cubic polynomial coefficients without regularization:
  double ac[4];
  getHermiteCoeffs1(y0, y1, ac);

  // compute quartic polynomial coefficients with regularization:
  double ar1[5];
  rsHermiteCoeffsS1R1(y0, y1, ar1);

  // compute quintic polynomial coefficients with 2nd order regularization:
  double ar2[6];
  rsHermiteCoeffsS1R2(y0, y1, ar2);

  // create linear, cubic, quartic and quintic interpolant:
  static const int N = 1000;
  double x[N], yl[N], yc[N], yr1[N], yr2[N];
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, 0.0, 1.0);
  for(int n = 0; n < N; n++)
  {
    yl[n]  = rsInterpolateLinear(0.0, 1.0, y0[0], y1[0], x[n]);
    yc[n]  = RAPT::rsPolynomial<double>::evaluate(x[n], ac,  3);
    yr1[n] = RAPT::rsPolynomial<double>::evaluate(x[n], ar1, 4);
    yr2[n] = RAPT::rsPolynomial<double>::evaluate(x[n], ar2, 5);
  }

  // create the running sums of the interpolants (which are approximations to integrals times N):
  double sl1[N], sl2[N], sc1[N], sc2[N], sr11[N], sr12[N], sr21[N], sr22[N];
  RAPT::rsArrayTools::cumulativeSum(yl,  sl1,  N);     // 1st order cumulative sum of linear interpolant
  RAPT::rsArrayTools::cumulativeSum(yl,  sl2,  N, 2);  // 2nd order cumulative sum of linear interpolant
  RAPT::rsArrayTools::cumulativeSum(yc,  sc1,  N);     // 1st order cumulative sum of cubic interpolant
  RAPT::rsArrayTools::cumulativeSum(yc,  sc2,  N, 2);  // 2nd order cumulative sum of cubic interpolant
  RAPT::rsArrayTools::cumulativeSum(yr1, sr11, N);     // 1st order cumulative sum of quartic interpolant
  RAPT::rsArrayTools::cumulativeSum(yr1, sr12, N, 2);  // 2nd order cumulative sum of quartic interpolant
  RAPT::rsArrayTools::cumulativeSum(yr2, sr21, N);     // 1st order cumulative sum of quintic interpolant
  RAPT::rsArrayTools::cumulativeSum(yr2, sr22, N, 2);  // 2nd order cumulative sum of quintic interpolant

  // plot:
  //plotData(N, x, yl,  yc,  yr1);   // linear, cubic and quartic interpolant
  plotData(N, x, yl,  yc,  yr1,  yr2);   // linear, cubic, quartic and quintic interpolant
  //plotData(N, x, sl1, sc1, sr11, sr21);  // 1st order running sums
  //plotData(N, x, sl2, sc2, sr12, sr22);  // 2nd order running sums

  // Observations:
  // Matching the 1st order integral seems to make sense because it reduces the tendency to
  // overshoot around the data points. When matching the 2nd as well, the overshoot is further 
  // reduced but the interpolating function shows additional wiggles in between the data points 
  // which would seem undesirable in most situations. So, matching even more integrals doesn't 
  // seem to make much sense.

  // todo: maybe the approach can easily be generalized to spline interpolation with arbitrary 
  //       order, i.e. an arbitrary number of matching derivatives
  //
  // todo: maybe instead of matching the integral, match the arithmetic average of the interpolated
  //       values to that of the linear interpolant - such that we don't do area normalization in
  //       the continuous domain but rather in the resampled discrete domain
  //
  // todo: maybe use other functions with a 5th parameter like:
  //
  //         a0 + a1*x + a2*x^2 + a3*x^3
  // f(x) = -----------------------------
  //         1  + bn*x^2
  //
  // where n is 1,2 or 3 ...but the derivative and antiderivative might become more complicated
  // and nonlinear - don't know, if that's practical
}

// ToDo: consider monotonic spline interpolation
// We must make sure that the derivative of our interpolating polynomial is nonnegative on our 
// interval 0..1. A polynomial p(x) is nonnegative on an interval a..b, iff it can be written as:
// p(x) = s(x) + (x-a)*(b-x)*t(x) for even degree
// p(x) = (x-a)*s(x) + (b-x)*t(x) for odd degree
// where s(x),t(x) are sums of squares. A polynomial is a sum of squares, iff it can be written as:
// p(x) = sum_k q_k^2(x)
// See here for more details:
// http://math.stackexchange.com/questions/60610/polynomial-fitting-where-polynomial-must-be-monotonically-increasing
// http://stellar.mit.edu/S/course/6/sp10/6.256/courseMaterial/topics/topic2/lectureNotes/lecture-10/lecture-10.pdf
// Presumably, to interpolate with 1st derivative matched, we will need an additional degree of 
// freedom (the monotonicity is kind of an additional constraint?), so we would have to use a 4th
// order polynomial.




// todo: implement trapezoidal rule a la numerical recipies and driver routine

template<class T>
class rsTrapezoidalStage
{

public:

  // maybe rename updateAndGetResult
  T update(const std::function<T(T)>& f, T a, T b, int n)
  {
    if( n == 1 )
      return s = T(0.5) * (b-a) * (f(a) + f(b));
    else 
    {
      //int i = 1;  // numIterations
      //for(int j = 1; j < n-1; j++)  
      //  i <<= 1;                   // i *= 2
      //// todo: use an integer-power function: i = 2^(n-2)

      int i = rsPowInt(2, n-2);
      T t   = T(i);
      T dx  = (b-a) / t;         // spacing of points to be added
      T x   = a + T(0.5) * dx;   // current evaluation point
      T sum = T(0);              // rename to s (rename member s to "result" or something
      for(int j = 1; j <= i; j++) {
        x   += dx;
        sum += f(x);
      }
      s = 0.5 * (s+(b-a)*sum/t); // update with refined value
      return s;
    }
  }

protected:

  T s = T(0); 

};

template<class T>
T integrateTrapezoidal(const std::function<T(T)>& f, T a, T b)
{
  T tol = T(100) * std::numeric_limits<T>::epsilon();
  int minNumStages = 5;
  int maxNumStages = 25;
  rsTrapezoidalStage<T> stage;
  T s, sOld = 0.0;

  for(int j = 1; j <= maxNumStages; j++) {
    s = stage.update(f, a, b, j);
    if(j > minNumStages) {
      if(rsAbs(s-sOld) < tol*rsAbs(sOld) || (s == T(0) && sOld == T(0)) ) // why the 2nd criterion?
        return s;
      sOld = s;
    }
  }

  //rsError("integrateTrapezoidal failed to converge");
  return s;
}
// not recommended for production code - too inefficient (takes many stages)

template<class T>
T integrateSimpson(const std::function<T(T)>& f, T a, T b)
{
  T tol = T(100) * std::numeric_limits<T>::epsilon();
  int minNumStages = 5;
  int maxNumStages = 25;
  rsTrapezoidalStage<T> stage;
  T s, st, ost = T(0), os = T(0);

  for(int j = 1; j <= maxNumStages; j++) {
    st = stage.update(f, a, b, j); // trapezoidal sum
    s  = (T(4)*st - ost)/T(3);     // simpson sum
    if(j > minNumStages) {
      if(rsAbs(s-os) < tol*rsAbs(os) || (s == T(0) && os == T(0)) ) // why the 2nd criterion?
        return s;
      os  = s;
      ost = st;
    }
  }

  //rsError("integrateTrapezoidal failed to converge");
  return s;
}
// this doesn't semm to be much of an improvement over integrateTrapezoidal
// try romberg integration and/or Ooura's routines: http://www.kurims.kyoto-u.ac.jp/~ooura/

template<class F>
double forwardDifference(const F& f, double x, double h)
{
  return (f(x+h) - f(x)) / h;
}

void numericDifferentiation()
{
  // When using numerical differentiation formulas, there are two sources of error: the error 
  // coming from the approximation itself and the roundoff error due to finite precision 
  // arithmetic. The first error decreases with decreasing approximation stepsize h, whereas the 
  // second increases with decreasing stepsize (with smaller h we tend subtract very close numbers,
  // in the formulas leading to cancellation of significant digits), so there will be a sweet spot
  // for h, where the accuracy of a given formula is best. Where that sweet spot is depends on the
  // x-scale of the function: for sin(2*x) we will have to choose a h of half the size compared to
  // sin(x) (..verify). In the range, where the finite precision effects are not yet severe, we 
  // observe an error increase with increasing h that follows a power rule, i.e. the error 
  // increases with h^n for some exponent n. In a 2nd order accurate formula, n=2, for example. 
  // Also, when h becomes too large, the power law ceases to hold because we are then not really 
  // in the neighbourhood of our evaulation point anymore, so the whole idea of approximating 
  // tangents by secants and stuff like that becomes invalid.

  // In this experiment, we investigate the behavior the error of the numerically computed 
  // derivative as function of the stepsize h.
  // todo: 
  // -investigate the optimal choice of h as function of the x-scaling of the function - i 
  //  expect inverse proportionality
  // -check various formulas - currently, we only look at the basic central difference formula
  //  -> we should see different slopes of the lines in the log-log plot of the error for the
  //     different formulas


  // h should scale like 2^k where k goes from , say -20 to -3
  static const int kMin = -25;
  static const int kMax =   5;
  static const int numK = kMax - kMin + 1;

  double x0 = 1.0;  
  // we compute the error at position x0 ..maybe compute at various positions and average? 
  // x0 = 1 seems to be a good choice for a "general" point of the sine function (no symmetries 
  // around or near it)


  //double xMin = 0.0;
  //double xMax = 2*PI;

  using NumDif = rsNumericDifferentiator<double>;

  double h[numK];
  double err[numK];
  double k[numK];
  double logErr[numK];


  // function to differentiate:
  auto f = [=](double x)
  {
    return sin(x);
  };

  for(int i = 0; i < numK; i++)
  {
    //double h = pow(2, kMin + i);
    k[i] = kMin + i;
    h[i] = pow(2.0, k[i]);
    double t = cos(x0);                         // true derivative
    double a = NumDif::derivative(f, x0, h[i]); // actually computed numerical derivative

    //a = forwardDifference(f, x0, h[i]);  // test

    err[i] = fabs(t-a);                         // error
    logErr[i] = rsLog2(err[i]);
    int dummy = 0;
  }

  double a, b;
  rsStatistics::linearRegression(numK, k, logErr, a, b);
  // a is equal to 2 - i think, that's what "2nd order in h" means, todo: try it with a forward
  // difference - this should be 1st order in h - a should come out as 1. then do it for the higher
  // order formulas and also for the higher derivatives...

  //rsPlotArraysXY(numK, h, err);
  rsPlotArraysXY(numK, k, logErr); // a nice straight line, as expected

  // todo: estimate the exponent by fitting a line to logErr as function of k



  // figure out the range of k for which this rule holds - it will be bounded below by roundoff 
  // error becoming a problem and above by the fact that the sine is periodic - when h is equal to
  // the period (or even just half of it) the derivative estimate will be zero even though the 
  // correct value is 1 - but what, if the example function is exp instead of sin? will the rule
  // then continue to hold with no upper limit (or well, overflow being the limit?), generally, the
  // approximation is valid only in the neighborhood of x0, so we will be bounded above by the
  // approximation becoming invalid - what counts as neighbourhood depends on the x-scale of the 
  // function - for sin(2*x), a "neighborhood" will only be half as wide as for sin(x)

  // Observations:
  // -choosing the test point as x0 = 0, we get a nice straight line and the computed order is 2
  //  even for the forward difference (may this have to do with the fact that sin(x) is symmetric
  //  around 0?)
  // -choosing a more general point such as x0 = 1, we see, that the accuracy is best at k = -17, 
  //  the genral, straight-line rule holds from k = -16 upward
  // -with x0 = 1, f = sin(x), plotting a range from k = -25..+5 shows, where the linear rule stops
  //  to hold (namely, for k = -16...+1)...that it goes up so high is actually surprising - k=1, 
  //  means h = 2, which seems way too large to count as neighborhood - what's going on?
}

void numericIntegration()
{
  std::function<double(double)> f;
  double a, b, I;

  // integral_0^pi sin(x) dx = 2:
  f = [=](double x) { return sin(x); };
  a = 0;
  b = PI;
  I = integrateTrapezoidal(f, a, b);
  I = integrateSimpson(    f, a, b);

  // integral_0^1 4/(pi*(1+x^2)) = 1:
  // wolfram alpha: integral 4/(pi*(1+x^2)) from 0 to 1
  f = [=](double x) { return 4/(PI*(1+x*x)); };
  a = 0;
  b = 1;
  I = integrateTrapezoidal(f, a, b);
  I = integrateSimpson(    f, a, b);


  int dummy = 0;
}

void nonUniformArrayDiffAndInt()
{
  // Test numerical differentiation and integration routines. We sample a sinewave at 
  // nonequidistant sample points and find the numeric derivative and  integral at these sample
  // points and compare them to the true derivative/integral values.

  static const int N = 100;   // number of sample points
  double p = 1.0;             // start-phase
  double w = 2.0;             // radian frequency 
  double xMax = 10.0;         // maximum x-axis value
  double x[N];                // x-axis values
  double y[N], yd[N], ydn[N]; // y, y' and numeric y'
  double       yi[N], yin[N]; // true and numeric integral

  // create x-axis:
  RAPT::rsArrayTools::fillWithRandomValues(x, N, 0.1, 1.5, 0);
  RAPT::rsArrayTools::cumulativeSum(x, x, N);
  double scaler = xMax/x[N-1];
  RAPT::rsArrayTools::scale(x, N, scaler);

  // compute sine and its derivative and integral analytically at the samples:
  int n;
  for(n = 0; n < N; n++)
  {
    y[n]  =        sin(w*x[n] + p);
    yd[n] =      w*cos(w*x[n] + p);
    yi[n] = -(1/w)*cos(w*x[n] + p);
  }

  // compute the numeric derivative and integral:
  rsNumericDifferentiator<double>::derivative(x, y, ydn, N, true);
  rsNumericIntegral(  x, y, yin, N, yi[0]);

  // plot function, true derivative and numeric derivative:
  //plotData(N, x, y, yd, ydn);
  plotData(N, x, y, yd, ydn, yi, yin);
}


/*
template<class T>
void integrateTrapezoidal(T* x, T* y, int N, T y0 = T(0))
{
  T x1 = x[0];
  y[0] = y0;
  for(int n = 1; n < N; n++)
  {
    y[n] = y[n-1]
  }
}
// test, if this works also in place - it should! ..maybe make a separate function to work in place
// which takes only one input array
*/

template<class T>
void integrateTrapezoidal(const T* x, T* y, int N, T y0 = T(0))
{
  rsAssert(y != x, "This function does not work in place");
  y[0] = y0;
  for(int n = 1; n < N; n++)
    y[n] = y[n-1] + 0.5*(x[n-1] + x[n]);
}

template<class T>
void differentiateTrapezoidal(const T* x, T* y, int N, T y0 = T(0))
{
  rsAssert(y != x, "This function does not work in place");
  y[0] = y0;
  for(int n = 1; n < N; n++)
    y[n] = 2*(x[n] - x[n-1]) - y[n-1];
}

void uniformArrayDiffAndInt()
{
  // todo: 
  // -create array of random numbers x[n]
  // -compute cumulative sum
  // -compute difference of cumulative sum - should give back original array
  // -compute trapezoidal integration of x[n]
  // -compute an inverse of trapezoidal integration ("trapezoidal differentiation"?)
  // ...maybe this should become a unit test...

  // trapezoidal integration: 
  //   y[0] = y0;                              // y0 is integration constant
  //   y[n] = y[n-1] + 0.5*(x[n-1] + x[n]);
  // trapezoidal differentiation:
  //   x[0] = ?
  //   x[n] = 2*(y[n] - y[n-1]) - x[n-1]

  static const int N = 100;

  using AT = RAPT::rsArrayTools;

  double x[N];  // input signal
  double y[N];  // integrator output signal
  double z[N];  // differentiator output signal
  double e[N];  // error

  // test cumulative sum and backward difference:
  AT::fillWithRandomValues(x, N, -1.0, +1.0, 0);
  AT::cumulativeSum(x, y, N);
  AT::copy(y, z, N);
  AT::difference(z, N);  // todo: make a function that works not in place
  AT::subtract(x, z, e, N);
  //rsPlotArrays(N, x, y, z, e);  // ok - looks good

  // test trapezoidal integration and differentiation:
  double y0 = 3;
  integrateTrapezoidal(x, y, N, y0);
  differentiateTrapezoidal(y, z, N, x[0]);
  AT::subtract(x, z, e, N);
  //rsPlotArrays(N, x, y, z, e);

  // plot cumulative sum vs trapezoidal integral:
  AT::cumulativeSum(x, y, N);
  integrateTrapezoidal(x, z, N);
  rsPlotArrays(N, x, y, z);

  // Observations:
  // -the tarpezoidal integral looks smoother and delayed compared to the cumulative sum - i think,
  //  the trapezoidal integrator is basically a cumulative sum of data that was passed through a 
  //  2-point moving average which explains both effects
  //  todo: try with an impulse input to figure out the amount of the delay - is it half a sample?
  // -AT::difference does indeed undo AT::cumulativeSum
  // -differentiateTrapezoidal undoes integrateTrapezoidal only, if we pass the correct value x[0]
  //  for the initial state - otherwise, we see an error that oscillates at the Nyquist freq 
  //  -passing x[0]-a for some a > 0 gives an error alternating between a and -a starting at +a and 
  //   passing x[0]+a gives alternating error starting at -a
  //  -whatever we pass as y0 to the integrator doesn't matter
  // -i think, the greater smoothness of the trapezoidal integral implies that the trapezoidal
  //  differentiator introduces more jaggies than the simple differencer - so it may actually be
  //  undesirable to use trapezoidal differentiation for obtaining instantaneous frequencies from
  //  instantaneous phases, for example in rsSingleSineModeler - i think, if we use it there, it
  //  should be optional


  int dummy = 0;
}

// ToDo: implement a numerical differentiation algorithm that is the inverse operation of 
// trapezoidal integration. Let x,y be the arrays of abcissa and ordinate value to be integrated
// and z the array of the integral values. Then, trapezoidal integration does the following:
//   z[0] = c; z[i] = z[i-1] + (x[i] - x[i-1]) * (y[i] + y[i-1]) / 2
// where c is the integration constant. The formula can be solved for y[i]:
//   y[i] = 2 * (z[i] - z[i-1]) / (x[i] - x[i-1]) - y[i-1]
// or for y[i-1]:
//   y[i-1] = 2 * (z[i] - z[i-1]) / (x[i] - x[i-1]) - y[i]
// so we have a nice recursion formula that we may either start at y[0] or y[N-1] - but how do we 
// get the boundary values? Let's unroll the recursion, starting at the end, at y[N-1]:
//   y[N-1] = 2*(z[N-1]-z[N-2])/(x[N-1]-x[N-2]) - y[N-2]   
// replace y[N-2] by its expression:
//   y[N-1] = 2*(z[N-1]-z[N-2])/(x[N-1]-x[N-2]) - 2*(z[N-2]-z[N-3])/(x[N-2]-x[N-3]) + y[N-3]
// replace y[N-3] by its expression:
//   y[N-1] = 2*(z[N-1]-z[N-2])/(x[N-1]-x[N-2]) - 2*(z[N-2]-z[N-3])/(x[N-2]-x[N-3]) + 2*(z[N-3]-z[N-4])/(x[N-3]-x[N-4]) - y[N-4]
// and so on (note the alternating signs) - so after unrolling everything, this leads to:
//   y[N-1] =   sum_{i odd}     2 * (z[N-i] - z[N-i-1]) / (x[N-i] - x[N-i-1])
//            - sum_{i>0, even} 2 * (z[N-i] - z[N-i-1]) / (x[N-i] - x[N-i-1])
// so - if i didn't make any mistakes (which is unlikely), this is the formula, we may use to 
// compute the boundary value y[N-1] to start the recursion from the end of the array. Once we have
// that, we can apply the recursion formula to compute all other values. ...maybe define arrays to
// hold the differences of the x and z arrays, like a[i] = x[i] - x[i-1], b[i] = z[i] - z[i-1]
// ...but what about a[0], b[0]?

// or - see it like a filter - the bilinear/trapezoidal integrator (using filter x,y notation) is:
//   y[n] = y[n-1] + 0.5 * (x[n] + x[n-1])
// and just solve this for x[n]:
//   x[n] = 2*(y[n] - y[n-1]) - x[n-1]   // x,y (in/out) have swapped roles
// of course, this works only, if the sample-spacing is 1 - but the formula can be generalized:
//   y[n] = y[n-1] + 0.5*dt * (x[n] + x[n-1])
//   x[n] = (y[n] - y[n-1]) / (0.5*dt) - x[n-1]
// where dt is the temporal spacing between sample x[n-1] and x[n], i.e. dt = t[n] - t[n-1], if we have a 
// t-array with time-stamps. The simplemost running-sum 
// integrator and the corresponding differentiatior would be:
//   y[n] = y[n-1] + dt * x[n]
//   x[n] = (y[n] - y[n-1]) / dt
// Maybe implement them and use them in a non-uniform version of the SVF. But what about this 
// "embedded integrator gain" thing? Also, maybe implement and expeeriment with PID controllers
// https://en.wikipedia.org/wiki/PID_controller ..maybe such a thing can be used for dynamics 
// processors?



void shiftPolynomial()
{
  static const int order = 6;
  double p[order+1]  = {2,1,-5,7,-3,2,-2}; // p(x) = -2x^6+2x^5-3x^4+7x^3-5x^2+1x^1+2x^0
  double x0          = 2.0;                  // shift value  

  double xMin = -1.0;
  double xMax = +1.0;
  static const int N = 1000;
  double x[N];
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, xMin, xMax);
  double y[N], ys[N], yst[N];   // y, stretched version of y, target for stretched version
  int n;
  for(n = 0; n < N; n++)
  {
    y[n]   = RAPT::rsPolynomial<double>::evaluate(x[n],    p, order);
    yst[n] = RAPT::rsPolynomial<double>::evaluate(x[n]-x0, p, order);
  }

  // establish coeffs of q(x) = p(x-x0):
  double q[order+1];
  RAPT::rsPolynomial<double>::coeffsForShiftedArgument(p, q, order, x0);

  // test - use composePolynomials using p(r(x)) where r(x) = x-x0, with x0 = 2.0
  //double r[2] = {-2.0, 1};
  //composePolynomials(r, 1, p, 6, q); 
    // ...yes - gives the same result - todo: find out, if polyCoeffsForShiftedArgument is actually
    // more efficient - otherwise, we may not need it

  // evaluate q from x-array:
  for(n = 0; n < N; n++)
    ys[n] = RAPT::rsPolynomial<double>::evaluate(x[n], q, order);

  plotData(N, x, y, yst, ys);
}

// void stretchPolynomial()
// {
// }






/*
void monotonicPolynomials()
{
  // Consider the general polynomial of order N:
  //
  // p(x)  = a0 + a1*x + a2*x^2 + a3*x^3 + ... + aN*x^N   with the derivative
  // p'(x) = a1 + 2*a2*x + 3*a3*x^2 + ... + N*aN*x^(N-1)
  //
  // we want the polynomial to be monotonically increasing, such that p'(x) >= 0 for all x. The 
  // transition between monotonic and nonmonotonic polynomials occurs at the critical case where
  // p'(x) = 0 for N-1 values of x (really? - verify this)


  // 
  // Try to construct a 4th order polynomial which has zero derivative at x=-1, x=0, x=1
  //

  static const int N = 6; // order of the polynomial
  double a[N+1];          // polynomial coefficients
  double x[N-1];          // x values where p'(x) = 0
  double yp[N];           // y': all zero vector of the derivatives at the selected x points
                          // and unity as last value
  double ap[N];           // a': coefficients of the derivative polynomial p'(x)

  double shift = 3.0;     // integration constant - determines function value at x=0
  double scale = 2.0;     // overall scaling - determines ap[N-1], the highest power coefficient
                          // in the derivative polynomial

  //rsFillWithRangeLinear(x, N-1, double(-(N-1)/2), double((N-1)-1-(N-1)/2)); // simplify

  rsFillWithRangeLinear(x, N-1, 0.0, double(N-2)); // simplify
  rsFillWithValue(yp, N-1, 1.0);
  yp[N-1] = scale;


  double **A = rsVandermondeMatrix(x, N);

  // modify last line...
  rsFillWithZeros(A[N-1], N-1);
  A[N-1][N-1] = 1;

  rsSolveLinearSystem(A, ap, yp, N);
  rsDeAllocateSquareArray2D(A, N);  

  polyIntegral(ap, a, N-1, shift);
    
  // a and ap are equal to zero - this is the trivial case that indeed satifies the constraints but
  // we need to ensure that ap[N-1] != 0
  // aahh - no - we need to establish another constraint: the leading coefficient of the derivative
  // polynomial should be nonzero (maybe, it makes most sense to fix it to unity): apN = 1.
  // This leaves us only with N-1 degrees of freedom to be used to set derivatives to zero.
  // If we want to prescribe N derivatives, the only polynomials that can satisfy this constraint
  // will have the same derivate everywhere, i.e. linear functions.

  // plot:
  static const int Np = 1000;
  double xPlt[Np], yPlt[Np];
  double xMin = rsMinValue(x, N-1);
  double xMax = rsMaxValue(x, N-1);
  rsFillWithRangeLinear(xPlt, Np, xMin, xMax);
  for(int n = 0; n < Np; n++)
    yPlt[n] = RSLib::evaluatePolynomialAt(xPlt[n], a, N);
  plotData(Np, xPlt, yPlt);

  // Observations: When choosing the points where the derivative should vanish at equidistant
  // x values, the polynomial is not monotonic.

  int dummy = 0;
}
*/


// Computes coefficients for an N-th order polynomial, whose N-1 values of the 1st derivative at
// x[0], ..., x[N-2] are given by yp[0],...,yp[N-2].
// a is of length N+1, x is of length N-1, yp is of length N where the additional last value of yp
// yp[N-1] which has no corresponding x-value is an overall scale factor for the whole polynomial. 
// It's a bit inconvenient but it's the way the data needs to be arranged for solving the linear 
// system.
void rsPolyWithDerivativeValues(double *a, int N, double *x, double *yp, double shift = 0)
{
  double **A = RAPT::rsPolynomial<double>::vandermondeMatrix(x, N);
  RAPT::rsArrayTools::fillWithZeros(A[N-1], N-1);
  A[N-1][N-1] = 1;
  RAPT::rsLinearAlgebra::rsSolveLinearSystem(A, a, yp, N);
  RAPT::rsArrayTools::deAllocateSquareArray2D(A, N);  
  RAPT::rsPolynomial<double>::integral(a, a, N-1, shift);
}
void monotonicPolynomials()
{
  // Consider the general polynomial of order N:
  //
  // p(x)  = a0 + a1*x + a2*x^2 + a3*x^3 + ... + aN*x^N   with the derivative
  // p'(x) = a1 + 2*a2*x + 3*a3*x^2 + ... + N*aN*x^(N-1)
  //
  // we want the polynomial to be monotonically increasing, such that p'(x) >= 0 for all x. The 
  // transition between monotonic and nonmonotonic polynomials occurs at the critical case where
  // p'(x) = 0 for N-1 values of x (really? - verify this)

  // With the function rsPolyWithDerivativeValues, we can construct polynomials the have arbitrary
  // values for the 1st derivative at N-1 selected points. Now, we have to find the points at which 
  // we should set the derivative zero in order to obtain a monotonic polynomial.

  static const int N = 3; // order of the polynomial
  double a[N+1];          // polynomial coefficients
  double x[N-1];          // x values where p'(x) = 0
  double yp[N];           // y': all zero vector of the derivatives at the selected x points
                          // and unity as last value
  //double ap[N];           // a': coefficients of the derivative polynomial p'(x)

  double shift = 0.0;     // integration constant - determines function value at x=0
  double scale = 1.0;     // overall scaling - determines ap[N-1], the highest power coefficient
                          // in the derivative polynomial

  //rsFillWithRangeLinear(x, N-1, double(-(N-1)/2), double((N-1)-1-(N-1)/2)); // simplify
  RAPT::rsArrayTools::fillWithRangeLinear(x, N-1, 0.0, double(N-2));

  // try to hand select values for x
  /*
  x[0] = 0;
  x[1] = x[0] + 1;
  x[2] = x[1] + 2;
  x[3] = x[2] + 3;
  x[4] = x[3] + 4;
  */

  x[0] = 0;
  x[1] = 2;


  RAPT::rsArrayTools::fillWithValue(yp, N-1, 0.0);
  yp[N-1] = scale;
  rsPolyWithDerivativeValues(a, N, x, yp, shift);

    
  // a and ap are equal to zero - this is the trivial case that indeed satifies the constraints but
  // we need to ensure that ap[N-1] != 0
  // aahh - no - we need to establish another constraint: the leading coefficient of the derivative
  // polynomial should be nonzero (maybe, it makes most sense to fix it to unity): apN = 1.
  // This leaves us only with N-1 degrees of freedom to be used to set derivatives to zero.
  // If we want to prescribe N derivatives, the only polynomials that can satisfy this constraint
  // will have the same derivate everywhere, i.e. linear functions.

  // plot:
  static const int Np = 1000;
  double xPlt[Np], yPlt[Np];
  double xMin = RAPT::rsArrayTools::minValue(x, N-1);
  double xMax = RAPT::rsArrayTools::maxValue(x, N-1);

  xMax = 3;
  RAPT::rsArrayTools::fillWithRangeLinear(xPlt, Np, xMin, xMax);
  for(int n = 0; n < Np; n++)
    yPlt[n] = RAPT::rsPolynomial<double>::evaluate(xPlt[n], a, N);
  plotData(Np, xPlt, yPlt);

  // Observations: When choosing the points where the derivative should vanish at equidistant
  // x values, the polynomial is not monotonic.

  // Links for study:
  // http://math.stackexchange.com/questions/60610/polynomial-fitting-where-polynomial-must-be-monotonically-increasing
  // https://en.wikipedia.org/wiki/Sturm's_theorem (recommended in the discussion above) 
  // http://math.ucsb.edu/~padraic/mathcamp_2013/root_find_alg/Mathcamp_2013_Root-Finding_Algorithms_Day_2.pdf (more about Sturm's theorem)
  // http://stellar.mit.edu/S/course/6/sp10/6.256/courseMaterial/topics/topic2/lectureNotes/lecture-10/lecture-10.pdf
  // https://www.physicsforums.com/threads/monotonic-polynomial.90551/
  // http://haralick.org/conferences/73102292.pdf
  // http://cran.r-project.org/web/packages/MonoPoly/MonoPoly.pdf
  // http://www.cse.ucla.edu/products/overheads/IMPS2013/Falk20130722.pdf

  int dummy = 0;
}

void parametricBell()
{
  static const int N = 1000;
  double center =  10.0;  // center of the bell
  double width  =  4.0;   // width
  double flat   =  0.5;   // relative length of flat top zone (between 0 and 1)

  // create and set up parametric bell object:
  rsParametricBellFunction<double> bell;
  bell.setCenter(center);
  bell.setWidth(width);
  bell.setFlatTopWidth(flat);

  // create x-axis and allocate y arrays:
  double xMin =  center - 1.2 * width/2;
  double xMax =  center + 1.2 * width/2;
  double x[N];
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, xMin, xMax);
  double yl[N], yc[N], yq[N], yh[N], yb[N]; // linear, cubic, quintic, heptic, bump
  int n;

  // create the family of curves (we look at different shapes for the prototype bell):
  bell.setPrototypeBell(&rsPositiveBellFunctions<double>::linear);
  for(n = 0; n < N; n++)
    yl[n] = bell.getValue(x[n]);
  bell.setPrototypeBell(&rsPositiveBellFunctions<double>::cubic);
  for(n = 0; n < N; n++)
    yc[n] = bell.getValue(x[n]);
  bell.setPrototypeBell(&rsPositiveBellFunctions<double>::quintic);
  for(n = 0; n < N; n++)
    yq[n] = bell.getValue(x[n]);
  bell.setPrototypeBell(&rsPositiveBellFunctions<double>::heptic);
  for(n = 0; n < N; n++)
    yh[n] = bell.getValue(x[n]);
  bell.setPrototypeBell(&rsPositiveBellFunctions<double>::bump);
  for(n = 0; n < N; n++)
    yb[n] = bell.getValue(x[n]);

  GNUPlotter plt;
  plt.addDataArrays(N, x, yl, yc, yq, yh, yb);
  plt.plot();
}

void partialFractionExpansion()
{
  // We recover the partial fraction expansion of the strictly proper rational function:
  // f(x) = 2/(x+1) - 3/(x-2) + 5/(x-5) 
  //      = (4x^2 - 7x + 25) / ((x+1)*(x-2)*(x-5))
  //      = (4x^2 - 7x + 25) / (x^3 - 6x^2 + 3x + 10)
  // via various algorithms. Example taken from "Hˆhere Mathematik f¸r Naturwissenschaftler und
  // Ingenieure", 2.Auflage (B‰rwolff), page 147

  static const int N = 2;      // numerator order
  static const int M = 3;      // denominator order
  double p[N+1] = {25,-7,4};   // p(x) =       4x^2 - 7x + 25
  double q[M+1] = {10,3,-6,1}; // q(x) = x^3 - 6x^2 + 3x + 10
  double r[M]   = {-1,2,5};    // roots of q(x)


  typedef RAPT::rsPolynomial<double> Poly;
  typedef std::complex<double> Complex;
  typedef RAPT::rsArrayTools Array;

  // probably, we have to divide the numerator coefficients by the leading coefficient of the 
  // denominator - in our example, it's unity, so it doesn't matter:
  Array::scale(p, N+1, 1.0/q[M]);

  // First, we use the approach via establishing and solving the linear system of equations that
  // results from equating coefficients:
  double A[M][M];
  double tmp[M+1];
  double dummy;
  int i, j;
  for(i = 0; i < M; i++) { // establish coefficient matrix
    Array::copy(q, tmp, M+1);
    Poly::divideByMonomialInPlace(tmp, M, r[i], &dummy);
      // todo: use a function that does not do it "in-place" - avoids copying and is probably 
      // simpler. perhaps, here, we have to do that division in a loop from 1 up to the 
      // multiplicity of the root r[i] - but where would the result go in the coefficient matrix?
    for(j = 0; j < M; j++)
      A[j][i] = tmp[j];
  }
  double x[3];
  RAPT::rsLinearAlgebra::rsSolveLinearSystem3x3(A, x, p);   // solve the linear system
  // x == {2,-3,5}, so: f(x) = 2/(x+1) - 3/(x-2) + 5/(x-5)
  // ok - this works - but it's inefficient for larger systems

  // now try a different approach - evaluating the numerator p(x) at the pole r_i and dividing by 
  // polynomial resulting from dividing out r_i from the denominator also gives the residuum a_i
  // corresponding to r_i, see page 148 (but i think, this only works when all poles are distinct, 
  // which is the case here:
  Complex rc[M];  // roots as complex numbers (as required by the evaluation routine)
  Complex xc[3];
  Array::convert(r, rc, M);
  Complex num, den;
  for(i = 0; i < M; i++) {
    num = Poly::evaluate(r[i], p, N);
    den = Poly::evaluateFromRootsOneLeftOut(Complex(r[i]), rc, M, i);
    double y, yd; // ...verify, if q'(r_i) is the same as den...
    Poly::evaluateWithDerivative(r[i], q, M, &y, &yd); // ...yes - it is!
    xc[i] = num/den;
  }
  // xc == {2,-3,5} - this works, too and is probably more desirable for production code - but it 
  // requires all poles to be distinct - can the approach be modified to also work for multiple 
  // poles?
  // I think, this is this method:
  // https://en.wikipedia.org/wiki/Heaviside_cover-up_method
  // here:
  // http://math.mit.edu/~jorloff/suppnotes/suppnotes03/h.pdf
  // it says, this method can only be used to find the coefficient for the highest power of a 
  // bunch of terms like a_i1/(x-r_i) + a_i2/(x-r_i)^2 + a_i3/(x-r_i)^3 + ..., so this means, we 
  // could find only a_i3 with this method if terms up to .../(x-r_i)^3 are present?

  dummy = 0;

  // ...and if this is not possible, then what about the approach taken here:
  // https://en.wikipedia.org/wiki/Partial_fraction_decomposition#Residue_method
  // is this the same? that would mean that the derivative of q at the pole r_i is the same as q 
  // with the linear factor (x-r_i) divided out ...right? verify this....it seems to be
  // but it seems, the residue method is useless for n-th poles since it includes taking 
  // the n-th derivative of a rational function - and taking derivatives of rational functions 
  // makes their order blow up exponentially (due to the g^2 term in the quotient rule) - each 
  // derivative doubles the order of the denominator


  // this here is also a good summary over all the relevant methods:
  // https://stanford.edu/~boyd/ee102/rational.pdf
  // which is part of this course:
  // https://stanford.edu/~boyd/ee102/

  // trying to find a way to figure out coeffs for multiple poles via the example
  // f(x) = 2/(x-1) + 3/(x-1)^2 + 4/(x-3)  with double-pole at x=1 and simple-pole at x=3
  // f.simplify_rational()
  // (6*x^2 - 13*x + 1)/(x^3 - 5*x^2 + 7*x - 3)

  // to recover the resiude 2, we need to consider:
  // f(x) = (6*x^2-13*x+1)/((x-1)*(x-3)) - 3/(x-1)
  // plot(f, 0, 2, gridlines=True)
  // it indeed goes through 2 at x=1, even though, it's actually not really defined there due to
  // the singularities - but apparently, they cancel each other...somehow subtractively inf-inf or
  // something, see here:
  // https://math.oregonstate.edu/home/programs/undergrad/CalculusQuestStudyGuides/SandS/lHopital/inf_minus_inf.html
  // -> yep - that works, returns the correct result and seems to lead to a feasible algorithm
  // -> todo: work out the algorithm in full generality


  // f(x) = 2/(x-1) + 3/(x-1)^2 + 1/(x-1)^3 + 1/(x-2) with triple-pole at x=1 and simple-root at x=2
  // f.simplify_rational()
  // (3*x^3 - 8*x^2 + 5*x - 1)/(x^4 - 5*x^3 + 9*x^2 - 7*x + 2)

  // papers on algorithms for partial fraction decomposition:
  // http://www.eecs.harvard.edu/~htk/publication/1977-siam-journal-on-computing-kung-tong.pdf
  // http://www-users.math.umn.edu/~am/Math5583-CompAnal/PartialFractions.pdf
  // https://arxiv.org/pdf/math/0408189.pdf
  // http://people.math.sfu.ca/~kya17/teaching/math343/3-343.pdf
}


// sanity checks for inputs to partialFractions
bool checkPartialFractionInputs(
  const RAPT::rsPolynomial<double>& B,
  const RAPT::rsPolynomial<double>& A,
  const std::vector<double>& p,
  const std::vector<int>& m)
{
  bool r = true;
  r &= B.getDegree() < A.getDegree();    // rational function must be strictly proper
  r &= A.getDegree() == RAPT::rsSum(m);  // pole multiplicities must sum to denominator degree
  r &= A.isMonic();                      // denominator must be monic
  //r &= A.hasRootsAt(p, m, tol); // should evaluate also derivatives for multiple roots
  // have function isRootWithMultiplicity(r, m, tol) ...or something
  return r;
}

std::vector<double> partialFractions(
  const RAPT::rsPolynomial<double>& B, 
  const RAPT::rsPolynomial<double>& A, 
  const std::vector<double>& p, 
  const std::vector<int>& m)  
{
  rsAssert(checkPartialFractionInputs(B, A, p, m)); // sanity checks for inputs
  typedef RAPT::rsPolynomial<double> Poly;
  int N = A.getDegree();                     // degree of denominator = number of residues
  int k, j0 = 0;                             // flat array index into r and base-index
  std::vector<double> r(N);                  // array of residues
  Poly Li, Bij, Aij, Lij, Cij;               // the involved polynomials
  double num, den;                           // numerator and denominator, evaluated at the pole
  for(int i = 0; i < (int)p.size(); i++) {   // loop over distinct poles
    Li  = std::vector<double>{ -p[i], 1.0 }; // linear factor (x - p[i])
    Bij = B;                                 // needed for computation of numerator
    Aij = A / (Li^m[i]);                     // parentheses needed, ^ has lower precedence than /
    den = Aij(p[i]);                         // we need to evaluate this only once
    for(int j = m[i]; j >= 1; j--) {         // loop over the exponents for pole i
      k    = j0+j-1;                         // index into r-array
      Lij  = Li^(m[i]-j);                    // Shlemiel the painter strikes again
      Cij  = Bij / Lij;                      // cancel common factor with denominator
      num  = Cij(p[i]);                      // evaluate numerator
      r[k] = num / den;                      // denominator stays the same inside the j-loop
      Bij  = Bij - Aij * r[k];               // establish B-polynomial for next iteration
      Aij  = Aij * Li;                       // establish A-polynomial for next iteration
    }
    j0 += m[i];                              // increment base-index
  }
  return r;
}
// production code should repeated use divideByMonomialInPlace to create the Lij polynomial to 
// avoid the Shlemiel 

// 2nd version, using l'Hospital instead of cancellation - the algorithm has the same overall 
// struture, so only the differences to the version above are commented:
std::vector<double> partialFractions2( 
  const RAPT::rsPolynomial<double>& B, 
  const RAPT::rsPolynomial<double>& A, 
  const std::vector<double>& p, 
  const std::vector<int>& m)  
{
  rsAssert(checkPartialFractionInputs(B, A, p, m));
  typedef RAPT::rsPolynomial<double> Poly;
  int N = A.getDegree();
  int k, j0 = 0;
  std::vector<double> r(N);
  Poly Li, Bij, Aij, Lij;                        // no Cij needed anymore
  double num, den;
  for(int i = 0; i < (int)p.size(); i++) { 
    Li  = std::vector<double>{ -p[i], 1.0 };
    Bij = B; 
    Aij = A / (Li^m[i]);
    for(int j = m[i]; j >= 1; j--) {
      k    = j0+j-1;
      Lij  = Li^(m[i]-j);
      num  = Bij.derivativeAt(p[i], m[i]-j);
      den  = Aij.derivativeAt(p[i], m[i]-j);     // den must now be evaluated for each j
      r[k] = num / den;
      Bij  = Bij - Aij * r[k];
      Aij  = Aij * Li;
    }
    j0 += m[i]; 
  }
  return r;
}

void partialFractionExpansion2()
{
  // Here, i implement a prototype of my new algorithm for a partial fraction expansion with 
  // multiple poles (see PartialFractions.txt in the Notes folder)

  // shorthands for convenience:
  typedef RAPT::rsPolynomial<double> Poly;
  typedef RAPT::rsRationalFunction<double> RF;
  typedef std::vector<double> VecD;
  typedef std::vector<int> VecI;

  // f(x) = 2/(x-1) + 3/(x-1)^2 + 1/(x-1)^3 + 1/(x-2)
  //      = (3x^3 - 8x^2 + 5x - 1) / ((x-1)^3 * (x-2))
  //      = (3x^3 - 8x^2 + 5x - 1) / (x^4 - 5x^3 + 9x^2 - 7x + 2)
  Poly B(VecD{ -1,  5, -8,  3   });         // numerator
  Poly A(VecD{  2, -7,  9, -5, 1});         // denominator (must be monic)
  VecD p =   {  1,  2 };                    // distinct poles
  VecI m =   {  3,  1 };                    // multiplicities of the poles
  VecD r = partialFractions( B, A, p, m);   // find residues - should be: 2,3,1,1
  r      = partialFractions2(B, A, p, m);   // 2nd version of algo

  // f(x) = 2/(x-1) + 3/(x-1)^2 + 4/(x-3)
  //      = (6x^2 - 13x + 1) / ((x-1)^2 * (x-3))
  //      = (6x^2 - 13x + 1) / (x^3 - 5x^2 + 7x - 3)
  B = VecD({ 1, -13,  6   });
  A = VecD({-3,   7, -5, 1});
  p =      { 1,   3 };
  m =      { 2,   1 };
  r = partialFractions( B, A, p, m); // 2,3,4
  r = partialFractions2(B, A, p, m);

  // f(x) = 1/(x-1)+2/(x-1)^2+3/(x-2)+4/(x-2)^2+5/(x-2)^3+6/(x-2)^4+7/(x-3)+8/(x-3)^2+9/(x-3)^3
  //
  //          11x^8 - 167x^7 + 1098x^6 - 4059x^5 + 9150x^4 - 12713x^3 + 10403x^2 - 4399x + 660
  // = ----------------------------------------------------------------------------------------
  //    x^9 - 19x^8 + 158x^7 -  754x^6 + 2273x^5 - 4483x^4 +  5776x^3 -  4680x^2 + 2160x - 432
  B = VecD({ 660, -4399, 10403, -12713,  9150, -4059, 1098, -167,  11   });
  A = VecD({-432,  2160, -4680,   5776, -4483,  2273, -754,  158, -19, 1});
  p = { 1, 2, 3 };
  m = { 2, 4, 3 };
  r = partialFractions( B, A, p, m); // 1,2,3,4,5,6,7,8,9
  r = partialFractions2(B, A, p, m);
  //r = RF::partialFractions(B, A, p, m); // doesn't work because parameters are complex


  // testing new algo to evaluate derivative of polynomial - move to unit tests:
  VecD testPoly = { 1,1,1,1,1,1 }; // 1 + x^2 + x^3 + x^4 + x^5
  double val;
  val = Poly::evaluateDerivative(1.0, &testPoly[0], 5, 0); //   6
  val = Poly::evaluateDerivative(1.0, &testPoly[0], 5, 1); //  15
  val = Poly::evaluateDerivative(1.0, &testPoly[0], 5);    //  15
  val = Poly::evaluateDerivative(1.0, &testPoly[0], 5, 2); //  40
  val = Poly::evaluateDerivative(1.0, &testPoly[0], 5, 3); //  90
  val = Poly::evaluateDerivative(1.0, &testPoly[0], 5, 4); // 144
  val = Poly::evaluateDerivative(1.0, &testPoly[0], 5, 5); // 120
  val = Poly::evaluateDerivative(1.0, &testPoly[0], 5, 6); //   0
  val = Poly::evaluateDerivative(1.0, &testPoly[0], 5, 7); //   0
  val = Poly::evaluateDerivative(1.0, &testPoly[0], 5, 8); //   0
  // this algo should go into class rsPolynomial and the .derivativeAt function should be based on
  // it - this is then used inside partialFractions2
  // why does it actually produce 0 for the derivatives of order higher than 5? because the product
  // starts at a negative integer and ends at a positive, so it includes 0 as factor - but we 
  // should treat the case order > degree separately



  int dummy = 0; 

  // todo: 
  // -move partialFractions/2 functions to Prototypes.h/cpp in rs_testing module
  // -turn this experiment into unit-test
  // -analyze scaling - i think, it's still quadratic in the degree of the num and/or denom due to
  //  the polynomial operations in the inner loop (which are themsleves linear in the degree)
  // -implement production versions that operate on raw arrays and avoid unnecesarry memory
  //  re-allocations, use accumulation, avoid shlemiel-the-painter
  // -in production code, we may check, if m[i] == 1 and if it is, do the simpler step like in the
  //  all-poles distinct case - avoids some overhead in the case of a mix of single and multiple 
  //  poles
  // -make performance and precision tests - also for the algo that uses the linear system approach
  // -maybe production code should be able to deal with non-monic denominators? we could easily do 
  //  it by initializing Aij,Bij to be c*A, c*B where c=1/leadingCoeff(A)

  // see also:
  // Partial Fraction Decomposition by Repeated Synthetic Division
  // https://pdfs.semanticscholar.org/57a6/8c9cee292db58a6c6f9deac76e8ef686085d.pdf
}

/*
void partialFractionExpansion2()
{
  // We try to recover the partial fraction expansion of the strictly proper rational function:
  // f(x) = 3/(x+5) - 4/(x+3) + 2/(x-1) + 5/(x-1)^2 - 3/(x-5)
  //      = P(x)/Q(x) = numerator(x) / denominator(x)
  //      = (-2*x^4-13*x^3+25*x^2-275*x-215)/(x^5+x^4-30*x^3-22*x^2+125*x-75)

  // todo: make a version of that test, where all poles are distinct...

  typedef std::complex<double> Complex;

  static const int Np = 4;  // numerator order
  static const int Nq = 5;  // denominator order
  static const int numPoles         = 5;
  static const int numDistinctPoles = 4;           // only 4, because of them is a double-root
  Complex p[Np+1] = {-215,-275,25,-13,-2};         // numerator
  Complex q[Nq+1] = {-75,125,-22,-30,1,1};         // denominator
  Complex poles[numDistinctPoles] = {-5,-3,+1,+5}; // positions of distinct poles
  int muls[numDistinctPoles] = { 1, 1, 2, 1};      // pole multiplicities

  // compute the results:
  Complex pfeCoeffs[numPoles];
  RAPT::rsRationalFunction<double>::partialFractionExpansion(
    p, Np, q, Nq, poles, muls, numDistinctPoles, pfeCoeffs);

  // todo: turn into unit test - done - this function is now obslete - delete it
  int dummy = 0;  
}
*/

void partialFractionExpansion3()
{
  // We test the RAPT::rsRationalFunction::partialFractionExpansion with a function that contains
  // a polynomial part and partial fractions with different multiplicities:
  // f(x) = 3*x^2-2*x+7 + 3/(x+5) - 4/(x+3) + 2/(x-1) + 5/(x-1)^2 - 3/(x-5)
  //      = 3*x^2-2*x+7 + (-2*x^4-13*x^3+25*x^2-275*x-215)/(x^5+x^4-30*x^3-22*x^2+125*x-75)
  //      = a(x) + b(x)/q(x) = polynomial(x) + numerator(x) / denominator(x)
  //      =        p(x)/q(x) =                 numerator(x) / denominator(x)
  // where p(x) = q(x)*a(x) + b(x)
  // sage:
  // f(x) = 3*x^2-2*x+7 + 3/(x+5) - 4/(x+3) + 2/(x-1) + 5/(x-1)^2 - 3/(x-5)
  // f.simplify_rational()
  // gives:
  // (3*x^7+x^6-85*x^5-x^4+196*x^3-604*x^2+750*x-740) / (x^5+x^4-30*x^3-22*x^2+125*x-75)

  typedef std::complex<double> Complex;
  typedef RAPT::rsPolynomial<Complex> Poly;

  static const int Na = 2, Nb = 4, Nq = 5, Np = Na+Nq;  // degrees of the various polynomials
  Complex a[Na+1] = { 7,-2,3 };
  Complex b[Nb+1] = {-215,-275,25,-13,-2};
  Complex q[Nq+1] = {-75,125,-22,-30,1,1};              // what, if q isn't monic?
  Complex p[Np+1] = {-740,750,-604,196,-1,-85,1,3};

  static const int numPoles         = 5;
  static const int numDistinctPoles = 4;                // only 4, because is a double-pole
  Complex poles[numDistinctPoles]   = {-5,-3,+1,+5};
  int muls[numDistinctPoles]        = { 1, 1, 2, 1};

  // compute the results:
  Complex pfeCoeffs[numPoles];
  Complex polyCoeffs[Np+1];     // can we get away with allocating only Na+1? ...nope!
  RAPT::rsRationalFunction<double>::partialFractionExpansion(
    p, Np, q, Nq, poles, muls, numDistinctPoles, pfeCoeffs, polyCoeffs);

  // ok - this seems to work - polyCoeffs come out as 7,-2,3 and pfeCoeffs as 3,-4,2,5,-3
  // todo: turn into unit test
  int dummy = 0;  
}
/*
void partialFractionExpansion4()
{
  // We try different algorithms for a partial fraction expansion of a strictly proper rational
  // function with all poles distinct

}
*/



/** For two complex numbers r1, r2 which are supposed to be roots of a polynomial, this function
decides, if these should be considered as distinct roots or a double root. When computing roots of 
polynomials, a double root may appear as 2 distinct roots with slightly different values due to 
numerical errors, so our equality check should include some error tolerance. 
\todo implement an error tolerant equality check - at the moment, it just checks for exact equality
*/
bool rsAreRootsDistinct(rsComplexDbl r1, rsComplexDbl r2)
{
  if( r1 == r2 )  
    return false;
  return true;
  // preliminary - use error-tolerant inequality check later. idea: take d = r2-r1 as the 
  // difference, if abs(d)^2 / (abs(r1)^2 + abs(r2)^2) < k*eps, consider them equal
  // abs(d)^2 = d.re*d.re + d.im*d.im, etc., the division normalizes by the length's of the
  // actual roots to get a relative error. the squaring means that we want the relative distance
  // to be of the order of sqrt(eps), which is reasonable according to "Numerical Recipies"
}

/** Given the parameters a0, a1, p1, p2 of the rational function with already factored denominator:

        a0 + a1*x
f(x) = ------------
       (x-p1)(x-p2)

this function computes the residues r1, r2 such that f(x) can be expressed as a sum of partial 
fractions:

        r1     r2
f(x) = ---- + ----
       x-p1   x-p2   */
void quadraticResiduesFromPoles(rsComplexDbl a0, rsComplexDbl a1, rsComplexDbl p1, rsComplexDbl p2,
                                rsComplexDbl *r1, rsComplexDbl *r2)
{
  if( rsAreRootsDistinct(p1, p2) )
  {
    *r1 = (a0 + a1*p1) / (p1-p2);
    *r2 = a1 - *r1;
  }
  else
  {
    *r1 = a1;
    *r2 = a0 + a1*p1;
  }
}

/** Computes the 2 solutions (a.k.a. "roots") of the quadratic equation: 
x^2 + p*x + q = 0 and stores them in r1 and r2. */
void quadraticRoots(rsComplexDbl p, rsComplexDbl q, rsComplexDbl *r1, rsComplexDbl *r2)
{
  rsComplexDbl t1 = -0.5*p;
  rsComplexDbl t2 = sqrt(0.25*p*p - q);
  *r1 = t1 + t2;
  *r2 = t1 - t2;
}

/** Given the parameters a0, a1, b0, b1 of the rational function:

           a0 + a1*x
f(x) = ----------------
        b0 + b1*x + x^2

this function computes the poles p1, p2 and residues r1, r2 such that f(x) can be expressed as a 
sum of partial fractions:

        r1     r2
f(x) = ---- + -----
       x-p1    x-p2   
*/
void quadraticPartialFractionExpansion(rsComplexDbl a0, rsComplexDbl a1, rsComplexDbl b0, 
   rsComplexDbl b1, rsComplexDbl *r1, rsComplexDbl *p1, rsComplexDbl *r2, rsComplexDbl *p2)
 {
   quadraticRoots(b1, b0, p1, p2);
   quadraticResiduesFromPoles(a0, a1, *p1, *p2, r1, r2);
 }

/** Given the parameters b0, b1, a1, a2 of the rational function:

            b0 + b1/z
H(z) = -------------------
        1 + a1/z + a2/z^2

this function computes the poles p1, p2 and residues r1, r2 such that f(x) can be expressed as a 
sum of partial fractions:

         r1       r2
H(z) = ------ + ------
       1-p1/z   1-p2/z   
       
It's just a convenience function that uses the partial fraction expansion for positive powers 
internally, but the numerator coefficients are swapped and in the denominator, the coefficient for 
z^-2 takes the role of the constant coefficient in the normal powers case. The swap of a- and 
b-variables is just a result of using a different notational convention here to be consistent with 
much of the DSP literature (b-coeffs are in the numerator, a-coeffs in the denominator).

Maybe for the general case with arbitrary order, we must just reverse the polynomial coefficient 
arrays and pass them to a regular partial fraction expansion routine when we are dealing with 
negative powers? -> quite possibly, but check the math
*/
void quadraticPartialFractionExpansionNegativePowers(rsComplexDbl b0, rsComplexDbl b1, 
  rsComplexDbl a1, rsComplexDbl a2, rsComplexDbl *r1, rsComplexDbl *p1, rsComplexDbl *r2, 
  rsComplexDbl *p2)
 {
   quadraticPartialFractionExpansion(b1, b0, a2, a1, r1, p1, r2, p2);
 }

/** Inverse of quadraticPartialFractionExpansion. See comments there. */
void quadraticPartialFractionComposition(rsComplexDbl r1, rsComplexDbl p1, rsComplexDbl r2, 
  rsComplexDbl p2, rsComplexDbl *a0, rsComplexDbl *a1, rsComplexDbl *b0, rsComplexDbl *b1)
 {
   if( rsAreRootsDistinct(p1, p2) )
   {
     *a0 = -(r1*p2 + r2*p1);
     *a1 = r1 + r2;
   }
   else
   {
     *a0 = r2 - r1*p1;
     *a1 = r1;
   }
   *b0 = p1*p2;
   *b1 = -(p1+p2);
}

/** Inverse of quadraticPartialFractionExpansionNegativePowers. See comments there. */
void quadraticPartialFractionCompositionNegativePowers(rsComplexDbl r1, rsComplexDbl p1, 
  rsComplexDbl r2, rsComplexDbl p2, rsComplexDbl *b0, rsComplexDbl *b1, rsComplexDbl *a1, 
  rsComplexDbl *a2)
{
  quadraticPartialFractionComposition(r1, p1, r2, p2, b1, b0, a2, a1);
}

// write function to compose quadratic function in inverse powers
                               
void partialFractionExpansionQuadratic()
{
  // may be moved to unit-tests later and move the functions that are tested in rsRationalFunction

  // We consider the rational function with 1st order numerator and 2nd order (quadratic) 
  // denominator:
  //
  //            a0 + a1*x
  // f(x) =  ----------------
  //         b0 + b1*x + x^2 
  //
  // and want to find its partial fraction expansion. Let the roots of the numerator be denoted
  // as p1, p2. If they are distinct, a partial fraction expansion of the form:
  //
  //          r1     r2
  // f(x) =  ---- + ----
  //         x-p1   x-p2
  //
  // exists. In the case of a pole p with multiplicity 2, the expansion takes the form:
  //
  //        r1      r2
  // f(x) = --- + -------
  //        x-p   (x-p)^2
  //
  // We test functions that convert back and forth between (a0, a1, b0, b1) and (r1, p1, r2, p2).

  bool testResult = true;
  rsComplexDbl i(0, 1);  // imaginary unit


  // two real poles at 5 and 1, with residues 3 and 2 respectively:
  //         3     2       5x - 13
  // f(x) = --- + --- = -------------
  //        x-5   x-1    x^2 - 6x + 5
  rsComplexDbl a0 = -13;
  rsComplexDbl a1 =   5;
  rsComplexDbl b0 =   5;
  rsComplexDbl b1 =  -6;
  rsComplexDbl r1, p1, r2, p2;
  quadraticPartialFractionExpansion(a0, a1, b0, b1, &r1, &p1, &r2, &p2);
  testResult &= r1 == 3.;
  testResult &= p1 == 5.;
  testResult &= r2 == 2.;
  testResult &= p2 == 1.;
  quadraticPartialFractionComposition(r1, p1, r2, p2, &a0, &a1, &b0, &b1);
  testResult &= a0 == -13.;
  testResult &= a1 ==   5.;
  testResult &= b0 ==   5.;
  testResult &= b1 ==  -6.;

  // complex conjugate poles and residues:
  //          7+5i       7-5i        14x - 58
  // f(x) = -------- + -------- = --------------
  //        x-(2+3i)   x-(2-3i)    x^2 - 4x + 13
  a0 = -58;
  a1 =  14;
  b0 =  13;
  b1 =  -4;
  quadraticPartialFractionExpansion(a0, a1, b0, b1, &r1, &p1, &r2, &p2);
  testResult &= r1 == 7.0 + 5.0*i;
  testResult &= r2 == 7.0 - 5.0*i;
  testResult &= p1 == 2.0 + 3.0*i;
  testResult &= p2 == 2.0 - 3.0*i;
  quadraticPartialFractionComposition(r1, p1, r2, p2, &a0, &a1, &b0, &b1);
  testResult &= a0 == -58.;
  testResult &= a1 ==  14.;
  testResult &= b0 ==  13.;
  testResult &= b1 ==  -4.;

  // complex conjugate poles and real residues:
  // ...

  // double pole at 3 with residues 2, 5 for 1st and 2nd order fraction
  //         2       5         2x - 1
  // f(x) = --- + ------- = ------------
  //        x-3   (x-3)^2   x^2 - 6x + 9
  a0 = -1;
  a1 =  2;
  b0 =  9;
  b1 = -6;
  quadraticPartialFractionExpansion(a0, a1, b0, b1, &r1, &p1, &r2, &p2);
  testResult &= r1 == 2.;
  testResult &= p1 == 3.;
  testResult &= r2 == 5.;
  testResult &= p2 == 3.;
  quadraticPartialFractionComposition(r1, p1, r2, p2, &a0, &a1, &b0, &b1);
  testResult &= a0 == -1.;
  testResult &= a1 ==  2.;
  testResult &= b0 ==  9.;
  testResult &= b1 == -6.;

  // test with totally arbitrary complex numbers for all values:
  a0 =  -70.0 +  1.0*i;
  a1 =    5.0 -  2.0*i;
  b0 =   59.0 +  2.0*i;
  b1 =  -10.0 +  2.0*i;
  quadraticPartialFractionExpansion(a0, a1, b0, b1, &r1, &p1, &r2, &p2);
  testResult &= r1 == 3.0 - 5.0*i;
  testResult &= p1 == 6.0 - 7.0*i;
  testResult &= r2 == 2.0 + 3.0*i;
  testResult &= p2 == 4.0 + 5.0*i;
  quadraticPartialFractionComposition(r1, p1, r2, p2, &a0, &a1, &b0, &b1);
  testResult &= a0 == -70.0 +  1.0*i;
  testResult &= a1 ==   5.0 -  2.0*i;
  testResult &= b0 ==  59.0 +  2.0*i;
  testResult &= b1 == -10.0 +  2.0*i;

  // test the function for inverse powers (relevant for digital filters):
  //
  //          3        2       5 - 13/z         
  // H(z) = ------ + ----- = ---------------
  //        1-5/z    1-1/z   1 - 6/z + 5/z^2 
  rsComplexDbl a2;
  b0 =   5;
  b1 = -13;
  a1 =  -6;
  a2 =   5;
  quadraticPartialFractionExpansionNegativePowers(b0, b1, a1, a2, &r1, &p1, &r2, &p2);
  rsComplexDbl z  = 0.5 + 0.2*i;
  rsComplexDbl Hd = (5. - 13./z) / (1. - 6./z + 5./(z*z)); // H(z) directly computed
  rsComplexDbl Hp = 3./(1.-5./z) + 2./(1.-1./z);           // H(z) via partial fractions
  rsComplexDbl d  = Hd-Hp;
  testResult &= abs(d) < 1.e-14;
  quadraticPartialFractionCompositionNegativePowers(r1, p1, r2, p2, &b0, &b1, &a1, &a2);
  testResult &= b0 ==   5.;
  testResult &= b1 == -13.;
  testResult &= a1 ==  -6.;
  testResult &= a2 ==   5.;



  int dummy = 0;
}


/** Definite integral from t=0 to infinity of f(t) = exp(-a*t) * sin(w*t + p). */
double rsDampedSineIntegral(double a, double w, double p)
{
  double s, c;
  rsSinCos(p, &s, &c);
  return (w*c + a*s) / (a*a + w*w);
}

/** Definite integral from t=0 to infinity of f(t) = exp(-a*t) * cos(w*t + p). */
inline double rsDampedCosineIntegral(double a, double w, double p)
{
  return rsDampedSineIntegral(a, w, p + PI/2);  // cos(x) = sin(x + pi/2);
}

/** Definite integral from t=0 to infinity of (f(t))^2 where 
f(t) = A * exp(-a*t) * sin(w*t + p). */
double rsDampedSineTotalEnergy(double A, double a, double w, double p)
{
  double s, c; rsSinCos(2*p, &s, &c);
  double a2   = a*a;
  double w2   = w*w;
  double a2w2 = a2+w2;
  return A*A*(a*w*s-a2*c+a2w2) / (4*a*a2w2);
}

/** Definite integral from t=0 to T of (f(t))^2 where 
f(t) = A * exp(-a*t) * sin(w*t + p). */
double rsDampedSineEarlyEnergy(double A, double a, double w, double p, double T)
{
  double a2   = a*a;
  double w2   = w*w;
  double a2w2 = a2+w2;
  double d    = 4*a*a2w2;
  double s, c; 
  rsSinCos(2*(p+w*T), &s, &c);
  double F_T  = (a*w*s-a2*c+a2w2) / (d*exp(2*a*T)); // F(T), up to scaling
  rsSinCos(2*p, &s, &c);
  double F_0  = (a*w*s-a2*c+a2w2) / d;              // F(0), up to scaling
  return -A*A * (F_T - F_0);
}



/** Definite integral from t=0 to infinity of (f(t))^2 where 
f(t) = A * (exp(-a1*t) - exp(-a2*t)) * sin(w*t+p). */
double rsAttackDecaySineTotalEnergy(double A, double a1, double a2, double w, double p)
{
  double s, c;  rsSinCos(2*p, &s, &c);
  double w2w2 = 2*w; w2w2 *= w2w2;
  double ws2  = 2*w*s;
  double a1a1 = a1+a1;
  double a2a2 = a2+a2;
  double a1a2 = a1+a2;
  double tmp  = 1/a1a1 + 1/a2a2 - 2/a1a2;
  tmp -=   (a1a1*c - ws2) / (a1a1*a1a1 + w2w2);
  tmp -=   (a2a2*c - ws2) / (a2a2*a2a2 + w2w2);
  tmp += 2*(a1a2*c - ws2) / (a1a2*a1a2 + w2w2);
  return 0.5*A*A*tmp;
}

/** Definite integral from t=0 to T of exp(-a*t)*(sin(w*t+p))^2. */
double rsDampedSineSquaredIntegral(double a, double w, double p, double T)
{
  double s, c; rsSinCos(2*(w*T+p), &s, &c);
  double ws4 = 4*w*w;
  return exp(-a*T) * (a*(a*c-a-2*w*s)-ws4) / (2*a*(a*a+ws4));
}

/** Definite integral from t=0 to T of (f(t))^2 where 
f(t) = A * (exp(-a1*t) - exp(-a2*t)) * sin(w*t+p). */
double rsAttackDecaySineEarlyEnergy(double A, double a1, double a2, double w, double p, double T)
{
  double FT, F0;  // F(T), F(0), F is the antiderivative

  // straightforward, but inefficient (lots of redundant calculations):
  FT  =   rsDampedSineSquaredIntegral(a1+a1, w, p, T);
  FT -= 2*rsDampedSineSquaredIntegral(a1+a2, w, p, T);
  FT +=   rsDampedSineSquaredIntegral(a2+a2, w, p, T);
  F0  =   rsDampedSineSquaredIntegral(a1+a1, w, p, 0);
  F0 -= 2*rsDampedSineSquaredIntegral(a1+a2, w, p, 0);
  F0 +=   rsDampedSineSquaredIntegral(a2+a2, w, p, 0); 

  // there are a lot of common subexpressions - optimize...
  //double x    = T;
  //double ea1  = exp(-a1*x);
  //double ea2  = exp(-a2*x);
  //double ea12 = ea1*ea2;     // e^(-(a1+a2)*x)
  // ....

  return A*A*(FT-F0);
}

void dampedSineEnergy()
{
  double A  =  2.0;
  double a1 =  2.3;
  double w  = 21.7;
  double p  =  0.6;

  double E = rsDampedSineTotalEnergy(A, a1, w, p);

  // find a value for E by approximating the integral with a sum:
  double t    = 0.0;
  double dt   = 0.00001;
  double En   = 0.0;
  double tau  = 1/a1;

  double dE;
  double f;
  while( t < 10*tau )
  {
    f   = A * exp(-a1*t) * sin(w*t + p);
    dE  = (f*f) * dt;
    En += dE;
    t  += dt;
  }
   // yes - the smaller dt, the closer En is to E, as expected

  // now, the same procedure for the attack/decay envelope:
  double a2 = 2*a1;  
  E  = rsAttackDecaySineTotalEnergy(A, a1, a2, w, p);
  t  = 0.0;
  En = 0.0; 
  while( t < 10*tau )
  {
    f   = A * (exp(-a1*t) - exp(-a2*t)) * sin(w*t + p);
    dE  = (f*f) * dt;
    En += dE;
    t  += dt;
  }

  static const int N = 1000;
  double tMax = 10*tau;
  double tAxis[N];
  RAPT::rsArrayTools::fillWithRangeLinear(tAxis, N, 0.0, tMax);
  double y[N], ySq[N], yI[N];
  int n;
  for(n = 0; n < N; n++)
  {
    t      = tAxis[n];
    y[n]   = A * exp(-a1*t) * sin(w*t + p);
    ySq[n] = y[n] * y[n];
    yI[n]  = rsDampedSineEarlyEnergy(A, a1, w, p, t);
  }
  double test = yI[N-1];



  for(n = 0; n < N; n++)
  {
    t      = tAxis[n];
    y[n]   = A * (exp(-a1*t)-exp(-a2*t)) * sin(w*t + p);
    ySq[n] = y[n] * y[n];
    yI[n]  = rsAttackDecaySineEarlyEnergy(A, a1, a2, w, p, t);
  }
  test = yI[200]; // 0.065141113152824xxx

  //plotData(N, tAxis, y, ySq, yI);
  plotData(N/2, tAxis, yI);
}

void sineIntegral()
{
  double test = rsSineIntegral(3.0);
  double tMin = -20.0;
  double tMax = +20.0;
  static const int N = 1000;
  double t[N];
  double y[N];
  RAPT::rsArrayTools::fillWithRangeLinear(t, N, tMin, tMax);
  for(int n = 0; n < N; n++)
    y[n] = rsSineIntegral(t[n]);
  plotData(N, t, y);
}


double lq(double x)
{
  if( fabs(x) < EPS )
    return -1.0;
  //rsComplexDbl xc(x, 0);
  //return (rsLogC(1.0-xc)/rsLogC(1.0+xc)).re;
  return log(1-x) / log(1+x);
}
void logarithmQuotient()
{
  // some investigations for the function lq(x) := log(1-x) / log(1+x)
  double xMin = -0.99;
  double xMax = +0.99;
  static const int N = 1000;
  double x[N];
  double y[N];
  RAPT::rsArrayTools::fillWithRangeLinear(x, N, xMin, xMax);
  RAPT::rsArrayTools::applyFunction(x, y, N, &lq);

  GNUPlotter p;
  p.plotFunctions(N, xMin, xMax, &lq);
}

void stirlingNumbers()
{
  static const int nMax = 10;
  int n, k;

  int s[nMax+1][nMax+1];   // Stirling numbers

  int **tmp;
  RAPT::rsArrayTools::allocateSquareArray2D(tmp, nMax+1);
  rsStirlingNumbersFirstKind(tmp, nMax);
  for(n = 0; n <= nMax; n++)
  {
    for(k = 0; k<= nMax; k++)
      s[n][k] = tmp[n][k];
  }
  RAPT::rsArrayTools::deAllocateSquareArray2D(tmp, nMax+1);
}


/** Returns the Bernoulli number for the given n (using the convention that B1 = -1/2). */
double rsBernoulliNumber(rsUint32 n)
{
  if( n == 1 )
    return -0.5;
  if( rsIsOdd(n) )
    return 0.0;
  double *A = new double[n+1];
  for(rsUint32 m = 0; m <= n; m++)
  {
    A[m] = 1.0 / (m+1);
    for(rsUint32 j = m; j >= 1; j--)
      A[j-1] = j * (A[j-1] - A[j]);
  }
  double result = A[0];
  delete[] A;
  return result;
  // hmm - this algorithm seems to be quite imprecise numerically (it's taken form wikipedia)
  // maybe we should use this algorithm with exact representations of rational numbers
}
void bernoulliNumbers()
{
  static const int nMax = 10;
  double b[nMax+1]; 
  for(int n = 0; n <= nMax; n++)
    b[n] = rsBernoulliNumber(n);
}

void sequenceSquareRoot()
{
  // The overall goal is to decompose the sequence C into two sequences q and r such that 
  // q^2 + r^2 = C (whe a sequence squared is understood as convolution with itself).
  // ideally, we yould want to get back q = a, r = b but due to ambiguities this is not possible.

  
  // two arbitrary sequences to be "squared" and added:
  static const int N = 5; // upper index in the sequence a0, a1, a2, ..., aN - likewise for b
  
  double a[N+1] = {2, 3, 4, 5, 6, 7};
  double b[N+1] = {9, 8, 7, 6, 5, 4};
  //double a[N+1] = {2, -3,  2, -1, 5,  3};
  //double b[N+1] = {3,  1, -4, -2, 2, -4};
  //double b[N+1] = {0,  0, 0, -2, 2, -4};  // not really arbitrary

  /*
  // two arbitrary sequences to be "squared" and added:
  static const int N = 4; // upper index in the sequence a0, a1, a2, ..., aN - likewise for b
  double a[N+1] = {2, -3,  2, -1, 5};
  double b[N+1] = {3,  1, -4, -2, 2};
  */


  // the "square-sequences": and its sum:
  double A[2*N+1];
  double B[2*N+1];
  RAPT::rsArrayTools::convolve(a, N+1, a, N+1, A);
  RAPT::rsArrayTools::convolve(b, N+1, b, N+1, B);

  // sum of the two square-sequences: 
  double C[2*N+1];
  RAPT::rsArrayTools::add(A, B, C, 2*N+1);

  // retrieve a sequence q that, when convolved with itself, gives a sequence in which the first
  // N values agree with C:
  double q[N+1];
  RAPT::rsArrayTools::sequenceSqrt(C, 2*N+1, q);

  // get Q, the full square sequence of q - it should agree with C up to index N:
  double Q[2*N+1];
  RAPT::rsArrayTools::convolve(q, N+1, q, N+1, Q);

  // get the difference D between C and Q:
  double D[2*N+1];
  RAPT::rsArrayTools::subtract(C, Q, D, 2*N+1);

  // find a sequence r that, when squared, gives D
  double r[N+1];
  int m = (N+1)/2;        // number of leading zeros in r
  RAPT::rsArrayTools::fillWithZeros(r, m);
  RAPT::rsArrayTools::sequenceSqrt(&D[2*m], N, &r[m]); 

  /*
  double Ds[2*N+1]; // D, shifted and zero padded
  rsCopyBuffer(&D[2*m], Ds, N);
  rsFillWithZeros(&Ds[N], N+1);
  rsSequenceSqrt(Ds, N, &r[m]); 
  */

  double R[2*N+1];
  RAPT::rsArrayTools::convolve(r, N+1, r, N+1, R);  // R should match D
                                  // ...yes but only the 1st 3 values - which makes sense because
                                  // r has only 3 degrees of freedom - hmmm.....
  int dummy = 0;
}




template<class TC, class TR> // TC: coefficients, TR: result
TR rsSolveDoubleSqrtEquation(TC p[10], TR x0, double tol = 1.e-13)
{
  // Finds a numeric solution to the equation:
  //
  // 0 = p0 + p1*x + p2*sqrt(p3+p4*x+p5*x^2) + p6*sqrt(p7+p8*x+p9*x^2)
  //
  // via Newton iteration using an initial guess x0. An equation of this type arises in the 
  // solution of a system of 2 conic equations.

  TR x, dx, f, fp, s1, s2;
  x = x0;
  do
  {
    s1  = sqrt(p[3] + p[4]*x + p[5]*x*x);
    s2  = sqrt(p[7] + p[8]*x + p[9]*x*x);
    f   = p[0] + p[1]*x + p[2]*s1 + p[6]*s2;
    fp  = p[1] + p[2]*(p[4]+2*p[5]*x)/(2.0*s1) + p[6]*(p[8]+2*p[9]*x)/(2.0*s2);
    dx  = -f/fp;
    x  += dx;
  } while(abs(dx/x) > tol);
  return x;

  // maybe include an iteration counter and leave loop if > maxNumIterations
}
template<class TC, class TR> // TC: coefficients, TR: result
void rsSolveConicSystem(TC a[6], TC b[6], TR *x, TR *y, rsUint8 n, TR x0 = 0)
{
  // We numerically solve a system of two quadratic equations of two variables x,y of the form:
  //
  // 0 = a0 + a1*x + a2*y + a3*x^2 + a4*x*y + a5*y^2
  // 0 = b0 + b1*x + b2*y + b3*x^2 + b4*x*y + b5*y^2
  //
  // The two equations each describe a conic section and these conic sections may intersect at
  // various points...
  // There are 3 places where we have to choose either a + or - sign when adding a value under a
  // square-root which gives rise to 8 possible solutions for the (x,y)-pair - one for each of the 
  // 8 possible combinations of +/-. The value n in the range 0..7 decides, which solution should 
  // be picked. For the internal Newton-iteration based solver, we need an initial guess for x 
  // given by x0.
  // Solving both equations for y gives:
  //
  // y = (-a2-a4*x +- sqrt((a2^2-4*a5*a0) + (2*a4*a2-4*a5*a1)*x + (a4^2-4*a5*a3)*x^2)) / (2*a5)
  // y = (-b2-b4*x +- sqrt((b2^2-4*b5*b0) + (2*b4*b2-4*b5*b1)*x + (b4^2-4*b5*b3)*x^2)) / (2*b5)
  //
  // Subtracting the 1st from the 2nd gives an equation suitable for rsSolveDoubleSqrtEquation
  //
  // 0 = p0 + p1*x + p2*sqrt(p3+p4*x+p5*x^2) + p6*sqrt(p7+p8*x+p9*x^2)
  //
  // using
  //
  // p0 = ...
  //

  // get the signs for the square-roots: 0:+++, 1:++-, 2:+-+, 3:+--, 4:-++, 5:-+-, 6:--+, 7:---
  double s1, s2, s3;
  s3   = 1 - 2 * ((n & 1) >> 0);
  s2   = 1 - 2 * ((n & 2) >> 1);
  s1   = 1 - 2 * ((n & 4) >> 2);

  // compute coefficients for rsSolveDoubleSqrtEquation:
  TC ka, kb, p[10];
  ka   = 1 / (2*a[5]);
  kb   = 1 / (2*b[5]);
  p[0] = ka*a[2] - kb*b[2];
  p[1] = ka*a[4] - kb*b[4];
  p[2] = s1 * ka;
  p[3] =   a[2]*a[2] - 4*a[5]*a[0];
  p[4] = 2*a[4]*a[2] - 4*a[5]*a[1];
  p[5] =   a[4]*a[4] - 4*a[5]*a[3];
  p[6] = s2 * kb;
  p[7] =   b[2]*b[2] - 4*b[5]*b[0];
  p[8] = 2*b[4]*b[2] - 4*b[5]*b[1];
  p[9] =   b[4]*b[4] - 4*b[5]*b[3];

  // solve:
  *x   = rsSolveDoubleSqrtEquation(p, x0);

  *y   = ka * (-a[2]-a[4]*(*x) + s3 * sqrt(p[3] + p[4]*(*x) + p[5]*(*x)*(*x)));
  // i think, the computation of y is wrong -  mabye we need to ty both signs for the sqrt and 
  // check each solution - only one of the two really is a solution


  //*y = kb * (-b[2]-b[4]*(*x) + s3 * sqrt(p[7] + p[8]*(*x) + p[9]*(*x)*(*x))); 
  // should be the same - maybe we can somehow choose the numerically more precise one based on some
  // criteria?


  // maybe a better algorithm for the case of real solutions is this:
  // https://en.wikipedia.org/wiki/Conic_section#Intersecting_two_conics
  // http://math.stackexchange.com/questions/316849/intersection-of-conics-using-matrix-representation
}

// evaluates a0 + a1*x + a2*y + a3*x^2 + a4*x*y + a5*y^2
template<class TC, class TR> // TC: coefficients, TR: result
TR rsEvaluateConic(TC a[6], TR x, TR y)
{
  return a[0] + a[1]*x + a[2]*y + a[3]*x*x + a[4]*x*y + a[5]*y*y;
}

// checks, if a0 + a1*x + a2*y + a3*x^2 + a4*x*y + a5*y^2 is 0, same for b-coeffs
template<class TC, class TR>
bool checkConicResult(TC a[6], TC b[6], TR x, TR y, double tol)
{
  TR   f;          // function value
  bool res;        // test result
  res  = true;
  f    = rsEvaluateConic(a, x, y);
  res &= abs(f) <= tol;
  f    = rsEvaluateConic(b, x, y);
  res &= abs(f) <= tol;
  return res;
}

void conicSystem()
{
  // move to unit tests...

  bool   res = true;        // test result
  double tol = 1.e-13;      // tolerance
  double p[10];             // parameters
  complex<double> x,  y;    // computed solution
  complex<double> xt, yt;   // target solution (found with wolfram alpha)

  // test rsSolveDoubleSqrtEquation:
  RAPT::rsArrayTools::fillWithValue(p, 10, 1.0);
  p[1] =  4;
  p[6] =  2;
  xt   = -1;
  x    = rsSolveDoubleSqrtEquation(p, 0.0, tol);
  res &= abs(x-xt) <= tol;
  p[8] = 2;
  xt   = -(1 + sqrt(3./35))/2; // 0.64639
  x    = rsSolveDoubleSqrtEquation(p, 0.0, tol);
  res &= abs(x-xt) <= tol;

  // test rsSolveConicSystem:
  complex<double> j(0.0, 1.0);   // imaginary unit
  double a[6], b[6];             // parameters
  RAPT::rsArrayTools::fillWithValue(a, 6, 1.0);
  RAPT::rsArrayTools::fillWithValue(b, 6, 1.0);
  b[4] = 2;

  rsSolveConicSystem(a, b, &x, &y, 4, complex<double>(0.1));
  xt = (-1.0 + j*sqrt(3)) / 2.0;
  yt = 0.0;
  res &= abs(x-xt) <= tol;
  res &= abs(y-yt) <= tol;
  res &= checkConicResult(a, b, x, y, tol);

  rsSolveConicSystem(a, b, &x, &y, 6, complex<double>(0.1));
  xt = (-1.0 - j*sqrt(3)) / 2.0;
  yt = 0.0;
  res &= abs(x-xt) <= tol;
  res &= abs(y-yt) <= tol;
  res &= checkConicResult(a, b, x, y, tol);

  // test:
  //rsSolveConicSystem(a, b, &x, &y, 0, complex<double>(-0.5, -0.5)); // no convergence


  // these give wrong results - how is this possible?:
  rsSolveConicSystem(a, b, &x, &y, 7, complex<double>(0.1, 0.1));
  res &= checkConicResult(a, b, x, y, tol);

  rsSolveConicSystem(a, b, &x, &y, 5, complex<double>(0.1, 0.1));
  res &= checkConicResult(a, b, x, y, tol);

  // we get two wrong results and two right results are missing - might this be a bug?

  // solutions 4,5,6,7 work.

  // todo: use the conics:
  // x^2    -  y^2 + xy + x   + 2y - 1 = 0    // hyperbola
  // x^2/2  + 2y^2 + xy + x/4 +  y - 2 = 0    // ellipse
  // should have 4 real solutions, if the -2 in the 2nd equations is replaced by -1/2, there are 2 
  // real solutions and with -1/8 no real solutions
  // there seem to ba always 4 solutions in general, some of which may be complex

  // maybe we should cast the problem in the more standardized form of conic section equations:
  // A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0

  // Factoring a degenerate conic into its 2 lines (required for the algorithm for the intersection
  // of 2 conics on wikipedia):
  // A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = (a*x + b*y + c) * (d*x + e*y + f)
  // A = a*d, B = a*e + b*d, C = b*e, D = a*f + c*d, E = b*f + c*e, F = c*f
  // -> solve(A=a*d,B=a*e+b*d,C=b*e,D=a*f+c*d,E=b*f+c*e,F=c*f;a,b,c,d,e,f) -> no result
  // use d = A/a, e = C/b, f = F/c in the other 3 equations:
  // B = a*C/b + b*A/a, D = a*F/c + c*A/a, E = b*F/c + c*C/b
  // solving each of these equations for a,b or c leads to 2 possible solutions (these equations
  // are quadratic equations in disguise), example solve(B = a*C/b + b*A/a, a)
  // the 1st one can be rewritten as: 0 = a^2 - (b*B/C)*a + b^2/C

  // Another possible algorithm to find the solution - write down the 2 conics as:
  // A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0
  // a*x^2 + b*x*y + c*y^2 + d*x + e*y + f = 0
  // For B=0, a general way (?) based on completing the square to solve this can be found here:
  // http://www.ck12.org/book/CK-12-Algebra-II-with-Trigonometry/section/10.5/
  // ...but, we have, in general B != 0 and b != 0, but we can solve the two equations for the
  // xy cross-term:
  // x*y = -(A*x^2 + C*y^2 + D*x + E*y + F) / B
  // x*y = -(a*x^2 + b*y^2 + d*x + e*y + f) / b
  // and plug the 1st solution in the 2nd original equation and vice versa and obtain a system with
  // modified coefficients in which in both equations the cross-term is zero, so the technique of 
  // completing the square can be applied to both equations (i'm actually not sure, if this really 
  // works out, it's just an idea).

  int dummy = 0;
}

void logisticMapNoise()
{
  // We generate a pseudo random sequence of numbers using the logistic map: xNew = r * x * (1-x)
  // where r is the growth rate parameter.

  int    N  = 8000;    // number of values
  double r  = 3.999;   // growth rate (should be < 4)
  double x0 = 0.5;     // initial x-value

  // create the time series:
  vector<double> x(N);
  x[0] = x0;
  for(int n = 1; n < N; n++)
    x[n] = r * x[n-1] * (1 - x[n-1]);

  writeToMonoWaveFile("LogisticMapNoise.wav", &x[0], N, 44100, 16);

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0]);
  plt.plot();
}

void variousFunctions()
{
  // For plotting various functions of the type y = f(x) to see if they look as expected.
}

inline double sin2(double x)
{
  return sin(x);
}

// move to rsRootFinder:
template<class T>
T findLeftBracket(const std::function<T(T)>& f, T y, T xL = T(0), T d = T(1))
{
  while(f(xL) > y) { xL -= d; d *= 2; }
  return xL;
}
template<class T>
T findRightBracket(const std::function<T(T)>& f, T y, T xR = T(0), T d = T(1))
{
  while(f(xR) < y) { xR += d; d *= 2; }
  return xR;
}
// passed xR is initial guess, d is the initial increment

template<class T>
inline std::function<T(T)> rsInverse(const std::function<T(T)>& f)
{
  std::function<T(T)> fi;  // inverse function
  fi = [=](T y) {
    // wrap these 3 lines into rsRootFinder::findRoot(f, y)
    T xL = findLeftBracket( f, y);
    T xR = findRightBracket(f, y);
    T x = rsRootFinder<T>::bisection(f, xL, xR, y); // use better algo
    return x;
  };
  return fi;

  // todo: 
  // -make it work for monotonically decreasing functions (it currently works only for 
  //  monotonically increasing function)
  // -allow the user to restrict the domain such that we can pick a chunk of a function where it
  //  is monotonic in case it is not monotonic over the whole real number line
  // -figure out, what it does when there are multiple solutions or no solution
  //  and make the behavior well defined in these cases - or maybe we should assume monotonic
  //  functions - otherwise, they are not uniqely invertible anyway
}

// computes the definite integral of f from a to b using the trapezoidal rule with N steps
// this is a formula for a closed interval, i.e. it evaluates the function at the endpoints a and b
template<class T>
T rsIntegral(const std::function<T(T)>& f, T a, T b, int N = 16)
{
  // make sure a < b:
  if(b == a) return T(0);
  if(b <  a) return -rsIntegral(f, b, a, N); // swap limits and invert sign

  T fL  = f(a);
  T y   = 0;
  T h   =  (b-a)/N;
  for(int i = 1; i <= N; i++) {
    T x  = a + i*h; // maybe optimize: x += h (avoids int/float conversion but accumulates error)
    T fR = f(x);
    y   += 0.5 * h * (fL + fR);
    fL   = fR;
  }
  return y;

  // very preliminary - uses fixed N - todo: refine until desired accuracy is obtained
  // NR has code that can re-use the results of previous computations
  // mayb implement the Romberg integration algorithm from NR
}

// Computes the indefinite integral "F" of the function "f" with lower integration limit "a" and 
// integration constant "c" where c is the value of the antiderivative at x=a, such that F(a) = c
template<class T>
std::function<T(T)> rsAntiDerivative(const std::function<T(T)>& f, T a, T c, int N)
{
  return [=](T x) { return rsIntegral(f, a, x, N) + c; };
}
// todo: remove N parameter - either replace with an accuracy parameter or (better) let the 
// function figure out the number of steps that gives maxium accuracy or even better make accuracy
// parameter optional


void functionOperators()
{
  // For testing various operators that take a function as input and return another functions as 
  // output

  std::function<double(double)> f;               // holds our function
  //double y;

  // periodicized polynomial:
  f = [=](double x) { return x*(x-PI)*(x+PI); }; // polynomial
  f = rsMakePeriodic(f, -PI, PI);                // sine-ish
  rsPlotFunction(f, -10.0, +10.0, 1000);

  // periodicized indentity:
  f = [=](double x) { return x; };               // identity
  f = rsMakePeriodic(f, -3.0, 2.0);              // sawtooth
  rsPlotFunction(f, -10.0, +10.0, 1000);

  // derivative:
  f = [=](double x) { return sin(2*x); };        // f(x)  =   sin(2*x)
  f = rsDerivative(f, 0.01);                     // f'(x) = 2*cos(2*x)
  rsPlotFunction(f, -10.0, +10.0, 1000);

  // antiderivative:
  f = [=](double x) { return sin(2*x); };        // f(x) = sin(2*x)
  f = rsAntiDerivative(f, 0.0, -0.5, 128);       // F(x) = -1/2 * cos(2*x) + c
  rsPlotFunction(f, -6.0, +6.0, 1000);

  // inverse function:
  f = [=](double x) { return x*x*x; };           // f(x)    = x^3
  f = rsInverse(f);                              // f^-1(x) = cubeRoot(x)
  rsPlotFunction(f, -5.0, +5.0, 1000);
  // inversion needs some more care - we need to support (monotonically) decreasing functions as 
  // well (we should automatically determine which case it is) - maybe have two versions
  // rsInverseIncreasing/Decreasing and let rsInverse dispatch bewteen these two
  // maybe the operator should also get a domain xL...xR as input to support functions that are
  // not monotonic over all real numbers - for example: rsInverse(sin, -pi/2, +pi/2) should give
  // the main branch of arcsin - it should trigger an assert when called with values that are not
  // within the range -1...+1

  // this will probably not work (hang):
  //f = [=](double x) { return -x*x*x; };           // f(x)   = -x^3
  //f = rsInverse(f);                              // f^-1(x) = -cubeRoot(x)
  //rsPlotFunction(f, -5.0, +5.0, 1000); 
  // yep - hangs because the function is decreasing and the bracket-search assumes an increasing 
  // function

  // even part:
  f = [=](double x) { return exp(x); };          // f(x)  = e^x
  f = rsEvenPart(f);                             // fe(x) = (e^x + e^-x) / 2 = cosh(x)
  rsPlotFunction(f, -5.0, +5.0, 1000);

  // odd part::
  f = [=](double x) { return exp(x); };          // f(x)  = e^x
  f = rsOddPart(f);                              // fo(x) = (e^x - e^-x) / 2 = sinh(x)
  rsPlotFunction(f, -5.0, +5.0, 1000);

  //f = &sin2;  // works
  //f = sin2;  // works
  //f = &sin;   // fails - why? has it to do with sin having many overloaded versions?

  // is it possible to find numerical solutions of differential equations using the derivative
  // opertator? try: y + a*y' + b*y'' = g(x) ...maybe throw a root-finder at...
  // ....‰‰‰hh - no - i think, that doesn't work if we don't know the input function y(x) already
  // and that's exactly what we try to find
}


// maybe move to vairiousFunctions
void hyperbolicFunctions()
{  
  GNUPlotter plt;
  int N = 1000;
  plt.addDataFunctions(N, -2.0, +2.0, rsSinh, rsCosh, rsTanh);
  //plt.addDataFunctions(N, -2.0, +2.0, rsCoth, rsCsch, rsSech);
  plt.plot();

  ///*
  //// test rsSinhCosh:
  //static const int N = 1000;
  //double x[N], ys[N], yc[N];
  //rsFillWithRangeLinear(x, N, -3.0, +3.0);
  //for(int n = 0; n < N; n++)
  //  SinhCosh(x[n], &ys[n], &yc[n]);
  //p.plotData(N, x, ys, yc);
  //*/

  //// test inverse functions:
}

void bigFloatErrors()
{
  // idea: the rsBigFloat should have a subclass rsBigFloatWithErrorBound that carries information 
  // about the worst-case error with it. it has an errorBound variable (an integer in units of 
  // "ULP"="unit-in-the-last-place") which is adjusted in each arithmetic operation. on 
  // contruction/initialization it's initiliazed to zero (or to 1, if rounding occurs during 
  // initialization). 

  // mulitplication: let d1, d2 be the relative error bounds of both operands as real value 
  // (not in ulps). then, the maxmimum error after multiplication d3 becomes:
  // d3 = (1+d1)*(1+d2)-1 (we assume that instead of multiplying 1*1, we multiply the approximate
  // values (1+d1), (1+d2), using 1 as the real value represents the worst case (1 is the smallest
  // nonzero digit in any base) and for d1,d2 > 0 we have (1+d1)*(1+d2) - 1 > 1 - (1-d1)(1-d2), so 
  // the worst case occurs when bot errors are upward. so:
  // d3 =  d1 + d2 + d1*d2. if u1, u2, u3 are the corresponding errors in ulps, we use 
  // d1 = u1 / (2*B^(N-1)) (for base B), for the resulting multiplication error in ulps, we get:
  // u3 = u1 + u2 + u1*u2*B^(3*N-1) / (2*B^(N^2+1))  // verify this
  // if we round the multiplication result, we increase the error by another ulp, so
  // u3 = 1 + u1 + u2 + ceil(k*u1*u2) where k = B^(3*N-1) / (2*B^(N^2+1))

  // division: the worst case occurs when error of the 1st operand is upward and the error
  // of the 2nd operand is downward: (1+d1)/(1-d2) - 1 > 1 - (1-d1)/(1+d2), so: 
  // d3 = (1+d1)/(1-d2) - 1 = d1 - d1*d2
  // u3 = u1 - k*u1*u2 - because the 2nd term is subtracted, the maximum absolute value of this is:
  // u3 <= max(u1, k*u1*u2) ...shoudln't it be max(u1, abs(u1 - k*u1*u2))?
  // so we may set after rounding:
  // u3 = 1 + max(u1, k*u1*u2)
    
  // addition/subtraction: worst case when both errors are in the same direction, then
  // u3 = u1 + u2 (before renormalization)

  // normalization:
  // if we shift the fraction to the left by N ( >= 1) digits (therby introducing N zeros at the 
  // right) we must increase the error u in ULPs by B^(N-1)

  // TODO: verify these calculations - this has been done quickly and i'm not completely sure, that
  // this analysis is valid...

  // base conversion: let x = Bx^Ex * (x1/Bx + x2/Bx^2 + x3/Bx^3 + ...)
  //                  and y = By^Ey * (y1/By + y2/By^2 + y3/By^3 + ...) 
  // and the goal is to find the new exponent By and the new fraction-digits y1, y2, y3
  // from the old exponent Ex and the old fraction digits x1, x2, x3,... and the old and new base 
  // Bx. By
  // we require x = y, so:
  // (y1/By + y2/By^2 + y3/By^3 + ...) = (Bx^Ex/By^Ey) * (x1/Bx + x2/Bx^2 + x3/Bx^3 + ...)
  // this means, we can use the rsDigitArithmetic::changeBaseForFraction function to convert the 
  // fractions and afterwards multiply the result by (Bx^Ex/By^Ey), represented as rsBigFloat
  // in the new base. but we could also first multiply by (Bx^Ex/By^Ey) in the old base and then
  // convert the digits. perhaps, it's best to do the multiplication in the representation with
  // the higher precision (which depends on the base as well as the number of digits - look at the
  // value of the ULP.

  // write a function that returns the largest integer that can be represented exactly
  // static rsBigInt rsBigFloat::largestExactInteger(rsUint32 numDigits, rsUint64 base);
  // this can be computed as follows: us a temporary rsBigInt in the passed base with the same 
  // number of digits as the fraction of the rsBigFloat - in this representation, the largest 
  // exact integer is given by an rsBigInt where all digits are equal to base-1. for example, for 
  // numDigits=2, base=10, the largest exact integer is 99_10 for numDigits=2, base=8, its 
  // 77_8 = 72_10. convert that temporary value into an rsBigInt in the default base and return 
  // that value.

  int dummy = 0;
}

void primeRecursion()
{
  //static const rsUint32 N = 1000;
  static const rsUint32 N = 168; // number of primes to generate
  rsUint32 p[N], tp[N];          // primes generated by recursion and true primes

  int dummy;

  rsFillPrimeTable(tp, N);

  // init first few primes, d and Pd:
  rsUint32 np = 11;
  RAPT::rsArrayTools::copy(tp, p, np);
  rsUint32 d  = 7;
  rsUint32 Pd = 2*3*5*7;
  while( np < N )
  {
    p[np] = (d*p[np-d-1] + Pd) / d;
    np++;
    dummy = 0;
  }


  /*
  rsUint32 i = 1;
  rsUint32 P;
  rsUint32 n;
  while( true )
  {
    P = rsProduct(tp, i);
    n = P + tp[i];
    i++;
  }
  */

  /*
  // find multiples of 11 which are not divisible by 2,3,5 or 7:
  rsUint32 m11 = 11;
  rsUint32 r11 = 1;
  rsUint32 f[100], m[100];
  rsUint32 i = 0;
  while(true)
  {
    if( m11 % 2 == 0 || m11 % 3 == 0 || m11 % 5 == 0 || m11 % 7 == 0 ) // rsIsDivisibeByAnyOf
    {

    }
    else
    {
      f[i] = r11;
      m[i] = m11;
      i++;
    }
    m11 += 11;
    r11++;
  }
  */


  dummy = 0;

/*
idea:
-cross multiples of 2: 
 stepsize 2, start: 2*2=4
-cross multiples of 3 which are no multiple of 2
 stepsize: 2*3=6, start: 3*3=9
-cross multiples of 5 which are no multiple of 2 or 3
 stepsize: 2*3*5=30, start: 5*5=25, 7*5=35
-cross multiples of 7 which are no multiple of 2,3 or 5:
 stepsize: 2*3*5*7=210, 
 start: 7*7=49, 11*7=77, 13*7=91, 17*7=119, 19*7=133, 
        23*7=161, 29*7=203, 31*7=217 
		the next candidates for the start are already covered: 
		37*7=259= 49+210, 41*7=287= 77+210, 43*7=301= 91+210, 
		47*7=329=119+210, 53*7=371=161+210, 59*7=413=203+210
		61*7=427=217+210, 
		67*7=469=259+210  -> here, it repeats again
		
this suggests, that we can predict the primes past 31 from 
previous known primes as:
P7 := product(primes-up-to-7) = 2*3*5*7 = 210
37 = 259/7 = ( 49+P7)/7 = (7* 7 + P7) / 7
41 = 287/7 = ( 77+P7)/7 = (7*11 + P7) / 7
43 = 301/7 = ( 91+P7)/7 = (7*13 + P7) / 7
47 = 329/7 = (119+P7)/7 = (7*17 + P7) / 7
53 = 371/7 = (161+P7)/7 = (7*23 + P7) / 7
59 = 413/7 = (203+P7)/7 = (7*29 + P7) / 7

does this suggest a recursion to generate primes????!!!!
that would be too good to be true!!!

let d be some prime, for example d = 7, 
Pd = product(primes-up-to-d) = 2*3*5*7=210, 
then: p[k] = (d*p[k-d] + Pd) / d  ?
*/
}

void primeSieveSchmidt1()
{  
  static const rsUint32 N = 1000;
  //bool isPrime[N+1];  // later
  rsUint32 a[N+1];
  rsUint32 p[168];      // there are 168 primes up to N=100

  RAPT::rsArrayTools::fillWithRangeLinear(a, N+1, rsUint32(0), rsUint32(N));

  rsUint32 i, j, k, k5, np, np1, pk, pk1;  
  int dummy;

  a[0] = 0;
  a[1] = 0;
  for(i = 4; i < N; i += 2)
    a[i] = 0;
  for(i = 9; i < N; i += 6)
    a[i] = 0;
  p[0] = 2; p[1] = 3; p[2] = 5; p[3] = 7;
  np = 4;
  k  = 1; 
  k5 = 2;  // index of prime "5"
  pk = p[k];
  while( true )  // we need some stopping criterion
  {
    // remember values from previous iteration:
    pk1 = pk;
    np1 = np;

    // step 1 - fetch new prime:
    k   = k+1;
    pk  = p[k];

    // step 2 - collect new primes:
    for(i = pk1*pk1; i < pk*pk; i++)
    {
      if( a[i] != 0 )
      {
        p[np] = a[i]; // later maybe use i itself (when using a bool-array), maybe we can also
                      // use i+start at some stage when the buffer does not start at 1
        np++;
      }
    }

    // step 3 - cross multiples of p[k]:
    for(i = k; i < np; i++)
    { 
      a[pk*p[i]] = 0;
      dummy = 0;
    }

    // step 4 - cross new available multiples of primes below pk
    for(j = k5; j < k; j++)
    {
      for(i = np1; i < np; i++)
        a[p[j]*p[i]] = 0;
      dummy = 0;
    }
    dummy = 0;
  }
  dummy = 0;

/*
(not yet working) algorithm for crossing every nonprime once (using 1-based array indexing):
write down natural numbers up to N an array a=1...N
cross 1 (it's not a prime)
cross all even numbers(start at 4, use stepsize of 2)
cross all odd mutiples of 3 (start at 9, use stepsize 6)
init prime array p with all primes <= 3^2: p = 2,3,5,7
init ip=2
init current prime pk=p[ip]=3
init number of known primes np=4

loop:
1: fetch next prime pk:
   pkOld = pk
   ip    = ip+1
   pk    = p[ip]
2: collect all nonzero values in a from pkOld^2 to pk^2-1 and 
   append them to our array of primes p, update the number of
   known primes np accordingly
3: cross pk*r where r is any of the primes from pk upward   
4: for all primes q from 5 up to but not including pk, cross q*r 
   in array a where r is any of the numbers that have just been 
   added to p in step 2
5: go back to 1  

1: fetch pk=5
2: append 11,13,17,19,23 to p
3: cross 5*5,5*7,5*11,5*13,5*17,5*19,5*23 
4: nothing to do

1: fetch pk=7
2: append 29,31,37,41,43,47 to p
3: cross 7*7,7*11,7*13,7*17,7*19,7*23,7*29,7*31,7*37,7*41,7*43,7*47 
4: cross 5*29,5*31,5*37,5*41,5*43,5*47
   
1: fetch pk=11
2: append 53,59,61,67,71,73,79,83,89,97,101,103,107,109,113 to p 
3: cross 11*11,11*13,11*17,...,11*109,11*113
4: cross 5*53,5*59,...,5*113
   cross 7*53,7*59,...,7*113
   
1: fetch pk=13
2: append 127,...,167
3: cross 13*13,...,13*167
4: cross  5*127,..., 5*167
   cross  7*127,..., 7*167
   cross 11*127,...,11*167
   
we need to establish that at any time, the greatest prime in the 
array times 5 is larger than the next prime squared, otherwise we 
cannot guarantee that all nonprimes up to pk^2 have been crossed 
already. an upper bound for the next prime after pk is 2*pk-1 
(http://en.wikipedia.org/wiki/Bertrand%27s_postulate). we have 
5*pk^2 > (2*pk-1)^2 = 4*pk^2 - 4*pk + 1, so the desired requirement 
is fulfilled. (we could actually have used the tighter formulation
with an upper bound of 2*pk-3)
   
loop-invariant: before and after each iteration, the array a has 
crossed out all nonprimes up to pk^2 and p contains all primes up
to pk^2 ...is that true? - noo, it misses to cross nonprimes with 
multiple factors (except for simple squares). the 1st such values 
are 5^2*5=125,5^2*7=175. an inner loop in steps 3 and 4 over the 
exponents could help, but then we would still miss to cross out 
numbers with more than one multiple factor like 5^2*7^2,... we 
could add another loop nesting level which runs over the exponent
of the second factor, but that would cover only nonprimes with 2
multiple factors and without any other factors: we would still miss 
5^2*7^2*11,...

maybe we sohould modify the algorithm as follows:
loop:
1: fetch next prime pk
2: collect all nonzero values in a from pkOld^2 to pk^2-1 and 
   append them to our array of primes p, update the number of
   known primes np accordingly
3: for all primes q with 5 <= q <= pk, and all multipliers r (not 
   required to be prime) with pkOld+1 <= r < pk, cross all 
   multiples q*r
4: go back to 1  

step 3 should now also cover composite numbers with multiple factors,
but we are not guaranteed anymore to cross out each nonprime exactly
once (although, redundant crossings occur more rarely than in the
sieve of Eratosthenes)
*/
}


void primeSieveSchmidt2()
{  
  static const rsUint32 N = 1000;
  //bool isPrime[N+1];  // later
  rsUint32 a[N+1];
  rsUint32 p[168];      // there are 168 primes up to N=100

  RAPT::rsArrayTools::fillWithRangeLinear(a, N+1, rsUint32(0), rsUint32(N));

  rsUint32 i, j, k, k5, np, np1, pk, pk1, q, r;  
  int dummy;

  a[0] = 0;
  a[1] = 0;
  for(i = 4; i < N; i += 2)
    a[i] = 0;
  for(i = 9; i < N; i += 6)
    a[i] = 0;
  p[0] = 2; p[1] = 3; p[2] = 5; p[3] = 7;
  np = 4;
  k  = 1; 
  k5 = 2;  // index of prime "5"
  pk = p[k];
  while( true )  // we need some stopping criterion
  {
    // remember values from previous iteration:
    pk1 = pk;
    np1 = np;

    // step 1 - fetch new prime:
    k   = k+1;
    pk  = p[k];

    // step 2 - cross out more non-primes:
    //for(j = k5; j <= k; j++)
    for(j = k5; j < np; j++)
    {
      q = p[j];
      //for(r = pk1+1; r <= pk; r++)
      for(r = pk1; r <= pk; r++)
        a[q*r] = 0;
        // can be done without multiplications in the loop: start at q*(pk1+1) with stepsize q
    }

    // step 3 - collect new primes:
    for(i = pk1*pk1; i < pk*pk; i++)
    {
      if( a[i] != 0 )
      {
        p[np] = a[i]; // later maybe use i itself (when using a bool-array), maybe we can also
                      // use i+start at some stage when the buffer does not start at 1
        np++;
      }
    }
    dummy = 0;
  }
  dummy = 0;

/*
1: fetch next prime pk
2: for all primes q with 5 <= q <= pk, and all multipliers r (not 
   required to be prime) with pkOld+1 <= r <= pk, cross all 
   multiples q*r
3: collect all nonzero values in a from pkOld^2 to pk^2-1 and 
   append them to our array of primes p, update the number of
   known primes np accordingly
4: go back to 1

1: fetch pk=5
2: q=5: cross 5*4,5*5
3: append 11,13,17,19,23 to p

1: fetch pk=7
3: q=5: cross 5*6,5*7
   q=7: cross 7*6,7*7;
2: append 29,31,37,41,43,47 to p

1: fetch pk=11
2: q= 5: cross  5*8, 5*9, 5*10, 5*11
   q= 7: cross  7*8, 7*9, 7*10, 7*11
   q=11: cross 11*8,11*9,11*10,11*11
3: append 53,59,61,67,71,73,79,83,89,97,101,103,107,109,113 to p 

1: fetch pk=13
2: q= 5: cross  5*12, 5*13
   q= 7: cross  7*12, 7*13
   q=11: cross 11*12,11*13
3: append 127,...,167 to p

optimization: in step 2, we could use an increment of 2 to avoid crossing out even values 
again and maybe with 2 passes with an increment of 6, we could avoid the multiples of 3 as 
well

no: misses to cross 65,85,91,95,115,119 when pk=11, perhaps, in step 2, q should not go only up 
to pk but up to p[np-1]

baah! that all doesn't work, in the current version, 125 is missed
*/
}

void primeSieveAtkin()
{
  static const rsUint32 limit = 1000;
  rsUint32 r = rsIntSqrt(limit);
  
  bool isPrime[limit+1];
  rsUint32 i;
  for(i = 0; i <= limit; i++)
    isPrime[i] = false;

  rsUint32 x, y, n;
  for(x = 1; x <= r; x++)
  {
    for(y = 1; y <= r; y++)
    {
      n = 4*x*x + y*y;
      if( n <= limit && (n % 12 == 1 || n % 12 == 5) )
        isPrime[n] = !isPrime[n];

      n = 3*x*x + y*y;
      if( n <= limit && n % 12 == 7 )
        isPrime[n] = !isPrime[n];

      n = 3*x*x - y*y;
      if( x > y && n <= limit && n % 12 == 11 )
        isPrime[n] = !isPrime[n];

      int dummy = 0;
    }
  }

  rsUint32 k, n2;
  for(n = 5; n <= r; n++)
  {
    n2 = n*n;
    k  = n2;
    while( k <= limit )
    {
      isPrime[k] = false;
      k += n2;
    }
  }

  rsUint32 np = 2;
  for(i = 0; i <= limit; i++)
  {
    if( isPrime[i] )
    {
      printf("%d %s", i, " ");
      np++;
    }
  }
  printf("%s %d %s %d", "\n Number of primes up to", limit, ": ", np);
}

void primeSieve()
{
  static const rsUint32 N = 5000;
  rsUint32 b[N], c[N];
  rsUint32 p[1000];  // array of known primes
  rsUint32 n;

  // list all numbers from 2 onwards:
  for(n = 0; n < N; n++)
    b[n] = c[n] = n+2;

  rsUint32 bMax = b[N-1];
  rsUint32 dMax = bMax;    // maximum divisor - to be used later for the stopping criterion
  rsUint32 ib;             // index in buffer b
  rsUint32 ip = 0;         // index of current prime
  rsUint32 pi = 2;         // current prime
  rsUint32 s;              // stepsize
  rsUint32 np;             // number of known primes
  bool test;               // to verify, the algorithm is correct


  // the simple sieve of Eratosthenes:
  while( pi <= dMax )
  {
    pi = b[ip];
    s  = pi;
    ib = ip + s;
    //ib = ip*ip;
    while( ib < N )
    {
      b[ib]  = 0;
      ib    += s;
    }
    ip   = RAPT::rsArrayTools::firstIndexWithNonZeroValue(&b[ip+1], N-ip-1) + ip + 1; 
    dMax = bMax / pi;
  }
  test = true;

  // 1: fetch new prime p_k (next nonzero entry) and add it to the list of known primes, 
  // 2: strike all entries q*r where q is any of the known odd primes r is between p_{k-1}+1 and 
  //    p_k 
  // 3: go back to 1
  // init: add 2 to the list, strike all even numbers
  // 1: add 3
  // 2: q=3: strike 3*3=9 ...nah - this should also be done as initialization
  // 1: add 5
  // 2: q=3: strike 3*4=12, 3*5=15
  //    q=5: strike 5*4=20, 5*5=25
  // 1: add 7
  // 2: q=3: strike 3*6=18, 3*7=21
  //    q=5: strike 5*6=30, 5*7=35
  //    q=7: strike 7*6=42, 7*7=49
  // 1: add 11
  // 2: q= 3: strike  3*8=24,  3*9=36,  3*10= 30,  3*11= 33
  //    q= 5: strike  5*8=40,  5*9=45,  5*10= 50,  5*11= 55
  //    q= 7: strike  7*8=56,  7*9=63,  7*70= 70,  7*11= 77
  //    q=11: strike 11*8=88, 11*9=99, 11*10=110, 11*11=121
  // 1: add 13
  // 2: q= 3: strike  3*12= 36,  3*13= 39
  //    q= 5: strike  5*12= 60,  5*13= 65
  //    q= 7: strike  7*12= 84,  7*13= 91
  //    q=11: strike 11*12=132, 11*13=143
  //    q=13: strike 13*12=156, 13*13=169
  // 1: add 17
  // 2: q= 3: strike  3*14= 42,  3*15= 45,  3*16= 48,  3*17= 51
  //    q= 5: strike  5*14= 70,  5*15= 75,  5*16= 80,  5*17= 85
  //    q= 7: strike  7*14= 98,  7*15=105,  7*16=112,  7*17=119
  //    q=11: strike 11*14=154, 11*15=165, 11*16=176, 11*17=187
  //    q=13: strike 13*14=182, 13*15=195, 13*16=208, 13*17=221
  //    q=17: strike 17*14=238, 17*15=255, 17*16=272, 13*17=289
  for(ib = 2; ib < N; ib += 2)
    c[ib] = 0;
  c[7] = 0; // strike 9 manually
  ip   = 1;
  p[0] = 2;
  p[1] = 3;
  np   = 2;
  while( true )
  {
    ip    = RAPT::rsArrayTools::firstIndexWithNonZeroValue(&c[ip+1], N-ip-1) + ip + 1; 
    pi    = c[ip];
    p[np] = pi;
    np    = np+1;
    for(rsUint32 i = 1; i < np; i++)
    {
      rsUint32 q = p[i];
      for(rsUint32 r = p[np-2] + 1; r <= pi; r++)
      {
        ib = q*r-c[0];
        if( ib < N )
        c[ib] = 0;
      }
    }
    int dummy = 0;
  }


  test = RAPT::rsArrayTools::equal(b, c, N);


  // idea: when (mutiples of) 2 were already sieved out, we can use and increment of 2p for all 
  // subsequent primes p so be sieved out. when 3 has already been sieved, we can use an increment
  // 3p for all subsequent primes (verify this). when 2 and 3 already have been sieved, we can use 
  // an increment of 6p (verify this) ...but we may need multiple passes with different start-values
  // generally: let m be the smallest multiple of p that is not also a multiple of any previous 
  // prime that has been sieved, then use m-p as increment (is that true?)
  // optimized (doesn't work yet):
  ip = 0;
  pi = 2;
  while( pi <= dMax )
  {
    pi = c[ip];
    s  = pi*(pi-1);

    // we may need mutiple passes of that loop when using larger stepsizes later:
    ib = ip + s;
    while( ib < N )
    {
      c[ib]  = 0;
      ib    += s;
    }

    ib = pi * c[ip+1] + 2;

    ip   = RAPT::rsArrayTools::firstIndexWithNonZeroValue(&c[ip+1], N-ip-1) + ip + 1; 
    dMax = bMax / pi;
  }
  test = RAPT::rsArrayTools::equal(b, c, N);









  // idea:
  // 1: fetch new prime p_k (next nonzero entry) and add it to the list of known primes
  // 2: strike all entries q*r where q is any of the known odd primes r is between p_{k-1}+1 and 
  //    p_k 
  // 3: go back to 1

  // init: add 2 to the list, strike all even numbers
  // 1: add 3 to the list
  // 2: q=3: strike 2*3=6
  // 1: add 5
  // 2: q=3: strike 3*3=9, 3*5=15





  // idea:
  // 1: fetch next prime p_k
  // 2: add all nonzero numbers between p_k (inclusive) and p_k^2 (exclusive) to the list of known
  //    primes ...hmm - maybe just add the next nonzero number
  // 3: strike out p_k^2
  // 4: strike out p_k^2 + s_i where s_i = p_k * (p_{k+i+1} - p_{k+i}) where i starts at 0 and goes 
  //    up such that p_{k+i+1} takes the values of all known primes between p_{k+1} and p_k^2
  // 5: for all primes that were not already known up to this iteration, do the same for p_{k-1},
  //    p_{k-2}, etc. - i.e. pick up again the striking of mutiples of p_{k-1},... again, where the 
  //    previous iteration had to stop do to running out of known primes
  // 6: go back to 1
  // this algorithm should strike out every non-prime exactly once

  // init: add 2 to the list, strike 8
  // 1: fetch 2 
  // 2: add 3 to list 
  // 3: strike 2^2 = 4
  // 4: strike 2^2 + 2*(3-2) = 6
  // 5: nothing to do
  // 1: fetch 3
  // 2: add 5 to list
  // 3: strike 3^2 = 9
  // 4: strike 3^2 + 3*(5-3) = 15
  // 5: strike 2^2 + 2*(3-2) + 2*(5-3) = 10
  // 1: fetch 5
  // 2: add 7 to list
  // 3: strike 5^2 = 25
  // 4: strike 5^2 + 5*(7-5) = 35
  // 5: strike 3^2 + 3*(5-3) + 3*(7-5) = 21
  //    strike 2^2 + 2*(5-3) + 2*(7-5) = 12
  // 1: fetch 7
  // 2: add 11 to list
  // 3: strike 7^2 = 49
  // 4: strike 7^2 + 7*(11-7) = 77
  // 5: strike 5^2 + 5*( 7-5) + 5*(11-7) = 55
  //    strike 3^2 + 3*( 5-3) + 3*(7-5) + 3*(11-7) = 33
  //    strike 2^2 + 2*( 5-3) + 2*(7-5) + 2*(11-7) = 20
  // 1: fetch 11
  // 2: add 13 to list
  // 3: strike 11^2 = 121
  // 4: strike 11^2 + 11*(13-11) = 143
  // 5: strike 

  // hmm....14 and 16 have slipped through

  // maybe instead of striking out p_k^2 in step 3, we should strike all multiples of p_k between
  // p_k^2 and p_{k+1}^2?
  // init: add 2 to the list
  // 1: fetch 2 
  // 2: add 3 to list 
  // 3: strike 4,6,8
  // 4: ?
  // 5: nothing
  // 1: fetch 3
  // 2: add 5 to list
  // 3: strike 9,12,15,18,21,24
  // 4: 
  // 5: strike 2^2 + 2*(5-3)




  //
  // 3: strike out 3^2 = 9
  // 4: strike 3^2 + 3*(5-3) = 15
  // 5: strike 2^2 + 2*(5-3) = 
  //    strike 2^2 + 2*(3-2) =



  int dummy = 0;
}
/*
let's have a look at the prime-numbers with a prepended 1:
1,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,
131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,
263,269,271

we see that
 3 =  2 + 1
 5 =  3 + 2
 7 =  5 + 2
11 =  7 + 3 + 1
13 = 11 + 2
17 = 13 + 3 + 1
19 = 17 + 2
23 = 19 + 3 + 1
29 = 23 + 5 + 1
31 = 29 + 2
37 = 31 + 5 + 1
41 = 37 + 3 + 1
43 = 41 + 2
47 = 43 + 3 + 1
53 = 47 + 5 + 1
59 = 53 + 5 + 1
61 = 59 + 2
67 = 61 + 5 + 1
71 = 67 + 3 + 1
...can a pattern be found which elements of the already-computed array we need to add to the 
most recently computed prime to obtain the next? this would be an awesome achievement because known
methods to generate primes involve prime sieves which reuqire a lot of memory. maybe continue this 
progression and try to find a pattern. there is probably none, otherwise someone would surely
alreday have dicovered it. but who knows?
maybe write a program that writes such a progression into a nicely formatted textfile. if we
look at prime at index i, we subtract the one from index i-1 and find the numbers to add by a 
greedy-algorithm: for p[19] = 71, we compute x = p[19] - p[18] = 71 - 67 = 4, then find the 
largest element of the table that is <= 4 (which is 3), then compute x = 4 - 3 = 1, find the 
largest element in the table that is <= 1 (which is 1), then compute x = 1 - 1 = 0 -> when x == 0,
we have all the summands. maybe cast this problem in terms of prime the differences instead of the 
primes themselves (this will make the output prettier, so it might be easier to see a pattern). 
maybe, we should also look at higher order differences
*/


// Returns the last index in the ascendingly sorted array "A", where the value is less-than or
// equal-to "key", if 0 is returned, and the 0th element does not equal "key", then all values in 
// the array are either less or all are greater than key -> check this
template<class T>
rsUint32 rsBinarySearch(T A[], T key, rsUint32 imin, rsUint32 imax)
{
  while( imin < imax )
  {
    rsUint32 imid = imin/2 + imax/2; // divide before add to avoid overflow

    rsAssert(imid < imax);

    if( A[imid] < key )
      imin = imid + 1;
    else
      imax = imid;
  }
  if(A[imin] == key || imin == 0)
    return imin;
  else
    return imin-1;
}
// remove - use function from rapt rsArrayTools::splitIndex

double rsLogIntegral(double x)
{
  // very quick and dirty implementation, realising:
  // Li(x) = gamma + log(log(x)) + sum_k=1^inf (log(x))^k / (k*k!)  
  // http://mathworld.wolfram.com/LogarithmicIntegral.html
  // it tends to underestimate the true value when too little terms are taken, we should probably
  // stop not after a constant number of terms but when the delta drops below relative epsilon
  // ...or somehting

  int numTerms = 20;
  static const double gamma = 0.577215664901532860606512090082402431;
  double lx  = log(x);
  double y   = gamma + log(lx);
  double kf  = 1.0;  // k!
  double lxk = lx;   // (log(x))^k
  for(int k = 1; k <= numTerms; k++)
  {
    kf  *= k;
    y   += lxk / (k*kf);
    lxk *= lx;
  }
  return y;
}
void primeDistribution()
{
  static const rsUint32 N = 10000;
  rsUint32 p[N];
  rsFillPrimeTable(p, N);

  double x[N], y[N], ya[N];

  x[0] = 0.0; x[1] = 1.0;
  y[0] = y[1] = ya[0] = ya[1] = 0.0;
  for(rsUint32 n = 2; n < N; n++)
  {
    x[n]  = n;
    y[n]  = rsBinarySearch(p, n, (rsUint32)0, N-1) + 1;
    ya[n] = rsLogIntegral(x[n]);
  }

  plotData(N, x, y, ya);
}


rsUint32 powModular(rsUint32 base, rsUint32 exponent, rsUint32 modulus)
{
  rsUint32 result = 1;
  for(rsUint32 p = 1; p <= exponent; p++)
    result = (result * base) % modulus;
  return result;
}
rsUint32 inverseElement(rsUint32 element, rsUint32 modulus)
{
  // the inverse element is ensured to exist, when the modulus is prime. otherwise, it may or may 
  // not exist. if it doesn't, we return 0
  rsAssert(false); // not yet implemented - we need the extended euclid-algorithm for gcd
  return 0; 
}
void numberTheoreticTransform()
{
  /*
  A number-theoretic transform (NTT) is similar to a Fourier transform, but instead of working on 
  the set of complex numbers, it works on the set of (nonnegative?) integer numbers modulo M, where
  M is called the modulus. For a transform of length N, we need a primitive N-th root of unity r, 
  i.e. a number with the property r^N = 1 (primitive means: r^k != 1 for k < N). In this 
  exponentiation, the modular definition of multiplication should be used - let's denote that by 
  %*, such that x %* y := (x*y) % M. For example, with M = 11, 3 is a root of order 5 
  because 3 %* 3 %* 3 %* 3 %* 3 = 1. The primitive N-th root of unity root of unity r is called 
  an "element of order N". It has a cyclic property: r^n+k = r^k which is required for the NTT to 
  work. For the inverse transform, we also need the inverse element of r and N, which we denote by 
  ri, Ni. These are ensured to exist, if gcd(r,M) = gcd(N,M) = 1 which is satisfied trivially for 
  prime moduli (gcd(x,M) = 1 for any x if M is prime). If the modulus is prime, we also have 
  ensured that every element (except 0) has an inverse (question: element in the sense of N-th 
  root of unity or in the sense of element of the set of numbers?)
  
  For a given modulus, there do not exist N-th roots for arbitrary N, but only for some maximum 
  value maxN and it's integer divisors. For M = 11, we have maxN = 10 which has divisors 2 and 5, 
  so modulus 11 has roots of order 2, 5 and 10. Generally, if M is a prime number P, then 
  maxN = P-1. For an NTT of length N, we also want N to be some power of 2, so in base 11, the 
  only possible NTT is of length 2 which is not very useful. If M = 17, we have maxN = 17 - 1 = 16 
  which is a power of two. The integer divisors of 16 are 2, 4, 8 - roots of these orders also 
  exist, so with M = 17, we can do NTTs of length 2, 4, 8, 16 which is more useful. Generally, if 
  we choose the modulus to be a prime P = v * 2^k + 1, for some v and k, we can do all NTTs up to 
  length N = 2^k. For P = 17, we have v = 1, k = 4, N = 2^k = 16. If we choose v = 3, k = 4 we have 
  P = 3 * 2^4 + 1 = 97 which is also a prime. With M = 97, we can also do only NTTs up to 
  length 16. (Observation: for P = 17 and P = 97, all numbers 1,2,3,...,P-1 are roots of unity - is
  this a general rule for prime-moduli of the form p = v * 2^k + 1?) . 
  
  For use in NTT-based multiplication in arbitrary-length integer and arbitrary-precision floating
  point arithmetic, it is desirable, to be able to do NTTs up to as high a length as possible and
  it's also desirable to choose a base that is as high as possible, the maximum base is limited
  by the underlying interger type that we use. For rsBigInt and rsBigFloat, it seems best to use 
  base 3489660929 which allows NTTs up to length 2^28 = 268435456. 

  prime moduli suitable for NTT-multiplication in rsBigInt/Float
  15 * 2^27 + 1             2013265921 <           2147483648 = 2^31
  13 * 2^28 + 1 =           3489660929 <           4294967296 = 2^32
  87 * 2^56 + 1 =  6269010681299730433 <  9223372036854775808 = 2^63
  27 * 2^59 + 1 = 15564440312192434177 < 18446744073709551616 = 2^64


  References:

  http://www.apfloat.org/ntt.html
  http://www.math.tuwien.ac.at/~melenk/teach/numerik_WS0708/project9.pdf
  http://en.wikipedia.org/wiki/Sch%C3%B6nhage%E2%80%93Strassen_algorithm
  http://domino.mpi-inf.mpg.de/intranet/ag1/ag1publ.nsf/c1469fafbc6d09dcc12569e30040d641/ca00677497561c7ec125763c0044a41a/$FILE/gpgpu_mul.pdf


  */

  static const rsUint32 M = 97;
  rsUint32 roots[M], inverseRoots[M], orders[M], inverseOrders[M];
    // we don't expect to find M roots in general, but it's an upper bound

  RAPT::rsArrayTools::fillWithZeros(roots,  M); // so it doesn't have undefined values when we search through it
  RAPT::rsArrayTools::fillWithZeros(orders, M);

  // find the roots of unity and their orders by brute force:
  rsUint32 i, j, k = 0; 
  for(i = 2; i < M; i++) // loop over potential roots
  {
    for(j = 2; j < M; j++) // loop over potential orders
    {
      bool isRoot      = powModular(i, j, M) % M == 1;
      bool isPrimitive = !RAPT::rsArrayTools::contains(roots, M, i);
      if( isRoot && isPrimitive )
      {
        roots[k]  = i;
        orders[k] = j;
        k++;
      }
    }
  }
  rsUint32 numRoots = k;

  // find inverse elements of roots and orders by the powering algorithm:
  rsUint32 maxOrder = RAPT::rsArrayTools::maxValue(orders, M);  // all other orders divide this value
  for(k = 0; k < numRoots; k++)
  {
    inverseRoots[k]  = powModular(roots[k],  maxOrder-1, M);
    inverseOrders[k] = powModular(orders[k], maxOrder-1, M);
  }

  int dummy = 0;

  rsUint32 test = 1;
  for(i = 1; i <= 5; i++)
    test = (test * 3) % M;  
  // test should be 1, because 3^5 = 1 in modular arithmetic with M = 11

  // 4 is the inverse of 3 - check this:
  for(i = 2; i < M; i++)
  {
    test = (i    * 3) % M;
    j    = (test * 4) % M;  // j should equal i
    rsAssert(i == j);
  }
}

