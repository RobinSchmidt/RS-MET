// Experiments and functions for dealing with meshes, such as computing numerical derivatives, 
// which is an importants step in numerical PDE solvers. These functions were previously located in
// MathExperiments.cpp, but because the amount of code grew so large that at some point, i moved 
// them into this separate file. Their declarations, however, are still in MathExperiments.h 
// because  didn't want to create yet another header file and at some point i want to geht of these
// .h files anyway and instead directly include the .cpp files in the test driver files. Less 
// files, less confusion. ..oh - i think, the declarations are actually in the wrong 
// MathExperiments.h...The experiments cod is all a bit messy anyway -> clean up!


void derivativeFormulas1D()
{

  int dummy = 0;
}

void derivativeFormulas()
{
  // Tests formulas for numerical estimation of a derivative on irregularly spaced data.

  derivativeFormulas1D();
}

/** Fills edges of a graph of 2D vectors (as vertices) with a user supplied function that takes as
input the source and target vector and returns a scalar that can be used as weight for the edge 
between the two vertices. */
template<class T>
void fillEdges(rsGraph<rsVector2D<T>, T>& g, 
  const std::function<T(rsVector2D<T>, rsVector2D<T>)>& f)
{
  using Vec = rsVector2D<T>;
  for(int i = 0; i < g.getNumVertices(); i++) {
    Vec vi = g.getVertexData(i);                 // vector stored at source vertex i
    for(int j = 0; j < g.getNumEdges(i); j++) {
      int k  = g.getEdgeTarget(i, j);            // index of target vertex
      Vec vk = g.getVertexData(k);               // vector stored at target vertex k
      T ed   = f(vi, vk);                        // compute edge data via user supplied function
      g.setEdgeData(i, j, ed); }}                // ...and store it at the edge
}
// move to MeshStuff.h

void vertexMeshGradient1()
{
  // We test the function rsNumericDifferentiator::gradient2D which operates on an irregular mesh 
  // of vertices and estimates the gradient from function values known only on such a mesh.
  // ...this experiment is now somwhat obsolete because down below we have now other, more 
  // elaborate ones in place that do the same and much more - so maybe delete it at some point

  using Vec2 = rsVector2D<float>;
  using VecF = std::vector<float>;
  using VecI = std::vector<int>;
  //using Mesh = rsGraph<Vec2, rsEmptyType>;  // later use float for the edge data
  using Mesh = rsGraph<Vec2, float>;
  using ND   = rsNumericDifferentiator<float>;

  // an (irregular) star-shaped mesh with a vertex P = (3,2) at the center and 4 vertices 
  // Q,R,S,T surrounding it that are connected to it:
  Mesh mesh;
  bool sym = true;                 // select, if edges should be added symmetrically
  mesh.addVertex(Vec2(3.f, 2.f));  // P = (3,2) at index 0
  mesh.addVertex(Vec2(1.f, 3.f));  // Q = (1,3) at index 1
  mesh.addVertex(Vec2(4.f, 2.f));  // R = (4,2) at index 2
  mesh.addVertex(Vec2(2.f, 0.f));  // S = (2,0) at index 3
  mesh.addVertex(Vec2(1.f, 1.f));  // T = (1,1) at index 4
  mesh.addEdge(0, 1, sym);         // connect P to Q
  mesh.addEdge(0, 2, sym);         // connect P to R
  mesh.addEdge(0, 3, sym);         // connect P to S
  mesh.addEdge(0, 4, sym);         // connect P to T

  // Create arrays of function values and (true) partial derivatives and their numerical estimates.
  // For the estimates, only vertices with neighbors are supposed to contain a reasonable value 
  // afterwards, all others are supposed to contain zero:
  int N = mesh.getNumVertices();
  VecF u(N), u_x(N), u_y(N);     // u(x,y) and its true partial derivatives with resp. to x,y
  VecF u_x0(N), u_y0(N);         // with weighting 0 (unweighted)
  VecF u_x1(N), u_y1(N);         // with weighting 1 (sum of absolute values, "Manhattan distance")
  VecF u_x2(N), u_y2(N);         // with weighting 2 (Euclidean distance)
  VecF e_x0(N), e_y0(N);         // error of u_x0, ...
  VecF e_x1(N), e_y1(N);         // ...etc.
  VecF e_x2(N), e_y2(N);

  // Define our test function u(x,y) and its partial derivatives:
  // u(x,y)   =    sin(wx * x + px) *    sin(wy * y + py)
  // u_x(x,y) = wx*cos(wx * x + px) *    sin(wy * y + py)
  // u_y(x,y) =    sin(wx * x + px) * wy*cos(wy * y + py)
  // and a function to fill the arrays of true partial derivatives:
  float wx = 0.01f, px = 0.3f;
  float wy = 0.02f, py = 0.4f;
  auto f  = [&](float x, float y)->float { return    sin(wx * x + px) *    sin(wy * y + py); };
  auto fx = [&](float x, float y)->float { return wx*cos(wx * x + px) *    sin(wy * y + py); };
  auto fy = [&](float x, float y)->float { return    sin(wx * x + px) * wy*cos(wy * y + py); };
  auto fill = [&]() 
  { 
    int N = mesh.getNumVertices();
    u.resize(N);
    u_x.resize(N);
    u_y.resize(N);
    for(int i = 0; i < N; i++) {
      Vec2 v = mesh.getVertexData(i);
      u[i]   = f( v.x, v.y);
      u_x[i] = fx(v.x, v.y);
      u_y[i] = fy(v.x, v.y); }
  };
  // todo: later compute also 2nd derivatives u_xx, u_yy, u_xy and Laplacian u_L

  // distance functions (constant, 1/Manhattan, 1/Euclidean)
  std::function<float(Vec2, Vec2)> d0, d1, d2;
  d0 = [&](Vec2 a, Vec2 b)->float { return 1.f; };
  d1 = [&](Vec2 a, Vec2 b)->float { Vec2 d = b-a; return 1.f / (rsAbs(d.x) + rsAbs(d.y)); };
  d2 = [&](Vec2 a, Vec2 b)->float { return 1.f / rsNorm(b-a); };

  // P = (3,2), Q = (1,3), R = (4,2), S = (2,0), T = (1,1)
  fill();
  fillEdges(mesh, d0); ND::gradient2D(mesh, u, u_x0, u_y0); e_x0 = u_x-u_x0; e_y0 = u_y-u_y0;
  fillEdges(mesh, d1); ND::gradient2D(mesh, u, u_x1, u_y1); e_x1 = u_x-u_x1; e_y1 = u_y-u_y1;
  fillEdges(mesh, d2); ND::gradient2D(mesh, u, u_x2, u_y2); e_x2 = u_x-u_x2; e_y2 = u_y-u_y2;

  // This is the regular 5-point stencil that would result from using a regular mesh:
  // P = (3,2), Q = (3,3), R = (4,2), S = (3,1), T = (2,2)
  mesh.setVertexData(0, Vec2(3.f, 2.f));   // P = (3,2)
  mesh.setVertexData(1, Vec2(3.f, 3.f));   // Q = (3,3)
  mesh.setVertexData(2, Vec2(4.f, 2.f));   // R = (4,2)
  mesh.setVertexData(3, Vec2(3.f, 1.f));   // S = (3,1)
  mesh.setVertexData(4, Vec2(2.f, 2.f));   // T = (2,2)
  fill();                                  // compute target values
  fillEdges(mesh, d0); ND::gradient2D(mesh, u, u_x0, u_y0); e_x0 = u_x-u_x0; e_y0 = u_y-u_y0;
  fillEdges(mesh, d1); ND::gradient2D(mesh, u, u_x1, u_y1); e_x1 = u_x-u_x1; e_y1 = u_y-u_y1;
  fillEdges(mesh, d2); ND::gradient2D(mesh, u, u_x2, u_y2); e_x2 = u_x-u_x2; e_y2 = u_y-u_y2;
  // Of interest are mostly the errors at index 0 because that's the only vertex which has
  // more than 1 neighbour. The other vertices are boundary vertices that either have 1 
  // neighbor (if sym == true) or 0 neighbors (if sym == false). It may be interesting to figure
  // out the accuracy when there's only 1 neighbor - it should give minimum x-error and maximum 
  // y-error when the edge is horizontal and maximum x-error and minimum y-error when the edge is
  // vertical. I think, in these cases, the formula reduces to the one-sided difference in x- and 
  // y-direction. The best compromise should be obtained, when the edge is at an angle of 45 
  // degrees. For horizontale or vertical edges, the max-error may go to infinity? Division by 
  // zero? ...figure out...

  int dummy = 0;

  // Observations: 
  // -The accuracy seems to be best with using the (inverse) Manhattan distance as weights. 
  //  Why is that? Shouldn't the Euclidean distance be better? ..the values are all very similar 
  //  though, so this single experiment may not really mean much - more tests needed...
  //  -maybe try the maximum norm, too
  // -In the case of the regular grid, all estimates are the same, as they should, since all 
  //  distances are unity.
}

/** Computes the weight for the k-th neighbor of vertex i.
p is the weighting exponent for the distance and q is the exponent for the dependency 
measure. ...tbc... */
template<class T>
void setupNeighbourWeight(rsGraph<rsVector2D<T>, T>& mesh, int i, int k, T p, T q)
{
  // under construction

  using Vec2 = rsVector2D<T>;

  Vec2 vi  = mesh.getVertexData(i);
  Vec2 vk  = mesh.getVertexData(k);
  Vec2 dik = vk - vi;                 // difference vector

  T r = rsNorm(dik);                  // distance or radius
  //T w = pow(radius, -p);        // weight determined by distance
  T w = pow(r, -p);                   // weight determined by distance

  // compute and set up additional weighting factor dertmined by the (in)depence of edge direction 
  // ik from the other edge direction in, n=0,...,N-1 where N is the number of neighbors of 
  // vertex i - the formula is heuristic and experimental:
  int N = mesh.getNumEdges(i);  // number of neighbors of vertex i
  T dependency(0);              // measure, how much vi is linearly dependent on the others
  for(int n = 0; n < N; n++)
  {
    int  j  = mesh.getEdgeTarget(i, n);
    Vec2 vn = mesh.getVertexData(j);

    // maybe, include a test if n == i and if so, do not accumulate the term - experiment with both
    // variants to figure out, which one works better in practice...and/or maybe try to justify it
    // theoretically

    if(n != i)  // todo: experiment with working without this condition
    {
      T c = rsDot(vi, vn) / (rsNorm(vi)*rsNorm(vn));  // correlation between vi and vn (right?)
      dependency += rsAbs(c);  // should we use the absolute value? - i think so
      // maybe the c could also be raised to a power in addition to the sum...or maybe just one or
      // the other
    }
  }
  dependency /= N;   // does this make sense?

  // ok - we have our (heuristic) measure of dependency - use it as additional factor for the 
  // weight:
  w *= (1 - pow(dependency, q));  // this formula is also heuristic

  // Set up the weight in the mesh:
  mesh.setEdgeData(i, k, w);
  int dummy = 0;
}

template<class T>
void addPolygonalNeighbours(rsGraph<rsVector2D<T>, T>& mesh, int i,
  int numSides, T radius, T angle = T(0), T p = T(0), bool symmetric = true)
{
  int N = (int)mesh.getNumVertices();
  using Vec2 = rsVector2D<T>;
  Vec2 vi = mesh.getVertexData(i);
  for(int j = 0; j < numSides; j++)
  {
    T a = T(angle + 2*j*PI / numSides); // angle
    if(numSides == 2) a /= 2;           // special case for "2-gon": use 2 perpendiculars
    T dx = radius*cos(a);               // x-distance
    T dy = radius*sin(a);               // y-distance
    T w  = pow(radius, -p);             // edge weight
    Vec2 vj(vi.x + dx, vi.y + dy);      // position of neighbor vertex
    mesh.addVertex(vj);                 // add the new neighbour to mesh
    mesh.addEdge(i, j+N, w, symmetric); // add edge to the new neighbor and back
  }
}


// soon obsolete - this stuff is now done in rsMeshWeightCalculator2D
template<class T>
void computeEdgeWeights(rsGraph<rsVector2D<T>, T>& mesh, T p, T q = 0)
{
  using Vec2 = rsVector2D<T>;
  for(int i = 0; i < mesh.getNumVertices(); i++)
  {
    Vec2 vi = mesh.getVertexData(i);
    int  K  =  mesh.getNumEdges(i);        // number of neighbors
    for(int k = 0; k < K; k++)
    {
      int  j  = mesh.getEdgeTarget(i, k);
      Vec2 vj = mesh.getVertexData(j);
      Vec2 dv = vj - vi;                   // difference vector
      T d     = rsNorm(dv);                // distance between Vi and vj
      T w     = pow(d, -p);                // edge weight


      // new - compute sum of distances between vj and all other neighbors vn of vi:
      T s(0);
      for(int m = 0; m < K; m++)
      {
        if(m != k)
        {
          int  n  = mesh.getEdgeTarget(i, m);
          Vec2 vn = mesh.getVertexData(n);
          s += rsNorm(vj - vn);
        }
      }
      w *= pow(s, q);
    
      mesh.setEdgeData(i, k, w);
    }
  }
}
template<class T>
void assignEdgeWeights(rsGraph<rsVector2D<T>, T>& mesh, T p)
{
  int N = (int)mesh.getNumVertices();
  using Vec2 = rsVector2D<T>;
  for(int i = 0; i < N; i++) {
    Vec2 vi = mesh.getVertexData(i);
    for(int k = 0; k < mesh.getNumEdges(i); k++) {
      int  j  = mesh.getEdgeTarget(i, k);
      Vec2 vj = mesh.getVertexData(j);
      Vec2 dv = vj - vi;
      T    d  = rsNorm(dv);                  // distance between vi and vj
      T    w  = pow(d, -p);                  // edge weight
      mesh.setEdgeData(i, k, w); }}
}
// maybe the convention of using pow(d, -p) is inconvenient in a more general setting - if this 
// code is moved to the library, we should probably just use pow(d, p) and the caller should set
// the minus if he wants an inverse relationship

template<class T>
void valueAndExactDerivatives(rsGraph<rsVector2D<T>, T>& mesh, 
  std::vector<T>& u, std::vector<T>& u_x, std::vector<T>& u_y,
  std::function<T(T, T)>& f, std::function<T(T, T)>& f_x, std::function<T(T, T)>& f_y)
{
  int N = (int)mesh.getNumVertices();
  //u.resize(N); u_x.resize(N); u_y.resize(N);    // are these needed?
  for(int i = 0; i < N; i++) {
    rsVector2D<T> vi = mesh.getVertexData(i);
    u[i]   = f(  vi.x, vi.y);
    u_x[i] = f_x(vi.x, vi.y);
    u_y[i] = f_y(vi.x, vi.y); }
}
template<class T>
rsVector2D<T> gradientErrorVector(rsGraph<rsVector2D<T>, T>& mesh, int i,
  std::function<T(T, T)>& f, std::function<T(T, T)>& f_x, std::function<T(T, T)>& f_y)  
{
  using Vec = std::vector<T>;
  int N = (int)mesh.getNumVertices();
  Vec u(N), u_x(N), u_y(N), U_x(N), U_y(N);
  valueAndExactDerivatives(mesh, u, U_x, U_y, f, f_x, f_y);    // U_* is exact value
  rsNumericDifferentiator<T>::gradient2D(mesh, u, u_x, u_y);    // u_* is numeric estimate
  Vec e_x = U_x - u_x;
  Vec e_y = U_y - u_y;
  return rsVector2D<T>(e_x[i], e_y[i]);
}
template<class T>
T gradientError(rsGraph<rsVector2D<T>, T>& mesh, int i,
  std::function<T(T, T)>& f, std::function<T(T, T)>& f_x, std::function<T(T, T)>& f_y)  
{
  rsVector2D<T> e = gradientErrorVector(mesh, i, f, f_x, f_y);
  return rsMax(rsAbs(e.x), rsAbs(e.y));
}
// todo: maybe instead of computing the error at a single vertex, return the error vector (error
// at all vertices) ...maybe also output the x- and y-error separately for more detailed analysis

template<class T>
void fillMeshValues(rsGraph<rsVector2D<T>, T>& mesh, const std::function<T(T, T)>& f, 
  std::vector<T>& u)  
{
  for(int i = 0; i < mesh.getNumVertices(); i++) {
    rsVector2D<T> vi = mesh.getVertexData(i);
    u[i] = f(vi.x, vi.y); }
}
template<class T>
void fillMeshGradient(rsGraph<rsVector2D<T>, T>& mesh, 
  const std::function<T(T, T)>& f_x, 
  const std::function<T(T, T)>& f_y,
  std::vector<T>& u_x, std::vector<T>& u_y)
{
  for(int i = 0; i < mesh.getNumVertices(); i++) {
    rsVector2D<T> vi = mesh.getVertexData(i);
    u_x[i] = f_x(vi.x, vi.y);
    u_y[i] = f_y(vi.x, vi.y); }
}
template<class T>
void fillMeshHessian(rsGraph<rsVector2D<T>, T>& mesh, 
  const std::function<T(T, T)>& f_xx, 
  const std::function<T(T, T)>& f_xy, 
  const std::function<T(T, T)>& f_yy,
  std::vector<T>& u_xx, std::vector<T>& u_xy, std::vector<T>& u_yy)
{
  for(int i = 0; i < mesh.getNumVertices(); i++) {
    rsVector2D<T> vi = mesh.getVertexData(i);
    u_xx[i] = f_xx(vi.x, vi.y);
    u_xy[i] = f_xy(vi.x, vi.y);
    u_yy[i] = f_yy(vi.x, vi.y); }
}


template<class T>
void taylorExpansion2D(rsGraph<rsVector2D<T>, T>& mesh, const T* u, 
  T* u_x, T* u_y, T* u_xx, T* u_xy, T* u_yy)
{


  int dummy = 0;
}

void meshGradientErrorVsDistance()
{
  // We plot the error between the estimated partial derivatives and true partial derivatives as
  // functions of the stepsize h for various numbers of neighbors. The neighbors are arranged as
  // a regular polygon around some center vertex. As example function, we use 
  // f(x,y) = sin(x) * exp(x).

  using Real = double;
  using Vec2 = rsVector2D<Real>;
  using Vec  = std::vector<Real>;
  using Mesh = rsGraph<Vec2, Real>;
  using ND   = rsNumericDifferentiator<Real>;

  // Settings:
  int minNumSides =  2;  // minimum number of sides/neighbors (todo: try to go down to 1)
  int maxNumSides =  8;  // maximum number of sides
  int Nh          = 10;  // number of stepsizes h
  double angle    = 0.0; // rotation angle of the polygon in radians
  //Vec2 x0(0.015, 0);     // position of center vertex (slightly off the x-symmetry center)
  Vec2 x0(1, 1);     // (1,1) is nicely general - no symmetries

  // Define functions for evaluating f(x,y) and its exact partial derivatives:
  std::function<Real(Real, Real)> f, f_x, f_y;
  f   = [&](Real x, Real y)->Real { return sin(x) * exp(y); };
  f_x = [&](Real x, Real y)->Real { return cos(x) * exp(y); };
  f_y = [&](Real x, Real y)->Real { return sin(x) * exp(y); };

  // Create measurement data:
  Vec h(Nh);
  for(int i = 0; i < Nh; i++)  // Create array of stepsizes
    h[i] = pow(0.5, i);
  int numMeshes = maxNumSides-minNumSides+1;
  rsMatrix<Real> err(numMeshes, (int)h.size());
  Mesh mesh;
  GraphPlotter<Real> meshPlotter;
  for(int numSides = minNumSides; numSides <= maxNumSides; numSides++)
  {
    for(int j = 0; j < (int)h.size(); j++)
    {
      // Create mesh for a particular setting for numSides and stepsize h:
      mesh.clear();
      mesh.addVertex(x0);
      addPolygonalNeighbours(mesh, 0, numSides, h[j], angle);  // unweighted
      //if(numSides >= 3 && j == 0) meshPlotter.plotGraph2D(mesh, {0});  // plot stencil for paper

      // Compute and the record the estimation error at vertex 0:
      Real e = gradientError(mesh, 0, f, f_x, f_y);
      err(numSides-minNumSides, j) = log10(e);
    }
  }

  // We use a log-log plot: the x-axis is the (negative) power of two (we use h = ...,0.25,0.5,1.0)
  // and the y-axis is the (negative) power of 10 that gives the order of magnitude of the error:
  Vec hLog = RAPT::rsApplyFunction(h, &rsLog2);  // does not compile
  plotMatrixRows(err, &hLog[0]);

  // Numerically estimate order of the errors from two successive errors (for two successive 
  // h-values) via the formula B.36 in "Finite Difference Computing with PDEs":
  //   r_{i-1} = log(R_{i-1} / R_i) / log(h_{i-1} / h_i)
  // where r_i is the index of the experiment. That is, we want to figure out by how much the 
  // error has increased when going from one h-value to the next higher one. We make the assumption
  // that the error R as function of h follows a power rule like R(h) ~ h^r for some exponent r and 
  // this exponent is what we are interested in.
  // Our err array is already logarithmized and division of arguments translates to subtraction
  // of results. Also, we have the h-values in ascending order, so the notation from the book 
  // translates to:
  rsMatrix<Real> errorOrder(numMeshes, (int)h.size()-1);
  for(int i = 0; i < numMeshes; i++)
    for(int j = 0; j < (int)h.size()-1; j++)
      errorOrder(i,j) = (err(i,j) - err(i,j+1)) / log10(h[j]/h[j+1]);
  plotMatrixRows(errorOrder, &hLog[0]);

  // Observations:
  // -We indeed see that the slope of error increases when the number of points is increased, so 
  //  more evaluation points do indeed lead to better orders of accuracy in this numerical 
  //  experiment (can this theoretically be proved?)
  // -the errors are actually really small for 5 sides upwards! the approximation seems to be very
  //  good! is the function too easy, i.e. too smooth? maybe try a more complicated function
  // -the 2-sides solution is exactly the same as the 4-sides solution when the function is 
  //  symmetrical around the valuation point, when we choose an evaluation point slighty off from
  //  such a (local) center of symmetry, the two errors are slightly different (black and green)
  // -the error of the triangular neighborhood seems to follow a h^1 rule, quadrilateral: h^2 rule,
  //  pentagonal: h^3, hexagonal: h^4, etc. - in general: h^(n-2) where n is the number of sides
  //  of the polygon. 
  // -The (black) special case for a 2-point estimate is interesting: for our particular choice of
  //  input point x0 = (0.015,0) - it follows a h^1 rule for small h but a h^2 rule for larger h, 
  //  i.e. it seemingly gets better for larger h. I think, this an artifact arising from the fact 
  //  that our function is symmetric in the x-direction around (0,0) - and if we sit exactly on a 
  //  symmetry center, a one-sided rule (that has only h^1 order) *seems* to follow the better h^2
  //  rule only because of the symmetry of the function. Here, we are slightly off the symmetry 
  //  center, and the larger h becomes, the less important that slight offset becomes - so it looks
  //  like a h^2 rule for larger h due to the "almost symmetry" around the particular evaluation 
  //  point but is actually in general the h^1 rule. When using a different evaluation like (1,0),
  //  the 2-neighbor stencil gives indeed the expected h^1 behavior. We also fall into the h^1
  //  behavior, when the rotation angle is not zero (or in general, a multiple of 90°?) because
  //  then, the direction of none of the edges is aligned with the direction along which f is 
  //  symmetric. So, in general, 3 neighbors are no better than 2, but from 3 upwards, we get
  //  the h^(n-2) rule - for two, the rule is still h^2 ...todo: figure out what happens when
  //  n = 1...
  //  ...the reason why around an odd(!) symmetry point, a first order method appears to be 2nd
  //  order can be understood by noting the the central difference can be seen as the arithmetic 
  //  mean of the forward and backward difference and in the case of odd symmetry, these two are 
  //  the same and therfore also equal to their average - so in case of odd symmetry: 
  //  central- = forward- = backward-difference
  // -between h = 2^-5 and h = 2^-6, the h^6 rule breaks down for the octagonal neighborhood and 
  //  the function becomes erratic. This indicates that at this point we have reached the numerical
  //  precision limit and choosing even smaller h will not give any benefits anymore. This is also 
  //  confirmed by the fact, that the error is numerially of order 10^-16 which is indeed the order 
  //  of the machine epsilon for double precision. The same thing also happens for the heptagonal 
  //  neighborhood between h = 2^-6 and h = 2^-7. So, for this particular function, using 
  //  h = 2^-5 = 0.03125 with an octagonal neighborhood seems to give results that are as accurate 
  //  as it gets. When using x0 = (1.015, 0) as evaluation point, the breakdown for the octagonal
  //  neighborhood happens already between h = 2^-4 = 0.0625 and h = 2^-5, so maybe such a choice
  //  of h is already fine enough in this case.
  
  //  Octagonal neighborhoods are actually convenient to create meshes for:
  //  we just make a rectangular mesh and take the 4 direct and 4 diagonal neighbors...well...of 
  //  course...with a regular mesh, we could also just use a standard-scheme except that we also 
  //  use diagonal neighbors. However, when doing so, we should probably take into account that the
  //  diagonal neighbours are further away from the evaluation point and we should take this 
  //  account in the edge weights. In another experiment below, it was found that taking
  //  d^(-n) as weighting gave most accurate results, where d is the distance, so d^(-8) should 
  //  probably be used in this case -> more research necessarry.

  // Notes:
  // -In this experiment, edge weighting by distance makes no differences because all the edges 
  //  have actually the same length.
  // -i really should write up a paper about this...maybe titled: 
  //  "Finite Differences on Irregular Meshes based on Directional Derivatives"

  // ToDo:
  // -try what happens when we just add more points in the same directions but further out than
  //  the already existing points - will these additional points also increase the accuracy order?
  // -todo: implement, for reference, the regular forward, backward and central difference methods 
  //  and plot the error with respect to those - basically, what we want to do is a sort of 
  //  least-squares approximation to those - maybe that 2nd level of approximation is what makes 
  //  things worse
  // -try using a weighting exponent equal to the number of sides - that turned out to be a good 
  //  choice in the experiments below - hmm - it seems to make no difference here
  // -plot the log of the absolute value of error as function of the log of the neighbor distance 
  //  for neighborhoods made from regular polygons with sides n = 2,3,4,5,6,7,8,... the n = 2 case 
  //  should just be two neighbors at right angles (needs to be treated outside the loop over n)
  // -also plot the error(s) as function of some exponent p when using the p-norm as weighting 
  //  function
  // -to figure out what the best weighting is, try it with a center vertex that has a concentric
  //  polygonal neighborhood, i.e. nodes at distance 2 and 2h away (maybe the outer polygon can be 
  //  at a different angle than the inner one) - with such a neighborhood, plot the error as 
  //  function of the p in a 1/p-norm weighting for some reasonably chosen h and n=3,4,5,6,7,8
  // -it could be, that the weighting is not yet optimally chosen in the implementation but we 
  //  don't really noticed it yet because we do not have a proper test in place
  // -maybe the code should be written such that eventually the user does not need to care about
  //  assigning the edge-weights - optimal accuracy should be obtained when weights are all 1 and 
  //  the option for additional weighting is only kept in for experimentation purposes

  // Compare the results to a 2D Taylor expansion of 2nd order (i.e. a conic section) around the
  // central point:
  //
  //    u = u(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y
  //
  // Estimate the 6 coeffs from the center point and 5 neighbors and compare it to the pentagonal
  // neighborhood results. The system, we need to solve is:
  //
  //    |u0|   |1 x0 y0 x0^2 y0^2 x0*y0|   |a0|
  //    |..| = |      ..........       | * |..|
  //    |u5|   |1 x5 y5 x5^2 y5^2 x5*y5|   |a5|
  //
  // for the 6 datapoints (xi,yi,ui), i = 0,..,5 where index 0 refers to the center point at which 
  // we want to estimate the partial derivatives and indices 1..5 are the 5 neighbors. After having 
  // found Taylor coeffs, we need to evaluate the partial derivatives with respect to x and y at 
  // index 0 as:
  //
  //    u_x = a1 + 2*a3*x0 + a5*y0
  //    u_y = a2 + 2*a4*y0 + a5*x0

  Vec errorTaylor(h.size());
  for(int j = 0; j < (int)h.size(); j++)
  {
    mesh.clear();
    mesh.addVertex(x0);
    addPolygonalNeighbours(mesh, 0, 5, h[j], angle);
    //meshPlotter.plotGraph2D(mesh);
    Vec u(6), u_x(6), u_y(6);
    valueAndExactDerivatives(mesh, u, u_x, u_y, f, f_x, f_y);
    rsMatrix<double> X(6, 6), a(6, 1), z(6, 1);
    for(int i = 0; i < 6; i++) {
      Vec2 vi = mesh.getVertexData(i);
      X(i, 0) = 1;
      X(i, 1) = vi.x;
      X(i, 2) = vi.y;
      X(i, 3) = vi.x * vi.x;
      X(i, 4) = vi.y * vi.y;
      X(i, 5) = vi.x * vi.y;
      z(i, 0) = u[i];
    }
    rsLinearAlgebraNew::solve(X, a, z);
    double u_x0 = a(1, 0) + 2*a(3, 0)*x0.x + a(5, 0)*x0.y; // numerical x-derivative at center point
    double u_y0 = a(2, 0) + 2*a(4, 0)*x0.y + a(5, 0)*x0.x; // numerical y-derivative
    double e_x0 = u_x0 - u_x[0];                           // x-error: numerical minus exact
    double e_y0 = u_y0 - u_y[0];                           // y-error
    errorTaylor[j] = rsMax(rsAbs(e_x0), rsAbs(e_y0));      // error taken to be the max-abs of both

    // Verify that the Taylor polynomial indeed passes through the 6 mesh points:
    Vec taylorValues(6);
    for(int i = 0; i < 6; i++) {
      Vec2 vi = mesh.getVertexData(i);
      taylorValues[i] = a(0, 0) + a(1, 0)*vi.x + a(2, 0)*vi.y
        + a(3, 0)*vi.x*vi.x + a(4, 0)*vi.y*vi.y + a(5, 0)*vi.x*vi.y;
    }
    Vec taylorValueError = taylorValues - u;  // should be zero up to roundoff -> yep, is zero
    int dummy = 0;
  }
  Vec errorOrderTaylor(h.size()-1);
  for(int j = 0; j < (int)h.size()-1; j++)
    errorOrderTaylor[j] = log(errorTaylor[j]/errorTaylor[j+1]) / log(h[j]/h[j+1]);
  Vec errorTaylorLog = RAPT::rsApplyFunction(errorTaylor, &log10);
  //rsPlotVectorsXY(hLog, errorOrderTaylor);
  // For some reason, the plotting fails, but by inspecting the errorOrderTaylor array, we see that
  // the order is 3. That's interesting: i expected 2 because we use a 2nd order bivariate 
  // polynomial. However, 3 fits the h^(n-2) rule where n is the number of neighbors. We use 5 
  // neighbors here to find the 6 2D Taylor coeffs. So, 2D Taylor is better than using 4 neighbors
  // and worse than using 6 neighbors and exactly just as good as using 5 neighbors with the 
  // directional derivative approach. In fact, the errors using the Taylor approach are numerically
  // exactly the same (up to roundoff) as the corresponding errors for the n = 5 case using the 
  // directional derivative approach. That's a pretty nice result! :-) In conclusion, i think, the 
  // directional derivative approach is favorable over the 2D Taylor approach because it's more 
  // general (works for any number of neighbors, not just 5) and is easier to compute (no 6x6 
  // system has to be solved - only a 2x2 system - and that stays the same, regardles of n). Also
  // the directional derivative approach readily generalizes to 3D. A 2nd order 3D Taylor 
  // polynomial would need 10 coeffs for 1,x,y,z,x^2,y^2,z^2,xy,xz,yz, so we would need meshes with 
  // exactly 9 neighbors for each mesh point. That number is inconvenient for 3D meshes, just like
  // 5 is inconvenient for 2D meshes. In 3D, convenient numbers are 4 (corners of tetrahedrons), 
  // 6 (sides or corners cubes), 12 (sides and corners of cubes), etc. Also, the boundary points 
  // may not have a full set of neighbors. The directional derivative approach does not care about 
  // that. Although, the global accuracy will suffer if some vertices have a lesser number of 
  // neighbors, so maybe the boundary points should somehow have "fanning-out" edges to some of 
  // their indirect neighbors as well. Fortunately, the algo seems to be quite oblivious of the
  // question whether the edges go all generally into the same directions (which they will have to 
  // at boundary points) or having some going into opposite directions (at least with respect to 
  // the error order...i think -> verify that...).

  int dummy = 0;
}

/*
// some example functions - todo: define them as lambdas and/or std::function where needed:
double sinX_times_expY(double x, double y) { return sin(x) * exp(y); }
double sinX_plus_expY(double x, double y)  { return sin(x) + exp(y); }
double sinX_times_cosY(double x, double y) { return sin(x) * cos(y); }
// i think, they are not used anywhere
*/

void meshGradientErrorVsWeight()  // rename to vertexMeshGradientWeighting
{
  // We plot the estimation error as function of the exponent p when using a p-norm as weighting
  // for the least-squares solver, i.e. the weight is proportional to 1/d^p where d is the distance
  // between the vertex under consideration and its current neighbor.

  using Vec  = std::vector<double>;
  using Vec2 = rsVector2D<double>;

  int    numSides = 5;
  int    Np = rsMax(2*numSides, 6); // number of integer p-values
  double h = 0.25;
  double s = sqrt(2);               // scale factor for the far away nodes
  double a = PI / numSides;         // rotation angle of the outer polygon
  Vec2   x0(1, 1);                  // position of center vertex

  s  = 2;    // test - doesn't make a difference qualitatively
  a *= 1.0;  // if a = PI / numSides (as as assigned above), we see a sharp minimum at 
             // p = numSides, scaling a down shifts the minimum up and makes it less sharp, scaling
             // it up seems to have less effect -> more research needed

  int fine = 10;
  Vec p = rsLinearRangeVector(fine*Np+1, 0.0, double(Np));

  std::function<double(double, double)> f, f_x, f_y;
  f   = [&](double x, double y)->double { return sin(x) * exp(y); };
  f_x = [&](double x, double y)->double { return cos(x) * exp(y); };
  f_y = [&](double x, double y)->double { return sin(x) * exp(y); };

  rsGraph<Vec2, double> mesh;
  GraphPlotter<double> meshPlotter;
  Vec err(p.size()), errX(p.size()), errY(p.size());
  for(size_t i = 0; i < p.size(); i++)
  {
    mesh.clear();
    mesh.addVertex(x0);
    addPolygonalNeighbours(mesh, 0, numSides, h,   0.0, p[i]);
    addPolygonalNeighbours(mesh, 0, numSides, s*h, a,   p[i]);
    //err[i] = gradientEstimationError(mesh, 0, f, f_x, f_y);      // scalar error
    Vec2 ev = gradientErrorVector(mesh, 0, f, f_x, f_y); // error vector
    errX[i] = ev.x;
    errY[i] = ev.y;
    err[i]  = rsMax(rsAbs(ev.x), rsAbs(ev.y));
    //err[i]  = rsNorm(ev);
  }

  // Plot mesh and estimation error as function of p:
  meshPlotter.plotGraph2D(mesh);
  rsPlotVectorsXY(p, err, errX, errY);
  rsPlotVectorsXY(p, err);

  // plot the example function (choose one - plotting one after another doesn't work):
  GNUPlotter plt;
  //plt.plotBivariateFunction(41, -10.0, 10.0, 11,  -1.0,  1.0, &sinX_times_expY);
  //plt.plotBivariateFunction(41, -10.0, 10.0, 11,  -1.0,  1.0, &sinX_plus_expY);
  //plt.plotBivariateFunction(41, -5.0, 5.0, 41, -5.0, 5.0, &sinX_times_cosY);
  // todo: plot the function *and* the stencil, allow a std::function object to be passed to
  // GNUPlotter -> avoid referring to a global function by a function pointer
  // maybe using sin(x)*cos(y) with (x,y) = (1,1) is a nice test example
  // ..but maybe try also functions that oscillate with very different frequencies along the
  // two axes...to get good estimates, the grid neighbours with direction vectors pointing more
  // toward the axis of faster oscillation should be closer to the current vertex...but that's
  // a question of the grid generation, not of the derivative estimation formula, i think


  // Observations:
  // -The sweet spot seems to be at p == numSides. The error minimum there is actually very sharp!
  //  ...this is a very unexpected result!
  //  -Maybe when we have more information from inner neighbors available, we can more and more 
  //   afford to discard the additional information from the outer neighbors - they just don't add
  //   much useful information anymore
  // -Seems to hold only for numSides >= 3, for 2, it seems to be slightly larger, like 2.2
  // -Seems to hold only when a = PI / numSides, choosing, for example, half that value, the 
  //  minimum is further to the right and less sharp (at least for numSides = 5)
  // -If the angle is zero (an unrealistic scenario - that would be a stupid grid!), the function
  //  is monotonically decreasing, but the decrease gets quite shallow after p > numSides, so 
  //  p = numSides looks like a reasonable choice in this case, too
  // -Seems the sharp dip occurs only when a = PI / numSides, when a = 0, the function 
  //  monotonically decreases - higher values of p tend to given more and more weight to the
  //  closest neighbor - in the limit, the closest neighbor alone will make the decision, which 
  //  makes sense because it has the most accurate estimate for the directional derivative
  // -When plotting the error for x- and y-coordinate separately and together with the abs-max 
  //  of both (defined as "the" error), we see that the minimum of the error occurs near but not 
  //  exactly at the point where the two separate errors cross
  //  -Using the Euclidean norm instead af the abs-max produces a very similar error shape, so that
  //   doesn't seem to make any importanz difference.
  // -The weights may get really large - especially when h is small and p is large. Maybe we 
  //  should renormalize the weights such that the maximum weight is 1. Maybe write a function
  //  normalizeWeights that is called *once after both* calls to addPolygonalNeighbours in the loop

  // ToDo:
  // -Figure out, if this rule also holds for less regular configurations of neighbor vertices or
  //  if it's an artifact of the specific geometry chosen here
  // -Try especially those geometric configurations that commonly occur in an actual mesh
  // -Maybe try a 3-edge stencil where the edges have different lengths (like 0.5, 1, 2)
  // -Maybe try it with a lot of random lengths and angles - collect statistical evidence 
  // -Figure out what the optimal weighting is when the edges do not form a regular polygon
  // -Maybe the optimal weights should also depend on the correlations between all the neighbour 
  //  edges - if two neighbors are in the same spot, they should count as one - try it with 3 
  //  neighbours, 2 of which are in the same spot...maybe the weight of a vector should be 
  //  inversely related to the sum of its projection onto all the others
  // -Maybe use as another factor in determining the weights something like:
  //  wi *= 1 - abs(sum_k <vi,vk>)^q / N ...the term after the minus should measure, how strong
  //  the vector vi is linearly dependent on all others (maybe the summand for k=i should not enter
  //  the sum)
  //  -vi: i-th direction vector, N: number of neighbours of the current node (not node i!), 
  //   wi: weight, q: user parameter
  //  -Maybe each term should be normalized: <vi,vk>/(norm(vi)*norm(vk))
  //  -Test this formula (and maybe variations of it) with 3 edges: v1,v2 in the x- and 
  //   y-directions and a third v3, sweeping a circle - plot accuracy vs angle with and without 
  //   this formula
  //   -when v3==v1, we want weights: w1=w3=0.5, w2=1 - check, if the formula produces this
  //  -Test it also with 3 fixed neighbors and another set of 3 rotating around
  // -the way, the weights are currently used: 
  //    A.a += w * dv.x * dv.x;  // or do we need to use w^2 here?
  //    A.b += w * dv.x * dv.y;  // ...and here
  //    A.d += w * dv.y * dv.y;  // ...and here
  //    b.x += w * dv.x * du;    // ...but not here
  //    b.y += w * dv.y * du;    // ...or here
  //  they just scale the contributions to the matrix coefficients and the right-hand-side values 
  //  in the same way...and it seems to work fine - but i think, when we use error-weighting in 
  //  the original least squares formula (not the result formula), the weights would end up as 
  //  w^2 for the contributions for the matrix because of the X^T * X thing with data-matrix X
  //  ...but i tried replacing w by w^2 for the A.a += ...etc. statements and it didn't work 
  //  -> figure out, what's going on - go back to the least-squares formula:
  //  X^T * X * beta = X^T * Y 
  //  from here: https://en.wikipedia.org/wiki/Linear_least_squares#Main_formulations 
  //  what is our matrix A?...i think, it's X^T * X

  // Conclusion:
  // -To optimize the accuracy of the estimated derivaties, we should use a grid where the edges 
  //  form regular polygon and use a weighting with p = numSides.
  // -For a more irregular distribtuion of neighbors, we still need to figure out a good formula 
  //  for the weights. Such a formula should probably take into account the position of a given 
  //  neighbor with respect to all other neighbors, maybe using correlations and/or mutual 
  //  distances.

  int dummy = 0;
}

void meshGradientErrorVsAngle()
{
  // We try to find a formula for optimal weights that take into account the correlations between
  // the direction vectors. This formula could be used in addition to the weighting by the lengths.
  // To this end, we create a mesh with 3 vertices, of which 2 form an orthonormal pair and the 
  // third is swept around in a circle. we plot the estimation accuracy f as function of the angle
  // for various choices of the weight-formula exponent.

  using Vec  = std::vector<double>;
  using Vec2 = rsVector2D<double>;
  using MWC  = rsMeshWeightCalculator2D<double>;

  int formula   = 0;    // formula for the weighting
  int numAngles = 721;  // stepping the angle in 0.5 degree steps
  double h = 1./16;     // approximation stepsize
  Vec2   v0(1, 1);      // position of center vertex

  // Define example function and its partial derivatives:
  std::function<double(double, double)> f, f_x, f_y;
  f   = [&](double x, double y)->double { return sin(x) * exp(y); };
  f_x = [&](double x, double y)->double { return cos(x) * exp(y); };
  f_y = [&](double x, double y)->double { return sin(x) * exp(y); };

  rsGraph<Vec2, double> mesh;
  mesh.addVertex(v0);
  mesh.addVertex(Vec2(v0.x+h, v0.y));  // fixed vector in x-direction
  mesh.addVertex(Vec2(v0.x, v0.y+h));  // fixed vector in y-direction
  mesh.addVertex(v0);                  // this is our rotating vector
  mesh.addEdge(0, 1);
  mesh.addEdge(0, 2);
  mesh.addEdge(0, 3);
  // todo: try different angles for the pair of fixed vectors - does the angle matter? ...maybe we
  // should also make a test in which we plot the accuracy as function of angle for regular polygon
  // neighborhoods - i hope that the accuracy is angle independent...but whether or nat that's the 
  // case that may also depend on the choice of the function and the evaluation point

  // test - add a 3rd fixed vertex and edge to it at 45° angle:
  //Vec2 v = v0 + (h/sqrt(2)) * Vec2(1,1);
  //mesh.addVertex(v);
  //mesh.addEdge(0, 4);
  // todo: maybe also experiment with the 3 vectors having different lengths - we want a formula 
  // that gives accurate results even for weird meshes - and then in practice actually use good 
  // meshes

  GraphPlotter<double> meshPlotter;
  //meshPlotter.plotGraph2D(mesh);
  Vec angles = rsLinearRangeVector(numAngles, 0.0, 2*PI);
  Vec errors(numAngles);
  Vec errX(numAngles), errY(numAngles);
  for(int i = 0; i < numAngles; i++)
  {
    double a  = angles[i];
    double dx = cos(a);
    double dy = sin(a);
    mesh.setVertexData(3, Vec2(v0.x + h*dx, v0.y + h*dy));

    // factor out into a single call to: weightCalculator.calcWeights(mesh):
    MWC::initEdgeWeights(mesh);
    MWC::weightEdgesByDistances(mesh);
    MWC::weightEdgesByPositions(mesh, formula);

    Vec2 err = gradientErrorVector(mesh, 0, f, f_x, f_y);
    errX[i] = err.x;
    errY[i] = err.y;
    errors[i] = rsMax(rsAbs(err.x), rsAbs(err.y));
    //errors[i] = rsNorm(err);
    //meshPlotter.plotGraph2D(mesh);
  }

  angles = angles * (180/PI);
  rsPlotVectorsXY(angles, errors, errX, errY);
  //rsPlotVectorsXY(angles, errX, errY);
  //rsPlotVectorsXY(angles, errors);

  // Observations:
  // -Without any weighting (i.e. formula == 0), the angular dependency of the error has a somewhat
  //  sharp minimum at around 135° degrees for v0 = (1,1). However, that minimum is somewhere else
  //  for v0 = (2,2) ...but: then it's near 315° which is off by 180° from 135°, so it seems that 
  //  adding a 3rd neighbor is most effective, when it's 45° outside the 90° wedge spanned by the 
  //  two fixed neighbors.
  //  -Using the Euclidean norm instead of the minimum of x,y-errors, the maximum is less sharp
  //   and shifted towards 140°  
  // -i actually expected to see error minima whenever v3 is at at a multiple of a 45° angle and 
  //  maxima, whenever v3 coincides with v1 or v2 - but that doesn't seem to be the case
  //  -the reason for expecting this is that when v3 is equal to one of the other 2 vectors, we 
  //   have effectively only 2 evaluation points. I thought, the further away the 3rd 
  //   evaluation point is from the other two, the more additional information it gives about the
  //   function and that would make the estimate more accurate. ...but it doesn't seem so....
  // -the x-error has 3 local maxima and minima, the y-error has 2 local maxima and minima and a 
  //  saddle with h = 1/16 and v0 = (1,1)
  // -the curves look generally sine-wavei'sh like a combination of 2 sines with f and 2*f?
  // -when changing the evaluation point v0, the error curves change wildly - there does not seem
  //  to be any particular angle that minimizes the error at all possible points
  // -maybe try adding a 3rd fixed vector and see, if that changes the behavior - the step from 
  //  2 to 3 is a step from critically a determined to overdetermined system, but the step from 
  //  3 to 4 is not -> done: nope, the curves look qualitatively the same
  // -the minima are sharp, notch-like. the maxima are smooth and wide..the whole function looks
  //  a bit like piecewise rectified sines
  //  -maybe these notches are due to taking the maximum of x- and y-error - todo: plot x- and
  //   y-error separately -> done: they look smooth
  //
  // Update - using formulas based on mutual distances, etc:
  // -formula 0: baseline, reference, all weights 1
  // -formula 1: doesn't seem to give a consistent improvement

  // Conclusion:
  // -Trying to take into account the angles of the neighbours with respect to one another does not 
  //  seem to be a promising idea to reduce the estimation error. For the time being, let's focus
  //  on the lengths of the individual edges and not about their interrelations...but maybe more 
  //  research into this at some point might be a good idea.

  // ToDo:
  // -Try it with 2 fixed neighbors at 0° and 120°. Probably, the best position for the 3rd 
  //  neighbor is 240° in this case
  // -Try it with 3 fixed neighbors, letting a 4th rotate. Also, try letting 3 other neighbors 
  //  rotate.
  // -Compare the error of the mesh using the additional neighbor(s) also with a reference mesh 
  //  that doesn't use the additional neighbors. Maybe, for some positions (like, when it coincides
  //  with an existing fixed neighbor), the additional neighbor may make the error even worse by 
  //  giving too much (i.e. twice as much) weight to that neighbor (unless compensated by an 
  //  appropriate weighting function, which is what we are trying to figure out). Maybe try this 
  //  with 4 fixed and 1 rotating neighbor.

  int dummy = 0;
}

void meshGradientErrorVsIrregularity()
{
  // We create a neighborhood of a given number of vertices and investigate, how the estimation 
  // error changes, when we randomize the positions of neighbors more and more. We intially start
  // with a regular polygon with neighbors at a fixed distance, then randomize them, then randomize
  // some more, etc.

  using Real = double;
  using Vec  = std::vector<Real>;
  using Vec2 = rsVector2D<Real>;

  // Setup:
  int numSides = 8;      // number of neighbors
  int randSeed = 1;      // seed for random number generator
  int numTests = 20;     // number of tests/steps
  Real h       = 1./16;  // approximation stepsize
  Real randMin = 0.0;    // minimum randomization (as fraction of h)
  Real randMax = 0.2;    // maximum randomization (as fraction of h)
  Real p       = 4.0;    // weighting exponent for center distance (i think, numSides is best)
  Real q       = 0.0;    // weighting exponent for separation
  Vec2 v0(1, 1);         // position of center vertex


  // Define example function and its partial derivatives:
  std::function<Real(Real, Real)> f, f_x, f_y;
  f   = [&](Real x, Real y)->Real { return sin(x) * exp(y); };
  f_x = [&](Real x, Real y)->Real { return cos(x) * exp(y); };
  f_y = [&](Real x, Real y)->Real { return sin(x) * exp(y); };

  // Create measurement data:
  Vec randAmount = rsLinearRangeVector(numTests, randMin, randMax);
  Vec error(numTests);   // actually log10 of error
  rsGraph<Vec2, Real> mesh;
  GraphPlotter<Real> meshPlotter;
  for(int i = 0; i < numTests; i++)
  {
    // Create mesh for a particular setting for randomization
    mesh.clear();
    mesh.addVertex(v0);
    addPolygonalNeighbours(mesh, 0, numSides, h, 0.0);
    randomizeVertexPositions(mesh, h*randAmount[i], h*randAmount[i], 0, randSeed);
    computeEdgeWeights(mesh, p, q);
    //meshPlotter.plotGraph2D(mesh, {0});

    // Compute and the record the estimation error at vertex 0:
    Real e = gradientError(mesh, 0, f, f_x, f_y);
    error[i] = log10(e);
  }
  rsPlotVectorsXY(randAmount, error);

  // Observations:
  // -For the regular mesh, the error is of the order of e-10, but even with slight randomization
  //  it sharply increases and then it continues to increase slowly.
  // -With numSides = 8, p = 4.0 seems to be a good value. This is in contrast to earlier findings
  //  with regular geometries, where a value of p = numSides = 8 seemed to be best.
  // -Using weighting does not seem to make much of a difference:

  // ToDo:
  // -Try using other weighting functions, perhaps we have just the wrong rule. We currently use an
  //  inverse power of the distance as weighting rule: w = pow(d, -p) where d is the distance and p
  //  is a parameter (here, our "weight" variable). Maybe something else works better? Maybe try
  //  an exponential function like w = exp(-p*d). Maybe we should also take into account the 
  //  relative positions of the neighbors. When two neighbors are at (almost) the same position, 
  //  they should get less weight. In the limit, when they are at exactly the same position, they 
  //  should count as one, i.e. their weights should be divided by two.
  // -Perhaps, the weight w_k of a vertex v_k should be inversely proportional to some power of 
  //  the distance d_k from the center point v_i (closer -> more weight) but also proportional to
  //  the sum of the distances to all other neighbors (if it has more far away nieghbors, it gets 
  //  more weight)
  // -define d_k as the distance between v_i and v_k and s_k as the sum (or average) of the 
  //  distances to all other neighbors: s_k is a measure of well the the neighbor v_k is separated
  //  from all other neighbors. maybe try something like s^q / d^p
  // -The weighting function should have the following features:
   // -if all neighbors are in the same spot, they all get the same weight
  //  -in a regular arrangement, they all get the same weights
  //  -if N are in the same spot and 1 is in another spot, the 1 should get relative weight of 1
  //   and in the the cluster of N, each should get relative weight of 1/N
  //  -when 2 neighbors are on the same radial line from the center point, the outer one may 
  //   actually need a negative weight (rationale: consider the 1D function f(x) = x^2 at x=2
  //   and two point at 2.1 and 2.2: both secants will overestimate the tangent...but wait: we 
  //   don't use secants, we use a least squares algo...hmm...we'll see
  //  -maybe not only the distances between the neighbors are relevant but also their scalar 
  //   product (after subtracting off the center)
  // -Test weighting function(s) numerically with:
  //  -2 points, 1 fixed, 1 rotating
  //  -2 points, 1 fixed, 1 going outward (on the same line, or opposite, or perpendicular,...)
  //  -only randomized angles (there's actually a formula for computing mesh Laplacians involving
  //   the angles - maybe it has to with it?)
  //   -all points have the same distance, so p becomes irrelevant which is convenient to 
  //    investigate only the part of the formula that deals with the mutual neighbor positions
  //  -only randomized radii
  //  -consider also the 1D case, maybe use it as starting point
  // -Compare the gradient estimation via directional derivatives to a full 2D fit of 6 points 
  //  (2D quadratic approximation)

  int dummy = 0;
}

void vertexMeshGradient()
{
  //vertexMeshGradient1();  // somewhat obsolete now - maybe delete at some point

  meshGradientErrorVsDistance();
  meshGradientErrorVsWeight();   // todo: try with geometries other than regular polygons
  meshGradientErrorVsAngle();
  meshGradientErrorVsIrregularity();
}

template<class T>
void exactHessian(rsGraph<rsVector2D<T>, T>& mesh, 
  std::vector<T>& u_xx, std::vector<T>& u_xy,
  std::vector<T>& u_yx, std::vector<T>& u_yy,
  std::function<T(T, T)>& f_xx, std::function<T(T, T)>& f_xy,
  std::function<T(T, T)>& f_yx, std::function<T(T, T)>& f_yy)
{
  int N = (int)mesh.getNumVertices();
  for(int i = 0; i < N; i++) {
    rsVector2D<T> vi = mesh.getVertexData(i);
    u_xx[i] = f_xx(vi.x, vi.y);
    u_xy[i] = f_xy(vi.x, vi.y); 
    u_yx[i] = f_yx(vi.x, vi.y); 
    u_yy[i] = f_yy(vi.x, vi.y); }
}
template<class T>
void meshHessianViaGradGrad(rsGraph<rsVector2D<T>, T>& mesh, std::function<T(T, T)>& f,
  std::vector<T>& u_xx, std::vector<T>& u_xy,
  std::vector<T>& u_yx, std::vector<T>& u_yy)
{
  // Computes the Hessian by first computing the gradient and then computing the gradients of the
  // two partial derivatives, i.e. the two elements of the gradient vector.
  int N = (int)mesh.getNumVertices();
  std::vector<T> u(N), u_x(N), u_y(N);           // temporaries
  for(int i = 0; i < N; i++) {
    rsVector2D<T> vi = mesh.getVertexData(i);
    u[i] = f(vi.x, vi.y); }
  rsNumericDifferentiator<T>::gradientAndHessian2D(
    mesh, &u[0], &u_x[0], &u_y[0], &u_xx[0], &u_xy[0], &u_yx[0], &u_yy[0]);
}
template<class T>
rsMatrix2x2<T> hessianErrorMatrix(rsGraph<rsVector2D<T>, T>& mesh, int i,
  std::function<T(T, T)>& f,
  std::function<T(T, T)>& f_xx, std::function<T(T, T)>& f_xy,
  std::function<T(T, T)>& f_yx, std::function<T(T, T)>& f_yy)
{
  using Vec = std::vector<T>;
  int N = (int)mesh.getNumVertices();
  Vec U_xx(N), U_xy(N), U_yx(N), U_yy(N);                               // exact values
  exactHessian(mesh, U_xx, U_xy, U_yx, U_yy, f_xx, f_xy, f_yx, f_yy);
  Vec u_xx(N), u_xy(N), u_yx(N), u_yy(N);                               // numeric estimates
  meshHessianViaGradGrad(mesh, f, u_xx, u_xy, u_yx, u_yy);
  Vec e_xx = U_xx - u_xx; Vec e_xy = U_xy - u_xy;                       // errors
  Vec e_yx = U_yx - u_yx; Vec e_yy = U_yy - u_yy;
  rsMatrix2x2<T> E;
  E.a = e_xx[i]; E.b = e_xy[i];
  E.c = e_yx[i]; E.d = e_yy[i];
  return E;
}
template<class T>
T hessianError(rsGraph<rsVector2D<T>, T>& mesh, int i,
  std::function<T(T, T)>& f, 
  std::function<T(T, T)>& f_xx, std::function<T(T, T)>& f_xy,
  std::function<T(T, T)>& f_yx, std::function<T(T, T)>& f_yy)
{
  rsMatrix2x2<T> E = hessianErrorMatrix(mesh, i, f, f_xx, f_xy, f_yx, f_yy);
  return rsMax(rsAbs(E.a), rsAbs(E.b), rsAbs(E.c), rsAbs(E.d));
}
template<class T>
void createMeshForHessianEstimation(rsGraph<rsVector2D<T>, T>& mesh, int numSides, T h, 
  rsVector2D<T> x0)
{
  // Creates a regular polygonal neighborhood around x0 and for all those created neighbors, it 
  // creates such a neighborhood as well.
  mesh.clear();
  mesh.addVertex(x0);
  addPolygonalNeighbours(mesh, 0, numSides, h, T(0), T(0), false);  // unweighted
  for(int k = 1; k <= numSides; k++)
    addPolygonalNeighbours(mesh, k, numSides, h, T(0), T(0), false);
}

void meshHessianErrorVsDistance()
{
  // We use the rsNumericDifferentiator<T>::gradient2D function on the two gradients again to 
  // estimate the Hessian matrix and measure the error in doing so. It's like 
  // meshGradientErrorVsDistance but for the Hessian instead of the gradient.

  using Vec2 = rsVector2D<double>;
  using Vec  = std::vector<double>;
  using ND   = rsNumericDifferentiator<double>;

  // Settings:
  int minNumSides =  2;  // minimum number of sides/neighbors (todo: try to go down to 1)
  int maxNumSides =  8;  // maximum number of sides
  int Nh          = 10;  // number of stepsizes h
  Vec2 x0(1, 1);         // (1,1) is nicely general - no symmetries

  // Define example function and its 2nd order partial derivatives:
  std::function<double(double, double)> f, f_xx, f_xy, f_yy;
  f    = [&](double x, double y)->double { return  sin(x) * exp(y); };
  f_xx = [&](double x, double y)->double { return -sin(x) * exp(y); };
  f_xy = [&](double x, double y)->double { return  cos(x) * exp(y); };
  f_yy = [&](double x, double y)->double { return  sin(x) * exp(y); };

  // Create measurement data:
  Vec h(Nh);
  for(int i = 0; i < Nh; i++)  // Create array of stepsizes
    h[i] = pow(0.5, i);
  int numMeshes = maxNumSides-minNumSides+1;
  rsMatrix<double> err(numMeshes, (int)h.size());
  rsGraph<Vec2, double> mesh;
  GraphPlotter<double> meshPlotter;
  for(int numSides = minNumSides; numSides <= maxNumSides; numSides++)
  {
    for(int j = 0; j < (int)h.size(); j++)
    {
      createMeshForHessianEstimation(mesh, numSides, h[j], x0);
      //if(j == 0) meshPlotter.plotGraph2D(mesh, {0});

      // Compute the and record the estimation error at vertex 0:
      double e = hessianError(mesh, 0, f, f_xx, f_xy, f_xy, f_yy);
      err(numSides-minNumSides, j) = log10(e);
      // todo: maybe recor errors for u_xx, u_xy, etc. separately
    }
  }

  // Plot:
  Vec hLog(h.size()); for(size_t i = 0; i < h.size(); i++) hLog[i] = rsLog2(h[i]);
  plotMatrixRows(err, &hLog[0]);


  rsMatrix<double> errorOrder(numMeshes, (int)h.size()-1);
  for(int i = 0; i < numMeshes; i++)
    for(int j = 0; j < (int)h.size()-1; j++)
      errorOrder(i,j) = (err(i,j) - err(i,j+1)) / log10(h[j]/h[j+1]);
  plotMatrixRows(errorOrder, &hLog[0]);


  int dummy = 0;

  // Observations:
  // -We see similar plots as in meshGradientErrorVsDistance, as expected.
  // -The error in the mixed derivatives u_xy, u_yx seems to be orders of magnitude less than the 
  //  error in u_xx, u_yy (this can't be seen in the plot because there, we plot only the maximum 
  //  of the 4 errors, but it can be seen in the debugger inspecting the E matrix in hessianError)

  // ToDo:
  // -Compute Laplacian and compare results to simplifed algorithm using the weighted average of 
  //  the neighborhood.

  // Notes:
  // -It's really important that we call addPolygonalNeighbours with false for the "symmetric" 
  //  parameter - otherwise, the results make no sense at all.

  // ..in the PDE solver context: can we take a DSP prespective and see each frame as a filtered
  // version of the previous frame with some suitable kernel? if so, what can the kernel tell us 
  // about numerical dispersion and diffusion? is this related to so specific way of accumulating
  // the local errors (some sort of interference maybe?)

  // Other ideas to estimate the Hessian:
  // The function meshHessian computes an estimate of the Hessian matrix by computing gradients of
  // gradients, using the same code for two levels of derivatives. Maybe it could be more efficient
  // and/or accurate to set up equations like:
  //   u_a ~= <g,a> + <a,H,a> = u_x*a_x + u_y*a_y + a_x^2*u_xx + 2*a_x*a_y*u_xy + a_y^2*u_yy
  //   u_b ~= <g,b> + <b,H,b> = u_x*b_x + u_y*b_y + b_x^2*u_xx + 2*b_x*b_y*u_xy + b_y^2*u_yy
  //   ...
  // to jointly estimate Hessian and gradient. That would require each vertex to have at least 5 
  // neighbors in order to not have an underdetermined 5x5 system. This is not required with the 
  // 2-level appraoch, because it indirectly incorporates information form 2nd order neighbors for 
  // computing the elements of the Hessian. Maybe we could also set up equations like:
  //   u_a ~= <g,a> + <a,H,a>  where the latter means the a^T * H * a
  //   u_a - <g,a> = <a,H,a>
  // where g is already known (computed as usual) and for determining the 3 independent values of
  // H, we would solve a 3x3 system. I really don't know which is best. Maybe it's not only about 
  // accuracy or efficiency but also about numerical stability, dispersion and diffusion when used
  // in a PDE solver. Maybe several possibilities should be checked out
}




template<class T>
void solveSymmetric3x3(
  T a11, T a12, T a13, T a22, T a23, T a33, T b1, T b2, T b3, T* x1, T* x2, T* x3)
{
  T a12_2 = a12*a12;
  T a13_2 = a13*a13;
  T a23_2 = a23*a23;
  T D = (a13_2*a22 - 2*a12*a13*a23 + a12_2*a33 + (a23_2 - a22*a33)*a11);
  T s = T(1) / D;
  *x1 =  s*((a33*b2 - a23*b3)*a12 - (a23*b2 - a22*b3)*a13 + (a23_2   - a22*a33)*b1);
  *x2 =  s*(a13_2*b2 - a12*a13*b3 - (a33*b2 - a23*b3)*a11 - (a13*a23 - a12*a33)*b1);
  *x3 = -s*(a12*a13*b2 - a12_2*b3 - (a23*b2 - a22*b3)*a11 - (a13*a22 - a12*a23)*b1);
}
// move to library...hmm..i this really any better than just doing Gaussian elimination?
// -> benchmark!

// Computes Hessian, when the gradient is already known - needs more testing - test it with the exact
// gradient, mabye using a quadratic function, i.e. f(x,y) = A + B*x + C*y + D*x^2 + E*y^2 + F*x*y
template<class T>
void hessian2DViaTaylor(const rsGraph<rsVector2D<T>, T>& mesh, 
  const T* u, const T* u_x, const T* u_y, int i,
  T* u_xx, T* u_xy, T* u_yy)
{
  using Vec2       = rsVector2D<T>;           // shorthand for convenience
  const Vec2& vi   = mesh.getVertexData(i);   // vertex i, at which we calculate the Hessian
  int numNeighbors = mesh.getNumEdges(i); 

  if(numNeighbors < 3) { *u_xx = *u_xy = *u_yy = T(0); return; } // preliminary

  T a11(0), a12(0), a13(0), a22(0), a23(0), a33(0);  // matrix elements
  T b1(0), b2(0), b3(0);                             // vector elements of right hand side
  for(int k = 0; k < numNeighbors; k++)              // loop over neighbors of vertex i
  {
    // Retrieve or compute intermediate variables:
    int j = mesh.getEdgeTarget(i, k);         // index of current neighbor of vi
    const Vec2& vj = mesh.getVertexData(j);   // current neighbor of vi
    Vec2 dv = vj - vi;                        // difference vector
    T vx  = dv.x,  vy  = dv.y;
    T vx2 = vx*vx, vy2 = vy*vy;

    // Accumulate matrix elements:
    a11 += vx2*vx2;     // vx^4
    a12 += vx*vx2*vy;   // vx^3 * vy
    a13 += vx2*vy2;     // vx^2 * vy^2
    a22 += vx2*vy2;     // vx^2 * vy^2 ..still the same as a13, scaled later (can be optimized)
    a23 += vy*vy2*vx;   // vy^3 * vx
    a33 += vy2*vy2;     // vy^4

    // Accumulate vector elements for right hand side:
    T q = T(2)*(u[j] - u[i] - (u_x[i]*vx + u_y[i]*vy));
    b1 += q*vx2;
    b2 += q*vx*vy;
    b3 += q*vy2;
  }
  // todo: maybe add weighting

  // Some elements must be scaled:
  a12 *= T(2);
  a22 *= T(4);
  a23 *= T(2);
  b2  *= T(2);

  // Solve the symmetric 3x3 system A*x = b (optimize later)
  rsMatrix<T> A(3,3), x(3,1), b(3,1);
  A(0,0) = a11;
  A(1,1) = a22;
  A(2,2) = a33;
  A(0,1) = A(1,0) = a12;
  A(0,2) = A(2,0) = a13;
  A(1,2) = A(2,1) = a23;
  b(0,0) = b1;
  b(1,0) = b2;
  b(2,0) = b3;
  rsLinearAlgebraNew::solve(A, x, b);
  *u_xx = x(0,0);
  *u_xy = x(1,0);
  *u_yy = x(2,0);

  // optimized:
  //solveSymmetric3x3(a11, a12, a13, a22, a23, a33, b1, b2, b3, u_xx, u_xy, u_yy);

  int dummy = 0;
}
template<class T>
void hessian2DViaTaylor(const rsGraph<rsVector2D<T>, T>& mesh, 
  const T* u, const T* u_x, const T* u_y, 
  T* u_xx, T* u_xy, T* u_yy)
{
  for(int i = 0; i < mesh.getNumVertices(); i++)
    hessian2DViaTaylor(mesh, u, u_x, u_y, i, &u_xx[i], &u_xy[i], &u_yy[i]);
}

void meshHessianViaTaylorErrorVsDistance()
{
  // Tests hessian2DViaTaylor in two ways: using the exact gradient as input and using the numeric
  // estimation of the gradient

  using Vec2 = rsVector2D<double>;
  using Vec  = std::vector<double>;
  using ND   = rsNumericDifferentiator<double>;

  // Settings:
  int minNumSides =  3;  // minimum number of sides/neighbors (todo: try to go down to 2 or 1)
  int maxNumSides =  9;  // maximum number of sides
  int Nh          = 10;  // number of stepsizes h
  Vec2 x0(1, 1);         // (1,1) is nicely general - no symmetries

  std::function<double(double, double)> f, f_x, f_y, f_xx, f_xy, f_yy;
  double a = 0.75;
  double b = 0.5;
  f    = [&](double x, double y)->double { return      sin(a*x) *     exp(b*y); };
  f_x  = [&](double x, double y)->double { return    a*cos(a*x) *     exp(b*y); };
  f_y  = [&](double x, double y)->double { return      sin(a*x) *   b*exp(b*y); };
  f_xx = [&](double x, double y)->double { return -a*a*sin(a*x) *     exp(b*y); };
  f_xy = [&](double x, double y)->double { return    a*cos(a*x) *   b*exp(b*y); };
  f_yy = [&](double x, double y)->double { return      sin(a*x) * b*b*exp(b*y); };


  // Create measurement data:
  Vec h(Nh);
  for(int i = 0; i < Nh; i++)  // Create array of stepsizes
    h[i] = pow(0.5, i);
  int numMeshes = maxNumSides-minNumSides+1;
  rsMatrix<double> err1(numMeshes, (int)h.size());  // error with exact gradients
  rsMatrix<double> err2(numMeshes, (int)h.size());  // error with numeric gradients
  rsGraph<Vec2, double> mesh;
  GraphPlotter<double> meshPlotter;
  int N = maxNumSides + 1;
  Vec u(N);                                       // mesh function values
  Vec U_x(N), U_y(N), U_xx(N), U_xy(N), U_yy(N);  // gradient and Hessian (exact)
  Vec u_x(N), u_y(N), u_xx(N), u_xy(N), u_yy(N);  // gradient and Hessian estimates
  for(int numSides = minNumSides; numSides <= maxNumSides; numSides++)
  {
    for(int j = 0; j < (int)h.size(); j++)
    {
      // Create mesh for a particular setting for numSides and stepsize h:
      mesh.clear();
      mesh.addVertex(x0);
      addPolygonalNeighbours(mesh, 0, numSides, h[j], 0.0);  // unweighted
      //if(j == 0) meshPlotter.plotGraph2D(mesh, {0});

      fillMeshValues(  mesh, f, u);
      fillMeshGradient(mesh, f_x, f_y, U_x, U_y);
      fillMeshHessian( mesh, f_xx, f_xy, f_yy, U_xx, U_xy, U_yy);

      // Compute error for estimating the Hessian with an exact gradient:
      hessian2DViaTaylor(mesh, &u[0], &U_x[0], &U_y[0], &u_xx[0], &u_xy[0], &u_yy[0]);
      double e_xx = U_xx[0] - u_xx[0];
      double e_xy = U_xy[0] - u_xy[0];
      double e_yy = U_yy[0] - u_yy[0];
      //double eMin = rsMin(rsAbs(e_xx), rsAbs(e_xy), rsAbs(e_yy));
      double eMax = rsMax(rsAbs(e_xx), rsAbs(e_xy), rsAbs(e_yy));
      err1(numSides-minNumSides, j) = log10(eMax);
      //err1(numSides-minNumSides, j) = log10(eMin);

      // Compute error for estimating the Hessian with a numeric gradient:
      ND::gradient2D(mesh, &u[0], &u_x[0], &u_y[0]); 
      hessian2DViaTaylor(mesh, &u[0], &u_x[0], &u_y[0], &u_xx[0], &u_xy[0], &u_yy[0]);
      e_xx = U_xx[0] - u_xx[0];
      e_xy = U_xy[0] - u_xy[0];
      e_yy = U_yy[0] - u_yy[0];
      eMax = rsMax(rsAbs(e_xx), rsAbs(e_xy), rsAbs(e_yy));
      err2(numSides-minNumSides, j) = log10(eMax);

      int dummy = 0;
    }
  }

  // ToDo: 
  // -plot errors as funtion of log2 of h
  // -compute and plot estimated error order

  // Plot:
  Vec hLog(h.size()); for(size_t i = 0; i < h.size(); i++) hLog[i] = rsLog2(h[i]);
  plotMatrixRows(err1, &hLog[0]);
  //plotMatrixRows(err2, &hLog[0]);


  // Observations:
  // -For numSides >= 5, it doesn't seem to make a difference, if we use the exact or numeric
  //  gradient.
  // -numSides = 6 is better than numSides = 5 and numSides = 7 is still better, but not in terms 
  //  of order, just in terms of a constant factor. Going higher than 7 gives no additional 
  //  improvement. This is unexpected! Using more information does not lead to further improvement.
  //  Why is this so? Could it be that the estimation of the Hessian has an uppr limit for the
  //  error orde when using only information from direct neighbors and to get further improvements
  //  one has to take indirect neighbors into account?
  // -numSides = 3 seems to have the same order as 5, but 4 seems to be of order 0, i.e. the error
  //  does not decrease at all with decreasing h. WTF? With 4, only a11 and a33 are nonzero - is 
  //  that the problem? could the linear solver have problems with a system like that? -> implement
  //  optimized symmetric 3x3 solver and try with that - done - it seems to produce NaNs or infs
  //  with 4 sides and the behavior for sides > 7 looks the same
  //  ...hmm - maybe the u_xy derivative can't be estimated when the neighbors are purely in the x-
  //  or y direction? how would we try to estimate u_xy "by hand" on such a grid, maybe like this:
  //  -estimate 1st derivatives by central, forward or backward diff at center, north, east, 
  //   south west
  // -from these 4 estimates of the 1st derivatives, ...
  // -...hmm...i really begin to think, that estimating the mixed 2nd derivative is problematic 
  //  using only direct neighbors - at least with sides = 4
  // -double check the derivation - maybe i've missed some term

  // ToDo:
  // -Compare performance with taking gradients of gradients. This seems to be better behaved with
  //  respect to error reduction. But is it more or less costly. It takes more memory, but what
  //  about CPU cycles?

  int dummy = 0;
}




void testHessianRectangularMesh()
{
  // Creates a somewhat more realistic mesh and tests the computation of the Hessian on it. As 
  // function, we use quadratic form: f(x,y) = A + B*x + C*y + D*x^2 + E*y^2 + F*x*y
  // i think, it should be possible to etsimate the gradient and Hessian exactly for this function
  // (-> verify!)

  // Setup:
  using Real = double;
  int densityX = 11;                    // number of x-samples
  int densityY = 11;                    // number of y-samples
  Real randX   = 0.2 / densityX;        // randomization amount of x-coordinates
  Real randY   = 0.2 / densityY;        // randomization amount of y-coordinates
  Real A = 1.f;
  Real B = 2.f;
  Real C = 3.f;
  Real D = 4.f;
  Real E = 5.f;
  Real F = 6.f;

  // Create the exact functions to compute u(x,y) and its various derivatives:
  std::function<Real(Real, Real)> f, f_x, f_y, f_xx, f_xy, f_yy;
  f    = [&](Real x, Real y)->Real { return A + B*x + C*y +   D*x*x +   E*y*y + F*x*y; };
  f_x  = [&](Real x, Real y)->Real { return     B         + 2*D*x             + F*y  ; };
  f_y  = [&](Real x, Real y)->Real { return           C             + 2*E*y   + F*x  ; };
  f_xx = [&](Real x, Real y)->Real { return                   2*D                    ; };
  f_xy = [&](Real x, Real y)->Real { return                                   + F    ; };
  f_yy = [&](Real x, Real y)->Real { return                           2*E            ; };

  // For convenience:
  using Vec  = std::vector<Real>;
  using Vec2 = rsVector2D<Real>;
  using AT   = rsArrayTools;
  using ND   = rsNumericDifferentiator<Real>;
  using MWC  = rsMeshWeightCalculator2D<Real>;

  // Create the mesh:
  rsMeshGenerator2D<Real> meshGen;
  meshGen.setNumSamples(densityX, densityY);
  meshGen.setTopology(rsMeshGenerator2D<Real>::Topology::plane);
  meshGen.setParameterRange(0.f, 1.f, 0.f, 1.f);             // rename to setRange
  meshGen.updateMeshes();                                    // get rid of this
  rsGraph<Vec2, Real> mesh = meshGen.getParameterMesh();     // rename mesh to graphMesh, getP.. to getMesh
  randomizeVertexPositions(mesh, randX, randY);            // looks wrong
  // todo: assign edge weights


  // Compute function values and exact derivtaives on the mesh. Maybe wrap these into a class and 
  // then provide functions to generate various mesh functions, such here, we can just do:
  // MeshData md = getMeshDataQuadraticForm(); ..SinCos, SinExp, etc such that we can conveniently
  // switch back and forth between various functions:
  int N = mesh.getNumVertices();
  Vec u(N);                                            // Mesh function values
  Vec u_x(N), u_y(N);                                  // Gradient (numerical)
  Vec u_xx(N), u_xy(N), u_yx(N), u_yy(N);              // Hessian (numerical)
  Vec U_x(N), U_y(N), U_xx(N), U_xy(N), U_yy(N), L(N); // Gradient, Hessian and Laplacian (exact)
  for(int i = 0; i < N; i++)
  {
    Vec2 v  = mesh.getVertexData(i);
    u   [i] = f(   v.x, v.y);
    U_x [i] = f_x( v.x, v.y);
    U_y [i] = f_y( v.x, v.y);
    U_xx[i] = f_xx(v.x, v.y);
    U_xy[i] = f_xy(v.x, v.y);
    U_yy[i] = f_yy(v.x, v.y);
    L   [i] = U_xx[i] + U_yy[i];
  }

  // Plot the mesh:
  GraphPlotter<Real> plt;
  plt.plotGraph2D(mesh);

  // Compute Hessian numerically via Taylor using the exact gradient:
  hessian2DViaTaylor(mesh, &u[0], &U_x[0], &U_y[0], &u_xx[0], &u_xy[0], &u_yy[0]);
  rsPlotVectors(u_xx-U_xx, u_xy-U_xy, u_yy-U_yy);
  // Looks generally good, but some values are wrong. Maybe these are the boundary values?

  // Compute Hessian using a direct 2nd order 2D Taylor expansion without estimating the gradient 
  // first - this should also give exact results:
  taylorExpansion2D(mesh, &u[0], &u_x[0], &u_y[0], &u_xx[0], &u_xy[0], &u_yy[0]);
  rsPlotVectors(u_x-U_x, u_y-U_y);
  rsPlotVectors(u_xx-U_xx, u_xy-U_xy, u_yy-U_yy);

  // Now the same thing but with a numeric gradient:
  //MWC::weightEdgesByDistances(mesh);     // test
  //MWC::weightEdgesByPositions(mesh, 1);
  ND::gradient2D(mesh, u, u_x, u_y);
  hessian2DViaTaylor(mesh, &u[0], &u_x[0], &u_y[0], &u_xx[0], &u_xy[0], &u_yy[0]);
  rsPlotVectors(u_x-U_x, u_y-U_y);
  rsPlotVectors(u_xx-U_xx, u_xy-U_xy, u_yx-U_xy, u_yy-U_yy);
  // -Error is very high, already for the gradient estimate -> seems like something is wrong
  //  -> try using a different edge weighting function -> uncommenting the MWC:: stuff doesn't seem
  //     to help
  //  -> try a direct Taylor expansion
  // -Seems to have also problems at the boundaries

  // Now using the gradient-of-gradient algorithm:
  meshHessianViaGradGrad(mesh, f, u_xx, u_xy, u_yx, u_yy);    // todo: pass u instead of f
  rsPlotVectors(u_xx-U_xx, u_xy-U_xy, u_yx-U_xy, u_yy-U_yy);
  // -Error is even higher
  // -Has also problems at the boundaries




  // ToDo: 
  // -randomize the vertex positions and check, if it still works

  int dummy = 0;
}

void vertexMeshHessian()
{
  //meshHessianErrorVsDistance();
  //meshHessianViaTaylorErrorVsDistance();
  testHessianRectangularMesh();
}

// First computes gradient, then Hessian, then Laplacian
template<class T>
void laplacian2DViaTaylor(const rsGraph<rsVector2D<T>, T>& mesh,
  const std::vector<T>& u, std::vector<T>& L)
{
  using Vec = std::vector<T>;
  using ND  = rsNumericDifferentiator<T>;
  int N = mesh.getNumVertices();
  Vec u_x(N), u_y(N), u_xx(N), u_xy(N), u_yy(N);
  ND::gradient2D(    mesh, &u[0], &u_x[0], &u_y[0]);
  hessian2DViaTaylor(mesh, &u[0], &u_x[0], &u_y[0], &u_xx[0], &u_xy[0], &u_yy[0]);
  for(int i = 0; i < N; i++)
    L[i] = u_xx[i] + u_yy[i];

  // ToDo: move to rsNumericDifferentiator, implement a version that generates a sparse matrix A to
  // compute L = A*u
  // ...hmm - produces large error - ToDo: test hessian2DViaTaylor
}

// still tyring to figure out the right formula:
template<class T>
void laplacian2D_1(const rsGraph<rsVector2D<T>, T>& mesh, 
  const std::vector<T>& u, std::vector<T>& L)
{
  using Vec2 = rsVector2D<T>;
  int N = mesh.getNumVertices();
  rsAssert((int) u.size() == N);
  rsAssert((int) L.size() == N);
  for(int i = 0; i < N; i++) {                     // loop over all vertices
    Vec2 vi = mesh.getVertexData(i);               // current vertex location
    T uSum = T(0);                                 // weighted sum of neighbors
    T wSum = T(0);                                 // sum of weights
    T dSum = T(0);                                 // sum of squared distances
    T dMax = T(0);
    T dMin = T(100000);
    for(int k = 0; k < mesh.getNumEdges(i); k++) { // loop over vi's neighbors
      int j = mesh.getEdgeTarget(i, k);            // index of current neighbor
      T w   = mesh.getEdgeData(  i, k);            // weight in weighted sum of neighbors
      Vec2 vj = mesh.getVertexData(j);             // location of current neighbor
      Vec2 dv = vj - vi;                           // difference vector
      T d  = dv.x*dv.x + dv.y*dv.y;                // squared distance between vi and vj

      //uSum += w*(u[j]-u[i])/d;                     // accumulate sum of ...
      uSum += w*(u[j]-u[i]);                       // accumulate sum of ...
      //uSum += w*(u[j]-u[i])/d;

      wSum += w;                                   // accumulate sum of weights
      //wSum += w/d;                                   // accumulate sum of weights

      //dSum += d;
      dSum += sqrt(d);
      dMax  = rsMax(dMax, d);
      dMin  = rsMin(dMin, d);
    }
    T dAvg = dSum / mesh.getNumEdges(i);           // average of squared distances
    dAvg *= dAvg;


    //L[i] = T(4)*uSum/wSum;
    L[i] = T(4)*uSum/(wSum*dAvg);
    //L[i] = T(4)*uSum/(wSum*dMax);
    //L[i] = T(4)*uSum/(wSum*dMin);
  }
}
// wait: this actually does not compute the difference of u[i] to the weighted mean of its local 
// neighborhood - instead, it computes a weighted mean of the differences of u[i] from its 
// neighbors - the difference is taken inside the loop rather than outside. does that make a 
// difference?

template<class T>
void laplacian2D_2(const rsGraph<rsVector2D<T>, T>& mesh, 
  const std::vector<T>& u, std::vector<T>& L)
{
  using Vec2 = rsVector2D<T>;
  int N = mesh.getNumVertices();
  rsAssert((int) u.size() == N);
  rsAssert((int) L.size() == N);
  for(int i = 0; i < N; i++) {                     // loop over all vertices
    Vec2 vi = mesh.getVertexData(i);               // current vertex location
    T uSum = T(0);                                 // weighted sum of neighbors
    T wSum = T(0);                                 // sum of weights
    T dMax = T(0);
    for(int k = 0; k < mesh.getNumEdges(i); k++) { // loop over vi's neighbors
      int j = mesh.getEdgeTarget(i, k);            // index of current neighbor
      T w   = mesh.getEdgeData(  i, k);            // weight in weighted sum of neighbors
      Vec2 vj = mesh.getVertexData(j);             // location of current neighbor
      Vec2 dv = vj - vi;                           // difference vector
      T d2 = dv.x*dv.x + dv.y*dv.y;                // squared distance between vi and vj
      w    /= d2;
      uSum += w*u[j];                              // accumulate weighted sum of neighbors
      wSum += w;                                   // accumulate sum of weights
      dMax  = rsMax(dMax, d2);
    }

    T mean = uSum/wSum;
    L[i]   = mean - u[i];

    L[i] *= T(12)/dMax;
    // found empricially, seems to work for numSides = 6
  }
}
// The basic idea is that the Laplacian measures, how far the value u[i] is away from its local 
// neighborhood, so we compute a weighted mean of this neighborhood and the difference to the 
// actual value u[i]. But this must be scaled by the (maximum? minimu? average?) distance to the 
// neighbors
// see:
// https://en.wikipedia.org/wiki/Discrete_Laplace_operator#Mesh_Laplacians
// https://people.eecs.berkeley.edu/~jrs/meshpapers/Sorkine.pdf
// https://cybertron.cg.tu-berlin.de/~philipp/SGP2015/files/HerholzSGP2015.pdf
// http://ddg.math.uni-goettingen.de/pub/Polygonal_Laplace.pdf

// unrelated, but could be useful for mesh-generation:
// https://en.wikipedia.org/wiki/Laplacian_smoothing


template<class T>
void laplacian2D_3(const rsGraph<rsVector2D<T>, T>& mesh,
  const std::vector<T>& u, std::vector<T>& L)
{
  // 3rd attempt, idea:
  // -We first estimate the gradient at the center point using rsNumericDifferentiator::gradient2D.
  // -From this gradient, we compute the directional derivatives into the directions of the 
  //  neighbors (via a dot product).
  // -From the neighbor's value together with the directional derivative, we try to predict the
  //  center value (or maybe the other way around: try to predict the neighbor value from the 
  //  center value and directional derivative).
  // -The error between prediction and actual value is computed.
  // -The (weighted) average over all neighbors of this prediction error is our estimate for the 
  //  Laplacian.

  rsFill(L, T(0));  // preliminary

  using ND   = rsNumericDifferentiator<T>;

  using Vec2 = rsVector2D<T>;
  int N = mesh.getNumVertices();
  rsAssert((int) u.size() == N);
  rsAssert((int) L.size() == N);
  for(int i = 0; i < N; i++)                       // loop over all vertices, i: vertex index
  {
    Vec2 vi = mesh.getVertexData(i);               // location of current vertex i
    Vec2 gi = ND::gradient2D(mesh, &u[0], i);      // gradient estimate at vertex i
    T wSum(T(0));                                  // sum of weights
    T eSum(T(0));                                  // sum of prediction errors
    int M = mesh.getNumEdges(i);                   // number of neighbors of vertex i
    for(int k = 0; k < M; k++)                     // loop over vi's neighbors
    {
      int  j   = mesh.getEdgeTarget(i, k);         // index of current neighbor
      //T    w   = mesh.getEdgeData(i, k);           // edge weight
      T    w   = 1;                                // test
      Vec2 vj  = mesh.getVertexData(j);            // location of current neighbor
      Vec2 dji = vj - vi;                          // difference vector
      T    gij = rsDot(gi, dji);                   // directional derivative in direction dji
      T    pj  = u[i] + gij;                       // prediction of u[j] from u[i] and gij
      T    uj  = u[j];                             // actual value of u[j]
      T    ej  = uj - pj;                          // error = actual - predicted

      // test:
      T h   = rsNorm(dji);                         // distance between vi and vj


      //eSum += w*ej;                                // accumulate error
      eSum += w*ej/h;                              // accumulate error


      wSum += w;                                   // accumulate weights
    }

    T eAvg = eSum / wSum;                          // weighted average of errors...
    L[i]   = eAvg;                                 // ...is estimate for Laplacian

    //L[i]  *= M;  // test
  }

  // doesn't work - maybe we need to divide by the distance h
}

/*
void meshLaplacianErrorVsDistance()
{

  // We test the error of the estimation of the Laplacian via the efficient function...
  // ...not true

}
*/

void meshLaplacianAlgorithms1()
{
  // Compares two methods of computing the Laplacian: by evaluating the Hessian and then taking its 
  // trace (that's the brute force, inefficient way) and using the difference of the value at the 
  // vertex and the (weighted) average of its neighborhood.

  using Vec2 = rsVector2D<double>;
  using Vec  = std::vector<double>;
  using ND   = rsNumericDifferentiator<double>;

  std::function<double(double, double)> f, f_xx, f_yy;
  double a = 0.75;
  double b = 0.5;
  f    = [&](double x, double y)->double { return      sin(a*x) *     exp(b*y); };
  f_xx = [&](double x, double y)->double { return -a*a*sin(a*x) *     exp(b*y); };
  f_yy = [&](double x, double y)->double { return      sin(a*x) * b*b*exp(b*y);  };
  // We have introduced factors a,b, because if they are both 1, the Laplacian happens to become 
  // identically zero which is no good test case

  int numSides = 5;
  double h = 1./8;
  Vec2 x0(1, 1);
  rsGraph<Vec2, double> mesh;

  createMeshForHessianEstimation(mesh, numSides, h, x0);
  GraphPlotter<double> meshPlotter;
  //meshPlotter.plotGraph2D(mesh, {0});

  int N = mesh.getNumVertices();

  // Compute Laplacian by brute force, i.e. computing the full Hessian and then summing the 
  // diagonal elements:
  Vec u_xx(N), u_xy(N), u_yx(N), u_yy(N);
  meshHessianViaGradGrad(mesh, f, u_xx, u_xy, u_yx, u_yy);
  Vec u_L1 = u_xx + u_yy;

  // Compute Laplacian by neighborhood average:
  Vec u(N), u_L2(N);
  fillMeshValues(mesh, f, u);
  ND::laplacian2D_2(mesh, u, u_L2); // this function works only for regular neighborhood geometries

  // Compute Laplacian by convenience function:
  Vec u_L3(N);
  ND::laplacian2D(mesh, &u[0], &u_L3[0]);

  // ToDo: compute Laplacian via Taylor

  // Compute true value and errors:
  double L  = f_xx(x0.x, x0.y) + f_yy(x0.x, x0.y);  // true value
  double e1 = L - u_L1[0];                          // error of first estimate
  double e2 = L - u_L2[0];                          // error of second estimate
  double e3 = L - u_L3[0];                          // same as e1 as it should be

  // Observations:
  // -u_L2[0] is more accurate than u_L1[0], so the more efficient algo is also more accurate
  // -unfortunately, that seems to be the case only if the neighbors are all the same distance
  //  away from the center vertex (see other experiment below) - can the formula be generalized to 
  //  work with more general neighborhood geometries?
  // -the basic idea is that the Laplacian measures by how much the value differs from its local
  //  neighborhood


  // ToDo:
  // -test both algorithms with less regular geometries - maybe create random neighborhoods and 
  //  check, if the fast algo always produces better results than the slow
  //  ...hmm - that raises the question what sort of neighborhoods the neighbors should have 
  //  - maybe instead of estimating 1st derivatives, assign them to exact values - this should make
  //  the Hessian estimate more accurate, so we may need a less restrictive requirement

  // see: https://en.wikipedia.org/wiki/Discrete_Laplace_operator

  int dummy = 0;
}






void meshLaplacianAlgorithms2()
{
  // We test various numerical schemes to estimate the Laplacian for a scalar function defined on a
  // mesh


  // pay special attention to using different distances to the neighbors - will the accuracy be 
  // determined by the distance to the largest neighbor? or maybe by the mean of the distances?
  // -use hexagonal neighborhood and then increase every other distance by factor 2 - will the 
  //  accuracy be in the middle between the two hexagonal neighborhoods with all distances the same


  using Vec2 = rsVector2D<double>;
  using Vec  = std::vector<double>;
  using ND   = rsNumericDifferentiator<double>;

  double h = 1./8;
  int numSides = 6;
  double p = 1.0;     // edge weight exponent/power

  std::function<double(double, double)> f, f_xx, f_yy;
  double a = 0.75;
  double b = 0.5;
  f    = [&](double x, double y)->double { return      sin(a*x) *     exp(b*y); };
  f_xx = [&](double x, double y)->double { return -a*a*sin(a*x) *     exp(b*y); };
  f_yy = [&](double x, double y)->double { return      sin(a*x) * b*b*exp(b*y);  };
  Vec2 x0(1, 1);
  double L = f_xx(x0.x, x0.y) + f_yy(x0.x, x0.y);  // true value of Laplacian at x0

  rsGraph<Vec2, double> mesh;
  GraphPlotter<double> meshPlotter;

  // Compute error for neighborhood at distance h:
  mesh.clear();
  mesh.addVertex(x0);
  addPolygonalNeighbours(mesh, 0, numSides, h, 0., p, false); 
  //meshPlotter.plotGraph2D(mesh, {0});
  int N = mesh.getNumVertices();
  Vec u(N), u_L(N);
  fillMeshValues(mesh, f, u);
  ND::laplacian2D_2(mesh, u, u_L);
  //laplacian2D(mesh, u, u_L);
  double eh1 = L - u_L[0];          // -0.00010707162148965166, -0.00010707162149820038

  // Compute error for neighborhood at distance 2*h:
  mesh.clear();
  mesh.addVertex(x0);
  addPolygonalNeighbours(mesh, 0, numSides, 2*h, 0., p, false); 
  //meshPlotter.plotGraph2D(mesh, {0});
  fillMeshValues(mesh, f, u);
  ND::laplacian2D_2(mesh, u, u_L);
  //laplacian2D(mesh, u, u_L);
  double eh2 = L - u_L[0];          // -0.00042702273039024741, -0.00042702273038902616
  // eh2 is roughly 4 times eh1 for numSides = 6...seems like the error increases with h^2.
  // shouldn't we expect the error to increase by h^4 = h^(numSides-2)?

  // Now contract half of the neighbor-distances by a factor of 0.5, such that half of the neighbor
  // vertices are at distance h and the other half at distance 2*h:
  for(int k = 1; k < N; k+=2)
  {
    Vec2 vk = mesh.getVertexData(k);
    Vec2 dv = vk - x0;
    mesh.setVertexData(k, x0 + 0.5*dv);
  }
  assignEdgeWeights(mesh, p);            // recompute mesh-weights according to new distances
  //meshPlotter.plotGraph2D(mesh, {0});  // plot the mesh
  fillMeshValues(mesh, f, u);
  ND::laplacian2D_2(mesh, u, u_L);
  //laplacian2D(mesh, u, u_L);
  double eh12 = L - u_L[0];   // -0.012610446761559257  ...way too high! something is wrong!
  // could it be that the error is (roughly) equal to taking only the inner neighbors into account?
  // ..hmm...not really - but the order of magnitude seems to fit

  // eh12 should be in between eh1 and eh2 but is much larger than both -> something is still wrong
  // in ND::laplacian2D - the uSum += w*(u[j]-u[i])/d2; in the loop is not the right way to account
  // for the distance(s) h
  // maybe we should form a sort of weighted average:
  // -instead of dividing each term by d2, divide at the very end by the average of d2
  //  -> should nor change behavior, if all d2 are the same
  // -then, divide each term by d2 again - this introdcuces weighting
  // -record the sum of 1/d2, i.e. the total sum of weights
  // -divide final result by total sum of weights
  // -...maybe figure out, if this can be simplified/optimized
  // -hmm...it seems to be better behaved for larger values of numSides - with 6, the error is 
  //  larger as it should be, but for 10 or 12 it is indeed in the expected range
  // -p=1 seems to work well with an even numSides >= 6
  // -for an odd numSides, the error is really large
  // -maybe instead of dividing each term by its own squared distance, we should divide the whole 
  //  sum at the end by the maximum of the squared distances - this will also save a lot of 
  //  divisions

  // try it with the new, experimental implementations:
  laplacian2D_1(mesh, u, u_L);
  double eh12_1 = L - u_L[0]; 

  laplacian2D_2(mesh, u, u_L);
  double eh12_2 = L - u_L[0];    // much better than the 1st

  laplacian2D_3(mesh, u, u_L);
  double eh12_3 = L - u_L[0];    // doesn't seem to work

  laplacian2DViaTaylor(mesh, u, u_L);
  double eh12_T = u_L[0] - L;    // large error :-O


  // shouldn't the error be u_L[0] - L, i.e. estimate minus true?

  //double ratio = L / u_L[0]; // why is this relevant?

  int dummy = 0;

  // -maybe instead of trying to estimate the Laplacian directly, try first to estimate the 
  //  diagonal elements of the Hessian u_xx and u_yy - can we do this via directional derivatives?
  // -maybe the center of gravity of the neighbor vertices should play a role in the formula?
  //  -maybe the weighted average applies the the (weighted) center of gravity and not to the 
  //   actual vertex position and maybe we should subtract a term that is proportional to the 
  //   directional derivative into the direction from the vertex to the COG
  //  -that may explain why the current formula seems to work well for even numSides - in this 
  //   case, the COG coincides with the vertex position
  //  -maybe, to implement it, it will make sense to factor out the body of the i-loop in
  //   rsNumericDifferentiator<T>::gradient2D, so we can compute the gradient u_x, u_y as sub-algo
  //   because we need it to compute the directional derivative into the direction of the COG

  // todo: compute error for h, compute error for 2*h, compute error for half of the neighbors at 
  // distance h and the other half at distance 2*h

  // -try it with a linear function - the Laplacian should be zero


  /*
  int minNumSides =   4;  // maybe go down to 3 or 2 later
  int maxNumSides =  10;
  int numTests    = 100;  // number of tests for each value of numSides
  for(int numSides = minNumSides; numSides <= maxNumSides; numSides++)
  {
    for(int i = 0; i < numTests; i++)
    {


    }
  }
  */
}

void vertexMeshLaplacian()
{
  meshLaplacianAlgorithms1();
  meshLaplacianAlgorithms2();
}



















