//=================================================================================================

template<class T>
void rsMeshWeightCalculator2D<T>::initEdgeWeights(rsGraph<rsVector2D<T>, T>& mesh)
{
  for(int i = 0; i < mesh.getNumVertices(); i++)
    for(int k = 0; k < mesh.getNumEdges(i); k++)
      mesh.setEdgeData(i, k, T(1));
}

template<class T>
void rsMeshWeightCalculator2D<T>::weightEdgesByDistances(rsGraph<rsVector2D<T>, T>& mesh, T p)
{  
  using Vec2 = rsVector2D<T>;
  T q = p;
  for(int i = 0; i < mesh.getNumVertices(); i++)
  {
    Vec2 v_i = mesh.getVertexData(i);
    if(rsIsNaN(p))
      q = (T)mesh.getNumEdges(i);
    for(int k = 0; k < mesh.getNumEdges(i); k++)
    {
      int  j   = mesh.getEdgeTarget(i, k);    // index of target vertex
      Vec2 v_j = mesh.getVertexData(j);       // position of target vertex
      T d = rsNorm(v_j - v_i);                // distance between v_i and v_j
      T w = mesh.getEdgeData(i, k);
      w  *= pow(d, -q);
      mesh.setEdgeData(i, k, w);
      int dummy = 0;
    }
  }
  // ToDo: 
  //  -instead of w =...; w *= ...; mesh.setEdge...; use a convenience function 
  //   mesh.scaleEdgeData(i, k, s) which does: setEdgeData(i, k, s*getEdgeData(i, k))
  //   where s = pow(d, -q);

  int dummy = 0;
}

template<class T>
T rsMeshWeightCalculator2D<T>::getNeighborSeparation(rsGraph<rsVector2D<T>, T>& mesh, int i, int k)
{
  // Computes a measure of seperateness of the k-th neighbor of node i.

  using Vec2 = rsVector2D<T>;

  int  j    = mesh.getEdgeTarget(i, k);
  Vec2 v_j  = mesh.getVertexData(j);
  T    s_ik = T(0);
  for(int m = 0; m < mesh.getNumEdges(i); m++)
  {
    if(m != k)       // maybe try it without this conditional
    {
      int  n    = mesh.getEdgeTarget(i, m);
      Vec2 v_n  = mesh.getVertexData(n);
      T    d_jn = rsNorm(v_n - v_j);
      s_ik += d_jn;  // maybe we should use some power p: pow(d_jn, p)
                     // maybe try multiplative accumulation (init with 1)
    }
  }

  return s_ik;
}

template<class T>
void rsMeshWeightCalculator2D<T>::weightEdgesByMutualDistance(rsGraph<rsVector2D<T>, T>& mesh)
{
  using Vec2 = rsVector2D<T>;
  T q(+1.00);  // todo: make parameter

  for(int i = 0; i < mesh.getNumVertices(); i++)
  {
    Vec2 vi = mesh.getVertexData(i);
    for(int k = 0; k < mesh.getNumEdges(i); k++)
    {
      T w_ik  = mesh.getEdgeData(i, k);
      T s_ik  = getNeighborSeparation(mesh, i, k);
      w_ik   *= pow(s_ik, q);
      mesh.setEdgeData(i, k, w_ik);
      int dummy = 0;
    }
  }

  int dummy = 0;
}

template<class T>
void rsMeshWeightCalculator2D<T>::weightEdgesByPositions(
  rsGraph<rsVector2D<T>, T>& mesh, int formula)
{
  if(formula == 0) {                                    return; }
  if(formula == 1) { weightEdgesByMutualDistance(mesh); return; }

  // ToDo: implement a formula using the mutual correlations of the neighbors c(k,m)
}

template class rsMeshWeightCalculator2D<double>;  // explicit instantiation for double

// ToDo:
// -Maybe move all functions that have to do with computing or manipulating edge-weighst into a 
//  class, maybe rsMeshWeightCalculator or something
// -Move the weightEdgesByDirection from the research codebase to here

//=================================================================================================

template<class T>
void randomizeVertexPositions(rsGraph<rsVector2D<T>, T>& mesh, T dx, T dy, 
  int minNumNeighbors, int seed)
{
  using Vec2 = rsVector2D<T>;
  rsNoiseGenerator<T> ng;
  ng.setSeed(seed);
  for(int k = 0; k < mesh.getNumVertices(); k++) 
  {
    if(mesh.getNumEdges(k) >= minNumNeighbors)
    {
      Vec2 v = mesh.getVertexData(k);
      v.x += dx * ng.getSample();
      v.y += dy * ng.getSample();
      mesh.setVertexData(k, v);
    }
  }
}
template void randomizeVertexPositions(rsGraph<rsVector2D<double>, double>& mesh, 
  double dx, double dy, int minNumNeighbors, int seed);


template<class T>
void taylorExpansion2D(const rsGraph<rsVector2D<T>, T>& mesh, const T* u, 
  T* u_x, T* u_y, T* u_xx, T* u_xy, T* u_yy)
{
  using Vec2 = rsVector2D<T>;
  using Mat  = rsMatrix<T>;
  static const T half = (T(1)/T(2));
  Mat A, x, b;                                     // for the matrix equation
  x.setShape(5, 1);                                // x := (g_x, g_y, H_xx, H_xy, H_yy)
  for(int i = 0; i < mesh.getNumVertices(); i++)
  {
    Vec2 vi = mesh.getVertexData(i);
    int  K  = mesh.getNumEdges(i);                 // number of neighbors of vertex i
    A.setShape(K, 5);
    b.setShape(K, 1);
    for(int k = 0; k < K; k++)
    {
      int  j  = mesh.getEdgeTarget(i, k);
      Vec2 vj = mesh.getVertexData(j);
      Vec2 dj = vj   - vi;
      b(k, 0) = u[j] - u[i];
      A(k, 0) = dj.x;
      A(k, 1) = dj.y;
      A(k, 2) = dj.x * dj.x * half;
      A(k, 3) = dj.x * dj.y;
      A(k, 4) = dj.y * dj.y * half;
    }
    solveOptimal(A, x, b);                         // allocates - ToDo: avoid that
    u_x[i]  = x(0, 0);
    u_y[i]  = x(1, 0);
    u_xx[i] = x(2, 0);
    u_xy[i] = x(3, 0);
    u_yy[i] = x(4, 0);
  }

  // -To introduce weighting, would be be sufficient to scale b(k,0), A(k,_) by w_ij? Or maybe
  //  we need to square the weight for the coeffs involving products like dj.x^2 ...but no - that
  //  makes no sense - we can only scale a matrix-row uniformly (together with the corresponding 
  //  rhs)
  // -I think, to construct the matrix, we can ignore u and instead need to form the pseudo-inverse
  //  of A. From that, we can exctract the coeffs that go into the matrices for computing u_x, u_y
  //  u_xx, u_xy, u_yy from u...oh but we may have to pre- or post multiply by A^T
  // -For the underdetermined cases (2,3,4 neighbors) try to split the Hessian in a different way
  //  involving Laplacian, etc ...see text in the private repo
  // -For the overdetermined case, try a weighted least squares using weights inversely 
  //  proportional to (some power of) the distance

  // ToDo:
  // -Make a similar function, but instead of estimating u_xx and u_yy, etimate u_xx + u_yy and
  //  u_xx - u_xy. I think this requires the change:
  //    A(k, 2) = (dj.x * dj.x + dj.y * dj.y) * half;  // not sure about the half...
  //    A(k, 3) = (dj.x * dj.x - dj.y * dj.y) * half;  // dito
  //  rationale: u_xx + u_yy is better interpretable (it's the Laplacian)
  // -Consider the function as a 3D surface: (x(u,v), y(u,v), z(u,v)) = (u, v, f(u,v)). Warning:
  //  in surface notation, u means our first independent parameter (i.e. 1st coordinate), in the 
  //  notation of PDE solvers, u means the function vaules, i.e. the height z ...maybe change 
  //  notation u, u_x, etc. to z, z_x etc. - on the other hand, we here have x(u,v) = u, 
  //  y(u,v) = v anyway, so we don't actually need a notation for u and v and could instead just 
  //  use x,y directly. We need to be careful when applying formulas from differential geometry, 
  //  then. We would then have for the 1st partial derivatives of the surface with respect to the
  //  coordinates: x_u = y_v = 1, x_v = y_u = 0, z_u = f_x, z_v = f_y. Whatever we end up doing: 
  //  watch out for this potoential source of confusion.
  // -Instead of using matrix entries a,b,c,d, express a 2x2 matrix in terms of its invariants, 
  //  such as trace, determinant, eigenvalues/vectors
  //    https://en.wikipedia.org/wiki/Invariant_(mathematics)
  //  maybe certain norms can be used, too:
  //    https://en.wikipedia.org/wiki/Matrix_norm#Matrix_norms_induced_by_vector_norms
  //  ...it says there that for symmetric matrices "the 2-norm is precisely the spectral radius".
  //  We have b=c ...so maybe instead of estimating a,b,d = H_xx, 2*H_xy, H_yy (2* could be /2)
  //  we should try to estimate Tr(A), det(A), rho(A). see also:
  //    https://nhigham.com/2021/02/02/what-is-a-unitarily-invariant-norm/
  //  The Frobenius norm is unitarily invariant. In general, for a 2x2 matrix A = (a,b; c,d), we 
  //  can compute the following "interesting numbers" to characterize the matrix in a somewhat
  //  coordinate independent way:
  //    T := a+d                         trace, twice the mean of eigenvalues
  //    D := a*d - b*c                   determinant
  //    F := a^2 + b^2 + c^2 + d^2       Frobenius norm
  //    V := a*a + 4*b*c - 2*a*d + d*d   twice the variance of eigenvalues, discriminant (my term)
  //    R := (T + sqrt(V))/2             larger eigenvalue, spectral radius, maximum grow factor
  //    S := (T - sqrt(V))/2             smaller eigenvalue, maximum shrink factor
  //  Maybe try to use T,D,V instead of a,b,d for a symmetric 2x2 matrix. The discriminant V says 
  //  something about both eigenvalues. In fact, if we also know T, we can comptute them both. Or
  //  try (T,D,F), (T,V,F)...try various things and figure out which combination gives better 
  //  estimations of values of a quadratic form...that may, of course, depend on the particular 
  //  function chosen. The Frobenius norm seems to be invariant only under unitary transformations,
  //  i.e. rotations...but maybe that's good enough. Are the others actually indeed invariant under
  //  any change opf basis? -> figure out, do also numerical tests for this. In a more general 
  //  setting where we don't assume b = c, maybe use (T,D,F,V) instead of (a,b,c,d). Solve the 4
  //  equations for a,b,c,d to find formulas for converting back
  // -What about the gradient vector? Can we characterize it in a coordinate independent way as 
  //  well? It's length is an invariant. Is there another invariant by which we can reconstruct its 
  //  components when the length is also known?
  // -What other interesting coordinate independent geometric properties could there be? What about 
  //  a measure of "twist" - maybe the difference in direction of the gradient: take gradient g0 at 
  //  v0, normalize it to get g0n, take gradient at v0+h*g0n, normalize that again, take difference
  //  to g0n, let h go to zero
  // -If we assume c=b, solve the T-equation for a, plug the expression for a into the D-equation 
  //  and solve for b^2 (=b*c) and then plug the results for a and b^2 into the V-equation, d drops 
  //  out and we are left with V = T^2 - 4*D. Similarly, using the F-equation instead of the 
  //  V-equation, again d drops out and we get F = T^2 - 2*D. That means, in case of a symmetric 
  //  matrix, the invariants F and V are not independent from T and D. So it seems, in order to be 
  //  able to reconstruct a,b,d, we need a 3rd, independent invariant. Is there such a thing? Maybe
  //  try these things:
  //  -angle between eigenvectors
  //  -norms of the columns: a^2 + c^2, b^2 + d^2
  //  -inner product of the 2 columns (seen as vectors): a*b + c*d
}
template void taylorExpansion2D(const rsGraph<rsVector2D<double>, double>& mesh, const double* u, 
  double* u_x, double* u_y, double* u_xx, double* u_xy, double* u_yy);



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
template void fillEdges(rsGraph<rsVector2D<double>, double>& g, 
  const std::function<double(rsVector2D<double>, rsVector2D<double>)>& f);
template void fillEdges(rsGraph<rsVector2D<float>, float>& g, 
  const std::function<float(rsVector2D<float>, rsVector2D<float>)>& f);

template<class T>
void addPolygonalNeighbours(rsGraph<rsVector2D<T>, T>& mesh, int i, int numSides, T radius, 
  T angle, T p, bool symmetric)
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

template<class T>
void createPolygonMesh(rsGraph<rsVector2D<T>, T>& mesh, int numSides, T radius, 
  rsVector2D<T> center, T angle)
{
  mesh.clear();
  mesh.addVertex(center);
  addPolygonalNeighbours(mesh, 0, numSides, radius, angle); 
}
template void createPolygonMesh(rsGraph<rsVector2D<double>, double>& mesh, int numSides, 
  double radius, rsVector2D<double> center, double angle);