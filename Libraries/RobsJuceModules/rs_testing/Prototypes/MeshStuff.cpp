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
void taylorExpansion2D(rsGraph<rsVector2D<T>, T>& mesh, const T* u, 
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
}
template void taylorExpansion2D(rsGraph<rsVector2D<double>, double>& mesh, const double* u, 
  double* u_x, double* u_y, double* u_xx, double* u_xy, double* u_yy);