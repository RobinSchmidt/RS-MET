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