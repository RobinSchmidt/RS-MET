template<class T>
void randomizeVertexPositions(rsGraph<rsVector2D<T>, T>& mesh, T dx, T dy, 
  int minNumNeighbors, int seed)
{
  using Vec2 = rsVector2D<T>;
  rsNoiseGenerator<T> ng;
  ng.setSeed(seed);
  T rnd;
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