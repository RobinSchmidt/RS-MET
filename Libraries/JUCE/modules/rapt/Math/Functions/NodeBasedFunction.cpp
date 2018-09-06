template<class T>
size_t rsNodeBasedFunction<T>::addNode(T x, T y)
{
  // maybe find insertion point and use vector::insert

  nodes.push_back(rsFunctionNode<T>(x, y));

  // keep arrays sorted according to ascending x:
  size_t i = nodes.size()-1;
  //while(i > 0 && isNextValueLess(i-1)) {
  //  rsSwap(nodes[i-1], nodes[i]);
  //  i--; }
  i = moveNodeToSortedIndex(i);
  return constrainNode(i);
}

template<class T>
bool rsNodeBasedFunction<T>::removeNode(size_t index)
{
  if(isNodeRemovable(index)) {
    rsRemove(nodes, index);
    return true; }
  return false;
}

template<class T>
size_t rsNodeBasedFunction<T>::moveNode(size_t i, T newX, T newY)
{
  nodes[i].x = newX;
  nodes[i].y = newY;
  i = moveNodeToSortedIndex(i);
  return constrainNode(i);
}

template<class T>
size_t rsNodeBasedFunction<T>::moveNodeToSortedIndex(size_t i)
{
  while(i > 0 && isNextValueLess(i-1)) { 
    rsSwap(nodes[i-1], nodes[i]); 
    i--; }
  while(i < nodes.size()-1 && isNextValueLess(i)) {
    rsSwap(nodes[i+1], nodes[i]);
    i++; }
  return i;
}

template<class T>
T rsNodeBasedFunction<T>::applyInverseFunction(T y) const
{
  //return y; // preliminary
  //return rsRootFinder<T>::bisection(std::function<T(T)>(*this), getMinX(), getMaxX(), y);
  
  std::function<T(T)> f = *this;
  return rsRootFinder<T>::falsePosition(f, getMinX(), getMaxX(), y);
  
  
  // return rsRootFinder<T>::falsePosition(std::function<T(T)>(*this), getMinX(), getMaxX(), y);
  // false-position method should need only one step in a linear mapping...but may be slow for
  // rational/exponential mappings with a strong bending...maybe a more robust algorithm is needed
  // that switches to bisection steps in cases of slow convergence of false-position - maybe
  // modified false position? see comments in rsRootFinder implementation

  //return rsRootFinder<T>::bisection(*this, getMinX(), getMaxX(), y);
}