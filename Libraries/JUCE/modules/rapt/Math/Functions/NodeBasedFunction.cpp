template<class T>
size_t rsNodeBasedFunction<T>::addNode(T x, T y)
{
  // maybe find insertion point and use vector::insert

  nodes.push_back(rsFunctionNode<T>(x, y));

  // keep arrays sorted according to ascending x:
  size_t i = nodes.size()-1;
  while(i > 0 && isNextValueLess(i-1)) {
    rsSwap(nodes[i-1], nodes[i]);
    i--; }
  return constrainNode(i);
}

template<class T>
void rsNodeBasedFunction<T>::removeNode(size_t index)
{
  if(isNodeRemovable(index))
    rsRemove(nodes, index);
}

template<class T>
size_t rsNodeBasedFunction<T>::moveNode(size_t i, T newX, T newY)
{
  nodes[i].x = newX;
  nodes[i].y = newY;

  // keep arrays sorted according to ascending x:
  while(i > 0 && isNextValueLess(i-1)) { 
    rsSwap(nodes[i-1], nodes[i]); 
    i--; }
  while(i < nodes.size()-1 && isNextValueLess(i)) {
    rsSwap(nodes[i+1], nodes[i]);
    i++; }
  return constrainNode(i);
}