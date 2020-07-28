
template<class T>
bool rsBinaryHeap<T>::isHeap(int i) const
{
  if(i >= size)
    return true;
  bool result = true;
  int l = left(i);
  int r = right(i);
  if(l < size) result &= !less(data[i], data[l]) && isHeap(l);
  if(r < size) result &= !less(data[i], data[r]) && isHeap(r);
  return result;
}

template<class T>
int rsBinaryHeap<T>::insert(const T& x)
{
  rsAssert(size < capacity, "Capacity exceeded");
  if( size == capacity ) return -1; 
  // Trying to insert an item when the heap is already full is a bug on client code side. We 
  // return -1 early here, to avoid an access violation if client code has this bug in a 
  // release version.

  data[size] = x;
  size++;
  return floatUp(size-1);
}

template<class T>
void rsBinaryHeap<T>::remove(int i) 
{ 
  rsAssert(i >= 0 && i < size, "Index out of range");
  if(i < 0 || i >= size)
    return;  // ...it's also a client code bug to try to remove a non-existent item.
  swap(data[i], data[size-1]); 
  size--; 
  floatIntoPlace(i);
}
// Could the return value of floatIntoPlace be interesting for the caller? It is the position
// where the fomerly last heap element ended up after all the re-ordering business - we'll see

template<class T>
int rsBinaryHeap<T>::floatUp(int i)
{
  while(i > 0) {
    int p = parent(i);
    if(less(data[p], data[i]))  {
      swap(data[i], data[p]);
      i = p; }
    else
      return i; }
  return i;
}

template<class T>
int rsBinaryHeap<T>::floatDown(int i)
{
  while(i < size-1)   // check if we should use size
  {
    int l = left(i);
    int r = right(i);  // == l+1  ->  optimize (but maybe parallel is actually better than serial?)
    int b = i; 
    if(l < size && less(data[i], data[l])) b = l;
    if(r < size && less(data[b], data[r])) b = r;
    if(b != i) { 
      swap(data[i], data[b]);
      i = b;  }
    else
      return i;
  }
  return i;
}

//=================================================================================================

template<class T>
int rsDoubleHeap<T>::replace(int key, const T& newValue)
{
  rsAssert(isKeyValid(key), "Key out of range");
  int k = key;

  // The actual replacement:
  if(isKeyInLargeHeap(k)) {
    k  = large.replace(toLargeHeapIndex(k), newValue);
    k |= firstBitOnly; }  
  else
    k = small.replace(k, newValue);

  // The potential swap:
  if(small.less(large[0], small[0])) {
    small.swap(small[0], large[0]);
    int ks = small.floatDown(0);
    int kl = large.floatDown(0);
    if(isKeyInLargeHeap(key)) // new value was in large and is now in small heap
      k = ks;
    else
      k = kl | firstBitOnly;  }

  return k; // return the new key where the replacement node ended up
}