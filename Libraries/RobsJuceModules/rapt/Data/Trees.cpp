
template<class T>
bool rsBinaryHeap<T>::isHeap(int i) const
{
  if(i >= this->size)
    return true;
  bool result = true;
  int l = this->left(i);
  int r = this->right(i);
  if(l < this->size) result &= !less(this->data[i], this->data[l]) && isHeap(l);
  if(r < this->size) result &= !less(this->data[i], this->data[r]) && isHeap(r);
  return result;
}

template<class T>
int rsBinaryHeap<T>::insert(const T& x)
{
  rsAssert(this->size < this->capacity, "Capacity exceeded");
  if( this->size == this->capacity ) return -1;
  // Trying to insert an item when the heap is already full is a bug on client code side. We
  // return -1 early here, to avoid an access violation if client code has this bug in a
  // release version.

  this->data[this->size] = x;
  this->size++;
  return floatUp(this->size-1);
}

template<class T>
void rsBinaryHeap<T>::remove(int i)
{
  rsAssert(i >= 0 && i < this->size, "Index out of range");
  if(i < 0 || i >= this->size)
    return;  // ...it's also a client code bug to try to remove a non-existent item.
  swap(this->data[i], this->data[this->size-1]);
  this->size--;
  floatIntoPlace(i);
}
// Could the return value of floatIntoPlace be interesting for the caller? It is the position
// where the fomerly last heap element ended up after all the re-ordering business - we'll see

template<class T>
int rsBinaryHeap<T>::floatUp(int i)
{
  while(i > 0) {
    int p = this->parent(i);
    if(less(this->data[p], this->data[i]))  {
      swap(this->data[i], this->data[p]);
      i = p; }
    else
      return i; }
  return i;
}

template<class T>
int rsBinaryHeap<T>::floatDown(int i)
{
  while(i < this->size-1)   // check if we should use size
  {
    int l = this->left(i);
    int r = this->right(i);  // == l+1  ->  optimize (but maybe parallel is actually better than serial?)
    int b = i;
    if(l < this->size && less(this->data[i], this->data[l])) b = l;
    if(r < this->size && less(this->data[b], this->data[r])) b = r;
    if(b != i) {
      swap(this->data[i], this->data[b]);
      i = b;  }
    else
      return i;
  }
  return i;
}

//=================================================================================================

#undef small

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

/*

Ideas:
In a binary max-heap, the property to be satisfied is that both children of a node must be <= the
node itself, but it doesn't specify the order of the children. Would it be easy and/or useful to
impose the additional condition, that, for example, the left child should be <= or >= the right
child? That would give even more structure to the data. Or maybe such a condition does already
hold due to the implementation details of floatUp/Down? If not, maybe a small change to those
function could have that effect? When we float a node down, we would always float it down the path
that starts with the smaller child...or something? hmm - one step in the floatDown algo just looks
at a node's two children and potentially swaps it with the larger of its children. If it holds as
precondition, that the larger child is the right child...the swap could change that - and to fix
it, we would have to swap whole subtrees, right? That would not be O(log(N)) anymore because due
to the implicitness of our tree (parent/child relationships defined by array-index), we would have
to move a lot of data around which would be O(N). However, if the tree is implemented as a linked
tree of nodes rather than an implicit tree, we could indeed do such a swap of subtrees in O(1).
..maybe such a structure could then be called linked heap or "leap" - sounds cool :-)

Maybe implement a k-ary heap:
-each node has k children
-all children are <= their parent
-children of node at index i are found at k*i+1, k*i+2, ..., k*i+k and the parent is at (i-1)/k
-advantage is that the height of the trees tend to get smaller when k gets larger (height is
 proportional to the base-k logarithm) but in floatUp/Down, we must at each node find the maximum
 out of k values



see:
https://en.wikipedia.org/wiki/Self-balancing_binary_search_tree
https://en.wikipedia.org/wiki/Heap_(data_structure)
https://en.wikipedia.org/wiki/Implicit_data_structure
https://en.wikipedia.org/wiki/Binary_tree#perfect
https://en.wikipedia.org/wiki/Order_statistics

https://opendatastructures.org/ods-cpp/10_1_Implicit_Binary_Tree.html
https://opendatastructures.org/ods-cpp/ods-cpp-html.html

https://algo.ics.hawaii.edu/pubs/18-ipdps.pdf

https://www.geeksforgeeks.org/binary-heap/
https://www.geeksforgeeks.org/heap-data-structure/

https://www.cs.csustan.edu/~john/Classes/CS4440/Notes/02_Heaps+PQs/heaps.html
http://staffwww.fullcoll.edu/aclifton/cs133/lecture-14-trees-and-heaps.html looks very comprehensive

https://cglab.ca/~morin/teaching/5408/refs/minmax.pdf

*/
