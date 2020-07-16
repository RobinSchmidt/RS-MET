

template<class T>
rsMovingMaximumFilter<T>::rsMovingMaximumFilter(size_t maxLength) 
  : delayLine(maxLength), maxDeque(maxLength)
{
  //greater = rsGreater;
  //greater = rsGreater<const T&, const T&>;
  //greater = std::function<bool(const T&, const T&)>(&rsGreater);
}


/*
ToDo:

try to implement a moving median filter, see here:

https://www.youtube.com/watch?v=VmogG01IjYc&list=PLI1t_8YX-Apv-UiRlnZwqqrRT8D1RhriX&index=11

but this implementation never forgets old samples

-at each sample, we must 
 -accept the new incoming sample and place it into either the bucket of values smaller than the 
  median or the bucket of values greater than the median (implemented as min- and max-heap), the 
  max heap contains the smaller values, so that we can access the largest of all smaller values in
  constant time
 -remove the oldest sample
 -i think, we need to keep an array of values sorted by age as circular buffer that somehow points
  to the nodes of the heap...the heap nodes also need to reference back to the values in this 
  array, so we can update the array, when we re-order the heap
-before and after the insertion and removal, the size of both buckets must be the same
-insertion/deletion should be done in log(N) time, i hope

*/