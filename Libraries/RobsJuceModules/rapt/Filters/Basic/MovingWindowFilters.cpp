

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

rsMovingMaximumFilter:
-maybe we could further reduce the worst case processing cost by using binary search to
 adjust the new tail-pointer instead of linearly popping elements one by one? search
 between tail and head for the first element that is >= in (or not < in) - but maybe that
 would destroy the amortized O(1) cost? if we would do the binary search for every sample, we would
 incur a log(L) cost for every sample where L is the current length of the queue. maybe the linear
 search is guranteed to encounter a length L whose average value, averged over many samples, is 
 constant? or that the element we look for is in within the first M samples at that M is const on
 the average? that seems to make sense. how about a search strategy that is a hybrid between linear 
 and binary search, like this:
 -look at the first k (say, k=8) elements by linear search
 -if the element was not found, look at the next 8 with binary search
 -if not found, look at the next 16 with binary search
 -if not found, look at the next 32 with binary search
 -and so on
 -in the worst case, the element is in the last segment that we look at and it takes us 
  O(log(8)) + O(log(16)) + O(log(32)) + O(log(64)) + ...to get to that segment and the search 
  itself within the segment is O(log(N)). let's for simplicity assume k = 1, so the sum starts at
  log(1) = 0 and looks like: 0+1+2+3+4+5+6+..+log(N/2) <= 2*log(N/2), so getting to the last 
  segment is also O(log(N))
 -it seems that such a search strategy could combine the advantage of linear search (amortized 
  O(1) cost) with the advantage of binary search (worst case O(log(N)) cost)
 -i have no idea, if that is really correct -> figure out


maybe using a "greater" function instead of a > operator is not such a good idea,
performance-wise - we should do performance tests and maybe for production code, provide
a version of getSample that just uses >

*/