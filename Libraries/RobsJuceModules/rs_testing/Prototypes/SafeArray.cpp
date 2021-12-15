



/*


rsNonReAllocatingArray:

ToDo:
-The container capacity should alway be a power of two, starting with some small'ish inital 
 capacity like 8.
-When growing the capacity, we always allocate a new chunk whose size is equal to the current
 total size.
-I think, due to the power of 2 chunk-sizes, we can do some bit-masking and shifting to figure 
 out the chunk-index j and element index k from the flat index i in []: first, we do a 
 right-shift by log2 of chunk-size 0, i.e. chunks[0].size. This should give the chunk index j.
 the element-index k is i minus the sum of all chunk-sizes of the chunks 0...j-1 which is 
 2^(j-m) where m is the log2 of the initial size?
 ...wait...really - should the right-shift give j? why? i think, we need the position l of the 
 leftmost nonzero bit. this is related to j. i think j = l - log2(C0) where C0 is the inital 
 capacity?

Let's look an example where the initial size, i.e. the size of chunk[0] is 4:
 chunk:  0   1   2   3   4   5
 size:   4   4   8   16  32  64  
 total:  4   8   16  32  64  128
and we want to figure out chunk index j and element index k for a flat index of i = 26 which is
11010 in binary. We have j=3, k=10: the element is at index 10 in chunk 3. Let l be the position,
counted from right, of the leftmost nonzero bit in i. we have l=5. From that, we get j as 
j = l-2 where the -2 comes from log2(4) where 4 is our initial capacity C0. The k is just 
k = 26-16 = 10 where the 16 comes from 2^(j+1) i think?

-try to make it STL compatible
-what about concatenationg two such arrays? maybe that would need a new datastructure? or maybe we 
 just appaend copies of the 2nd to the 1st where the 1st may grow?




SafeArray (maybe rename to rsDoublyBufferedArray):

The idea is that we maintain two copies of the data which are always in sync except during 
mutation. During mutation, a client will see the old copy which is kept in a valid state while the
other, new copy of the data is being prepared. When the mutation is complete, we atomically swap
the numElems member and the data-pointer (each of the two swaps chouls be atomic). If the mutation
has increased the array size, the data pointer is updated first and numElems is updated second. If 
the mutation has decreased the array size, the order is the other way around. This ensures that the
client always sees the shorter length in cases when a context switch occurs during mutation, so it 
will never try to access the underlying vectors beyond their end. It may happen that the client 
sees a too short length in which case it may miss the last element while looping over the array but
it will never try to loop beyond the end. The intended use is in situations where it's acceptable
to miss looping over the last element - for example in EasyQ, when the user adds a new band
...tbc...

ToDo:
-Implement two variants of a thread-safe dynamically allocated array with the following proerties:
 -read-access to the array's size and to an element is lock-free
 -adding or removing elements may acquire a lock
 -replacing elements may acquire a lock
-One variant should be based on std::vector, the other on new/delete. The former is for debugging, 
 the latter may or may not be swapped in in performance critical production code
-It should be a drop-in replacement for std::vector, but thread-safe.





rsArrayView:

The idea is to wrap a raw C-array into an STL compatible container to make the STL algorithms 
available for it. But actually, for many algorithms, we can just use pointers instead of more 
general iterators when calling the std::algorithms, so thi si not necessary. ..but anyway, it may
be a good excercise in writing an STL compliant container

https://medium.com/@vgasparyan1995/how-to-write-an-stl-compatible-container-fc5b994462c6
https://github.com/vgasparyan1995/prefix_tree/tree/master/include


https://stackoverflow.com/questions/5454127/stl-compliant-container
https://www.codeproject.com/Articles/4736/STL-compliant-container-example


https://github.com/electronicarts/EASTL/
Implementation of the STL that was called "fairly readable" somewhere
...hmm - vector.h just typedefs iterator as T* and does not actually implement an internal
class for iterators:
  https://github.com/electronicarts/EASTL/blob/master/include/EASTL/vector.h
maybe we should do the same

https://www.fluentcpp.com/2018/05/08/std-iterator-deprecated/

-implement a non-reallocating dynamic array here

*/