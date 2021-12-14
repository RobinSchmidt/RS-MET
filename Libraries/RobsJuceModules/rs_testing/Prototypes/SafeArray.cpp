



/*

ToDo:
-Implement two variants of a thread-safe dynamically allocated array with the following proerties:
 -read-access to the array's size and to an element is lock-free
 -adding or removing elements may acquire a lock
 -replacing elements may acquire a lock
-One variant should be based on std::vector, the other on new/delete. The former is for debugging, 
 the latter may or may not be swapped in in performance critical production code
-It should be a drop-in replacement for std::vector, but thread-safe.

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



*/