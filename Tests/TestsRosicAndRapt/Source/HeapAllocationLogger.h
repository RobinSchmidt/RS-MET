#ifndef HEAPALLOCATIONLOGGER_H_INCLUDED
#define HEAPALLOCATIONLOGGER_H_INCLUDED

/*

Unfotunatly, this doesn't work. At least not when trying to include this file before everything 
else. In MSVC, I get compiler errors linked to the lines:

using _CSTD free;
using _CSTD malloc;

cstdlib

Also, the functions rsLoggingMalloc/Free don't compile because malloc/free are not yet defined.
I think, it is not right to include this file as the very first thing. malloc and free have to be
already defined. If this is to work at all, I think, it has to included in between the files that
define malloc/free and eveything else. 


Maybe it could be made to work by including it after that file?

*/



/** A singleton object for logging heap allocations. */

class rsHeapAllocationLogger
{

public:

  static rsHeapAllocationLogger* getInstance()
  {
    if(theObject == nullptr)
      theObject = new rsHeapAllocationLogger;
    return theObject;
  }

  static void deleteInstance()
  {
    delete theObject;
  }

  void logAllocation() {  numAllocs++;  }

  int getNumAllocations() const { return numAllocs; }


  void logDeallocation() { numDeallocs++; }

  int getNumDeallocations() const { return numDeallocs; }


  int getNumAllocatedChunks() const { return getNumAllocations() - getNumDeallocations(); }


  void reset()
  {
    numAllocs   = 0;
    numDeallocs = 0;
  }


private:

  rsHeapAllocationLogger(){};

  static rsHeapAllocationLogger* theObject; // = nullptr;

  int numAllocs   = 0;
  int numDeallocs = 0;

};

rsHeapAllocationLogger* rsHeapAllocationLogger::theObject = nullptr;


// https://stackoverflow.com/questions/1008019/how-do-you-implement-the-singleton-design-pattern


// Maybe also log the allocated size. But this requires us to keep track of all the addresses of 
// the allocated chunks and their sizes, so it would compicate the implementation a lot. That's 
// overkill at the moment. 




/*
void* rsLoggingMalloc(size_t* size)
{
  rsHeapAllocationLogger::getInstance()->logAllocation();
  return malloc(size);

  // See: https://en.cppreference.com/w/c/memory/malloc
}

void rsLoggingFree(void* ptr)
{
  rsHeapAllocationLogger::getInstance()->logDeallocation();
  free(ptr);

  // See: https://en.cppreference.com/w/c/memory/free
}
*/



//#define malloc(size) (rsLoggingMalloc(size))
//#define free(ptr)    (rsLoggingFree(ptr))
// I think, this also covers new, new[], delete, delete[], because they use malloc and free 
// internally...but cant we really count on that? Also, what about realloc and calloc?







// https://stackoverflow.com/questions/438515/how-to-track-memory-allocations-in-c-especially-new-delete
// https://en.wikipedia.org/wiki/Electric_Fence
// https://stackoverflow.com/questions/9702292/overriding-malloc-to-log-allocation-size
// https://stackoverflow.com/questions/262439/create-a-wrapper-function-for-malloc-and-free-in-c


// https://valgrind.org/
// 

#endif