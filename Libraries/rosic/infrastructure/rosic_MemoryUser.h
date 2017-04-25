#ifndef rosic_MemoryUser_h
#define rosic_MemoryUser_h

// rosic-indcludes:
#include "../basics/GlobalFunctions.h"

namespace rosic
{

  /**

  This class serves as a baseclass for objects that need some memory to be allocated and used. Its 
  intended purpose is to be able to share memory areas between different objects when not all of 
  them are needed simultaneously (as for example in Quadrifex).

  */

  class MemoryUser
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. Optionally allocates some memory area to be used. */
    MemoryUser(int sizeToAllocateInBytes = 0);   

    /** Destructor. */
    ~MemoryUser();  

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets the memory area to be used. If this object already has some memory aallocated, this 
    will either be freed or not, depending on the state of the wasAllocatedHere flag. */
    void setMemoryAreaToUse(void *newStartAddress, int newSizeInBytes);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the start address of the used memory area. */
    void* getStartAddress() const { return startAddress; }

    /** Returns the size (in bytes) of the used memory area. */
    int getSizeInBytes() const { return sizeInBytes; }

    //=============================================================================================

  protected:

    /** Allocates some area of memory and sets the wasAllocatedHere flag true. */
    void allocateMemory(int sizeToAllocateInBytes);

    /** Frees the allocated memory area but only if it was allocated by this object. */
    void freeMemoryIfItWasAllocatedHere();

    void *startAddress;
    int  sizeInBytes;
    bool wasAllocatedHere;

  private:

    /** Frees the allocated memory area. */
    void freeMemory();

  };

} // end namespace rosic

#endif // rosic_MemoryUser_h
