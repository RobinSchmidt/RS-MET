//#include "rosic_MemoryUser.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

MemoryUser::MemoryUser(int sizeToAllocateInBytes)
{
  startAddress     = NULL;
  wasAllocatedHere = false;
  if( sizeToAllocateInBytes >= 0 )
  {
    sizeInBytes = sizeToAllocateInBytes;
    if( sizeToAllocateInBytes > 0 )
      allocateMemory(sizeInBytes);
  }
  else
  {
    sizeInBytes = 0;
    DEBUG_BREAK;
  }
}

MemoryUser::~MemoryUser()
{
  freeMemoryIfItWasAllocatedHere();
}

//-------------------------------------------------------------------------------------------------
// memory setup:

void MemoryUser::setMemoryAreaToUse(void *newStartAddress, int newSizeInBytes)
{
  freeMemoryIfItWasAllocatedHere();
  startAddress     = newStartAddress;
  sizeInBytes      = newSizeInBytes;

  wasAllocatedHere = false;
}

void MemoryUser::allocateMemory(int sizeToAllocateInBytes)
{
  freeMemoryIfItWasAllocatedHere();

  if( sizeToAllocateInBytes < 0 )
  {
    sizeToAllocateInBytes = 0;
    DEBUG_BREAK;
  }

  if( sizeToAllocateInBytes > 0 )
  {
    sizeInBytes  = sizeToAllocateInBytes;
    startAddress = malloc(sizeInBytes);
    if( startAddress == NULL )
    {
      DEBUG_BREAK;    // memory allocation failed
      sizeInBytes = 0;
    }
  }
  else
  {
    startAddress = NULL;
    sizeInBytes  = 0;
  }

  wasAllocatedHere = true;
}

void MemoryUser::freeMemoryIfItWasAllocatedHere()
{
  if( startAddress != NULL && wasAllocatedHere == true )
    freeMemory();
}

void MemoryUser::freeMemory()
{
  free(startAddress);  
  startAddress = NULL;
  sizeInBytes  = 0;
}

