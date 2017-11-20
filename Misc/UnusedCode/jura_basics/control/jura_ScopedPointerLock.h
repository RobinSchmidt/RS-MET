#ifndef jura_ScopedPointerLock_h
#define jura_ScopedPointerLock_h

/**

This class works analogously to juce::ScopedLock, except that it works with a 
pointer-to-CriticalSection instead of a const reference. That makes it suitable to use with 
pointers to CriticalSection objects, where the pointer may also be NULL. In this case, a regular
ScopedLock would try to dereference a NULL-pointer whereas this object here does nothing (i.e. it 
verifies that the pointer is not NULL before dereferencing it).

*/

class JUCE_API ScopedPointerLock
{
public:

  inline ScopedPointerLock(const CriticalSection *lock) throw() : lock_(lock)
  {
    if(lock_ != NULL)
      lock->enter();
  }

  inline ~ScopedPointerLock() throw()
  {
    if(lock_ != NULL)
      lock_->exit();
  }

private:

  const CriticalSection *lock_;

  ScopedPointerLock(const ScopedPointerLock*);
  const ScopedPointerLock& operator=(const ScopedPointerLock*);

};

//=================================================================================================

/**

This class works analogously to juce::ScopedUnlock but uses a pointer instead of a const reference
in the same way as ScopedPointerLock.

*/

class JUCE_API ScopedPointerUnlock
{
public:

  inline ScopedPointerUnlock(const CriticalSection *lock) throw() : lock_(lock)
  {
    if(lock_ != NULL)
      lock->exit();
  }

  inline ~ScopedPointerUnlock() throw()
  {
    if(lock_ != NULL)
      lock_->enter();
  }

private:

  const CriticalSection *lock_;

  ScopedPointerUnlock(const ScopedPointerLock*);
  const ScopedPointerUnlock& operator=(const ScopedPointerUnlock*);

};

#endif   
