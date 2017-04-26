#ifndef rosic_MutexLock_h
#define rosic_MutexLock_h

#if (defined (_WIN32) || defined (_WIN64))
#define USE_WINAPI_CRITICAL_SECTION
#elif (defined (MAC) || defined (LINUX))
#define USE_PTHREADS_MUTEX
#endif

#if defined (USE_WINAPI_CRITICAL_SECTION)
#include <windows.h>
#elif defined (USE_PTHREADS_MUTEX)
#include <pthread.h>
#endif

namespace rosic
{

  /**

  This is a class for conveniently dealing with critical sections for classes
  that require special care that variables are not accessed concurrently by 
  different threads (for example an user-interface- and an audio-thread). This 
  is of importance for variables which must be modified only in a consistent way
  such as the coefficient set of an IIR-filter.

  credits: this class is a platform independent mutex-class derived from two 
  platform-specific implementations by kvr-fellow aciddose. it switches between 
  the two implementations by means of conditional compilation.

  */



  class MutexLock
  {

  public:

    MutexLock(char* mutexName = NULL);    
    /**< Constructor. */

    ~MutexLock();   
    /**< Destructor. */

    void lock();
    /**< Acquire the mutex-lock for some thread. */

    void unlock();
    /**< Release the mutex-lock for some thread. */

  protected:


#if defined (USE_WINAPI_CRITICAL_SECTION)

    CRITICAL_SECTION cs;

#elif defined (USE_PTHREADS_MUTEX)

    pthread_mutex_t mutex;

#endif

    char* name;

  };







} // end namespace rosic

#endif // MutexLock_h
