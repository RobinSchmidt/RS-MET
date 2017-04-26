#include "rosic_MutexLock.h"
using namespace rosic;

MutexLock::MutexLock(char* mutexName)
{
#if defined (USE_WINAPI_CRITICAL_SECTION)

  InitializeCriticalSection(&cs);

#elif defined (USE_PTHREADS_MUTEX)

  pthread_mutex_init(&mutex, NULL);

#endif

  name = NULL;
  if( mutexName != NULL )
  {
    int length = strlen(mutexName);
    name = new char[length+1];
    strcpy(name, mutexName);
  }
}

MutexLock::~MutexLock()
{
#if defined (USE_WINAPI_CRITICAL_SECTION)

  DeleteCriticalSection(&cs);

#elif defined (USE_PTHREADS_MUTEX)

  pthread_mutex_destroy(&mutex);

#endif

  delete[] name;
}

void MutexLock::lock()
{
#if defined (USE_WINAPI_CRITICAL_SECTION)

  EnterCriticalSection(&cs);

#elif defined (USE_PTHREADS_MUTEX)

  pthread_mutex_lock(&mutex);

#endif
}

void MutexLock::unlock()
{
#if defined (USE_WINAPI_CRITICAL_SECTION)

  LeaveCriticalSection(&cs);

#elif defined (USE_PTHREADS_MUTEX)

  pthread_mutex_unlock(&mutex);

#endif
}



