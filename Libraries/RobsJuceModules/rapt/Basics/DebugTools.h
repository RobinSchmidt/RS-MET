#pragma once

// move to somewhere else:
template<class T>
bool doNothing(T x)
{
  return x == x;
}

inline void rsPrintLine(const std::string& message)
{
  std::string str = message + "\n";
  //std::cout << str;        // APE (Audio Programming Environment) does not support std::cout...
  printf("%s", str.c_str()); // ...so we resort to the old-skool C way of doing it
}

/** This function should be used to indicate a runtime error. */
inline void rsError(const char *message = nullptr)
{
#ifdef RS_DEBUG
  //doNothing(errorMessage);
  // fixes "unreferenced formal parameter" warning - we don't do anything with the message, but 
  // we may want to see it in the debugger 

  if(message) rsPrintLine("Error " + std::string(message));
  else        rsPrintLine("Error");
  RS_DEBUG_BREAK;
  // \todo have some conditional compilation code based on the DEBUG macro (trigger a break),
  // maybe open an error message box, etc.
#endif
}
// -Maybe have an else branch that writes to some error log file and/or maybe throws and exception
// -Maybe have a stronger version rsFatalError that also triggers some action in release builds
// -Maybe call the error function that only affects debug builds rsDebugError



inline void rsWarning(const char* message = nullptr)
{
#ifdef RS_DEBUG
  if(message) rsPrintLine("Warning " + std::string(message));
  else        rsPrintLine("Warning");
#endif
}

/** This function should be used for runtime assertions. */
inline void rsAssert(bool expression, const char *errorMessage = nullptr)
{
#ifdef RS_DEBUG
  if( expression == false )
    rsError(errorMessage);
#endif
}
// Maybe we should have another version that fires an error in release-builds, too? But maybe the
// way the error is fired should be different? Maybe open an error message box? ...but that 
// requires to drag in system libraries, which i'd rather not here. like so:
// http://www.cplusplus.com/forum/beginner/53571/
// Maybe call this rsAssertDbg and the other rsAssertRls?
// Maybe we should throw an exception in release builds? see here:
// https://docs.microsoft.com/en-us/cpp/cpp/errors-and-exception-handling-modern-cpp?view=msvc-160

inline void rsAssertFalse(const char *errorMessage = nullptr) { rsAssert(false, errorMessage); }
// LOL! This is stupid! We should use rsError directly instead! JUCE had (or still has?) such a 
// thing and I just mimicked the idea from there.

// todo: have functions for logging and printing

inline void rsStaticAssert(bool expression, const char* errorMessage = nullptr)
{
#ifdef RS_CPP17
  static_assert(expression, errorMessage);
#else
  rsAssert(expression, errorMessage);
#endif
}


/*
Some Notes on hard to catch bugs

Heap Corruptions:

I think, heap corruptions may occur erratically when trying to write beyond the end of a 
std::vector when the function, that accesses the vector, does so via a pointer to its first 
element. In this case, the range-check in std::vector itself is bypassed, so we don't trigger the 
debug assertion there. Whether or not the heap actually gets corrupted depends on how much memory 
the vector has allocated which may be more than its current size. The heap corruption only occurs
when we try to write beyond its capacity, which we have no control over, because the growing is 
controlled by the standard library - which is why it's so erratic and hard to reproduce reliably 
(...i think...although, the capacity should actually be deterministic...hmmm - but i sometimes get
or don't those heap corruptions depending on some change in a totally unrelated piece of code)

*/
