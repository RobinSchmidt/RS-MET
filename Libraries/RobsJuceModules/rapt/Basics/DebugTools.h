#pragma once

// move to somewhere else:
template<class T>
bool doNothing(T x)
{
  return x == x;
}

/** This function should be used to indicate a runtime error. */
inline void rsError(const char *errorMessage = nullptr)
{
  doNothing(errorMessage);
  // fixes "unreferenced formal parameter" warning - we don't do anything with the message, but 
  // we may want to see it in the debugger

  //printf("%s", errorMessage);
  RS_DEBUG_BREAK;
  // \todo have some conditional compilation code based on the DEBUG macro (trigger a break),
  // maybe open an error message box, etc.
}

inline void rsWarning(const char* message = nullptr)
{
  std::cout << "Warning: " << message << "\n";
}


/** This function should be used for runtime assertions. */
inline void rsAssert(bool expression, const char *errorMessage = nullptr)
{
  if( expression == false )
    rsError(errorMessage);
}

inline void rsAssertFalse(const char *errorMessage = nullptr) { rsAssert(false, errorMessage); }
// lol! this is stupid! remove and use rsError directly!

// todo: have functions for logging and printing




