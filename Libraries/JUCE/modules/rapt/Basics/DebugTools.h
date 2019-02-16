#pragma once



/** This function should be used to indicate a runtime error. */
inline void rsError(const char *errorMessage = nullptr)
{
  //printf("%s", errorMessage);
  RS_DEBUG_BREAK;
  // \todo have some conditional compilation code based on the DEBUG macro (trigger a break),
  // maybe open an error message box, etc.
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




