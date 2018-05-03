#ifndef RS_ERRORHANDLING_H
#define RS_ERRORHANDLING_H

namespace RSLib
{

  /** This function should be used to indicate a runtime error. */
  RSLib_API void rsError(const char *errorMessage);

  /** This function should be used for runtime assertions. */
  RSLib_API void rsAssert(bool expression, const char *errorMessage = NULL);

}

#endif
