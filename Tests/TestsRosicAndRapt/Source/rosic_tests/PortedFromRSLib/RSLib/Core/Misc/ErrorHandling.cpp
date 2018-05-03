namespace RSLib
{

  void rsError(const char *errorMessage)
  {
    printf("%s", errorMessage);
    RS_DEBUG_BREAK;
      // \todo have some conditional compilation code based on the DEBUG macro (trigger a break),
      // maybe open an error message box, etc.
  }

  void rsAssert(bool expression, const char *errorMessage)
  {
    if( expression == false )
      rsError(errorMessage);
  }

}
