#ifndef RS_SYSTEMINFORMATION_H
#define RS_SYSTEMINFORMATION_H

/** This file contains functions to retrieve informations about the current machine, operating 
system etc.

\todo maybe wrap these functions into a class "Platform" and make them static member functions.

*/

namespace RSLib
{

  /** Returns the currently running application as a file. */
  RSLib_API rsFile getCurrentApplicationFile();
    // prepend rs-prefix

  /** Returns the directory in which this application sits as string (including the final 
  backslash). */
  RSLib_API rsString rsGetCurrentApplicationDirectory();

  // \todo File getStandardApplicationForExtension(const String& fileExtension); 
  // String getCurrentSystemTime(); 
  // int getCurrentThreadID();
  // getCurrentProcessID();

}

#endif
