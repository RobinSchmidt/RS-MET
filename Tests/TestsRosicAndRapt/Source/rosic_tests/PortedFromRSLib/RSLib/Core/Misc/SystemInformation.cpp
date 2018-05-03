using namespace RSLib;

namespace RSLib
{

rsString rsGetCurrentApplicationDirectory()
{
  rsString dir = getCurrentApplicationFile().getAbsolutePath();
  int end = dir.rsFindLastOccurrenceOf('\\');
  return dir.getSubString(0, end);
}

//=================================================================================================
// windows-specific implementation:

#if defined RS_SYSTEM_WIN32

#include <windows.h>

rsFile getCurrentApplicationFile()
{
  static const int bufferLength = MAX_PATH + 256;
  char buffer[bufferLength];  // maybe we need WCHAR to support unicode?
  GetModuleFileNameA(GetModuleHandle(0), buffer, bufferLength);
  return rsFile(rsString(buffer));
}

//=================================================================================================
// Linux-specific implementation:

#elif defined RS_SYSTEM_LINUX

#include <unistd.h>
//#include <linux/limits.h> // why this? seems, it's not needed

rsFile getCurrentApplicationFile()
{
  static const int bufferLength = PATH_MAX + 256;
  char buffer[bufferLength];
  readlink("/proc/self/exe", buffer, bufferLength);
  return rsFile(rsString(buffer));
}

//=================================================================================================
// RS_SYSTEM_OSX-specific implementation:

#elif defined RS_SYSTEM_OSX

rsFile getCurrentApplicationFile()
{
  return rsFile(); // preliminary
}

#endif

}


