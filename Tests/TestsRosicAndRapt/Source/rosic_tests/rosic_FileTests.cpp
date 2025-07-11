using namespace rotes;
using namespace rosic;

//#include "../RSLib/Core/RSCore.h"  // get rid
#include "PortedFromRSLib/RSLib/Core/RSCore.h"  // get rid - it should go to rosic

namespace rotes
{

rsString createStringWithAllCharacters()
{
  char cString[256];
  for(int i=0; i<256; i++)
    cString[256-i-1] = i;  // backwards for compliance with C strings
  return rsString(cString);
}
rsString createStringWithAllPrintableCharacters()
{
  static const int firstPrintableIndex = 32;   // whitespace ' '
  static const int lastPrintableIndex  = 127;  // tilde '~'
  static const int numPrintables       = lastPrintableIndex-firstPrintableIndex;

  char cString[numPrintables+1];
  for(int i=0; i<numPrintables; i++)
    cString[i] = i+firstPrintableIndex;
  cString[numPrintables] = '\0';
  return rsString(cString);
}

bool testFileText()
{
  // Tests, if we can write a string into a file and retrieve it again, the string must not 
  // contain non-printable characters.

  bool ok = true;

  // Create a test string:
  rsString stringOriginal = createStringWithAllPrintableCharacters();

  // Write the string to a file and read it back again - the two strings must match. We use 
  // rsFile::appendText and rsFile::readFileAsString for this.
  //rosic::rsFile testTextFile("E:\\TmpData\\testTextFile.txt"); // absolute path
  rosic::rsFile testTextFile("testTextFile.txt"); // path relative from current directory
  testTextFile.appendText(stringOriginal);
  rsString stringReconstructed = testTextFile.readFileAsString();
  ok &= stringOriginal == stringReconstructed;

  // Write and read the string again. This time we use the convenience functions 
  // rsWriteStringToFile and rsReadStringFromFile which operate on C-style strings.
  rosic::rsWriteStringToFile("testTextFile2.txt", stringOriginal.getRawString());
  char* stringReconstructed2 = rosic::rsReadStringFromFile("testTextFile2.txt");
  int cmp = strcmp(stringOriginal.getRawString(), stringReconstructed2);
  ok &= cmp == 0;
  free(stringReconstructed2);

  return ok;

  // ToDo:
  // 
  // 
  // - Write testTextFile.txt and testTextFile2.txt into a subdirectory TempFiles or something like
  //   that and then add that directory to the .gitignore file. This is better than directly adding
  //   the files to .gitignore. BUT: I'm not sure if it works. Check, if we are correctly handling
  //   the path separators on all platforms.  // 
  // 
  // - Use a subdirectory using relative paths
  // 
  // - Use paths for which subdirectory does and does not exist beforehand
  // 
  // - Use different path seperators (i.e. forward slash insetad of backslash)
  // 
  // - Test rsWriteStringToFile, rsReadStringFromFile
}

bool testFileWave()
{
  bool ok = true;

  using WF = RSLib::rsWaveFile;

  // Test roundtrip of 16 bit numbers through conversion to 32 bit float:
  float f;
  int   i;
  f = WF::int16ToFloat32(-32768);    // maps to -1.00004578
  f = WF::int16ToFloat32(+32767);    // maps to +0.999984741
  i = WF::float32ToInt16(-1.f);      // maps to -32767
  i = WF::float32ToInt16(+1.f);      // maps to +32767
  for(i = -32768; i <= 32767; i++)
  {
    f = WF::int16ToFloat32((rsInt16)i);
    rsInt16 j = WF::float32ToInt16(f);
    ok &= i == j;
  }

  // The same for 24 bit numbers:
  f = WF::int24ToFloat32(-8388608);  // maps to -1.00000012
  f = WF::int24ToFloat32(+8388607);  // maps to +1.00000000
  i = WF::float32ToInt24(-1.f);      // maps to -8388607
  i = WF::float32ToInt24(+1.f);      // maps to +8388607
  for(i = -8388608; i <= 8388607; i++)
  {
    f = WF::int24ToFloat32(i);
    rsInt32 j = WF::float32ToInt24(f);
    ok &= i == j;
  }

  return ok;

  // ToDo:
  // -For the sake of completeness, implement and test 8 bit conversion and maybe also 32-bit int
  //  ...but for that, exhaustive testing my take too long...although, 4 billion may still be
  //  managable
  // -test actually writing and reading files in various formats (different bit-depths, 
  //  sample-rates, numbers of channels, etc.)
  // -maybe use constants min16 = -32768, max16 = +32767, etc.
  // -Test the roundtrip in parentheses in: int -> (float -> int -> float) to check if the final
  //  float matches the first.
}

}