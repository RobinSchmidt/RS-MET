using namespace rotes;
using namespace rosic;

//#include "../RSLib/Core/RSCore.h"  // get rid
#include "PortedFromRSLib/RSLib/Core/RSCore.h"  // get rid - it should go to rosic

bool rotes::testFileText()
{
  // Tests, if we can write a string into a file and retrieve it again, the string must not 
  // contain non-printable characters..

  bool ok = true;
  rosic::rsFile testTextFile("E:\\TmpData\\testTextFile.txt"); // absolute path
  rsString stringOriginal = createStringWithAllPrintableCharacters();
  testTextFile.appendText(stringOriginal);
  rsString stringReconstructed = testTextFile.readFileAsString();
  ok &= stringOriginal == stringReconstructed;
  return ok;

  // ToDo: 
  // -use relative path 
  // -use paths where subdirectory does and does not exist
  // -use different path seperators (i.e. forward slash insetad of backslash)
}

bool rotes::testFileWave()
{
  bool ok = true;

  using WF = RSLib::rsWaveFile;

  // Test roundtrip of 16 bit numbers through conversion to 32 bit float:
  float f;
  f = WF::int16ToFloat32(-32768);  // maps to -1.00004578
  f = WF::int16ToFloat32(+32767);  // maps to +0.999984741
  for(rsInt16 i = -32768; i <= 32767; i++)
  {
    f = WF::int16ToFloat32(i);
    rsInt16 j = WF::float32ToInt16(f);
    ok &= i == j;
    if(i == 32767)
      break;          // without that, we enter an infinite loop
  }

  // ToDo:
  // -do the same test for all possible 24bit values


  return ok;
}

