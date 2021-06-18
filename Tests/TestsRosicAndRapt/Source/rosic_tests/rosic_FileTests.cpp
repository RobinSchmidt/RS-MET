//#include "rosic_FileTests.h"
using namespace rotes;

//#include "rosic/rosic.h"
using namespace rosic;

bool rotes::testFileTextReadWrite()
{
  bool ok = true;
  rosic::rsFile testTextFile("E:\\TmpData\\testTextFile.txt"); // absolute path
  rsString stringOriginal = createStringWithAllPrintableCharacters();
  testTextFile.appendText(stringOriginal);
  rsString stringReconstructed = testTextFile.readFileAsString();
  ok &= stringOriginal == stringReconstructed;
  return ok;

  // Todo: 
  // -use relative path 
  // -use paths where subdirectory does and does not exist
  // -use different path seperators (i.e. forward slash insetad of backslash)
}


