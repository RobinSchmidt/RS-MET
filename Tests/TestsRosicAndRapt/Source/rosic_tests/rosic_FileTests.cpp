#include "rosic_FileTests.h"
using namespace rotes;

//#include "rosic/rosic.h"
using namespace rosic;

void rotes::testFileTextReadWrite()
{
  rosic::File testTextFile("d:\\TmpData\\testTextFile.txt");
  rsString stringOriginal = createStringWithAllPrintableCharacters();
  testTextFile.appendText(stringOriginal);
  rsString stringReconstructed = testTextFile.readFileAsString();
  rassert( stringOriginal == stringReconstructed );
}


