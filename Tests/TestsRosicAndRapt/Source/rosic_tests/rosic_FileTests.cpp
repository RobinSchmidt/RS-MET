#include "rosic_FileTests.h"
using namespace rotes;

//#include "rosic/rosic.h"
using namespace rosic;

void rotes::testFileTextReadWrite()
{
  rosic::File testTextFile("d:\\TmpData\\testTextFile.txt");
  String stringOriginal = createStringWithAllPrintableCharacters();
  testTextFile.appendText(stringOriginal);
  String stringReconstructed = testTextFile.readFileAsString();
  rassert( stringOriginal == stringReconstructed );
}


