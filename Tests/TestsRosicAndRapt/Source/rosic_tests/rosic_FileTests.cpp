#include "rosic_FileTests.h"
using namespace rotes;

void rotes::testFileTextReadWrite()
{
  rosic::File testTextFile("d:\\TmpData\\testTextFile.txt");
  String stringOriginal = createStringWithAllPrintableCharacters();
  testTextFile.appendText(stringOriginal);
  String stringReconstructed = testTextFile.readFileAsString();
  rassert( stringOriginal == stringReconstructed );
}


