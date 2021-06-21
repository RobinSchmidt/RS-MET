using namespace rotes;
using namespace rosic;

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


  return ok;
}

