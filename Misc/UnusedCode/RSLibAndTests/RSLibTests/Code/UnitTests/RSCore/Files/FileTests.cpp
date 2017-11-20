#include "FileTests.h"

bool testFile(std::string &reportString)
{
  std::string testName = "rsFile";
  bool testResult = true;

  testResult &= testFileTextReadWrite(reportString);
  testResult &= testFileStreamReadWrite(reportString);
  testResult &= testFileWaveReadWrite(reportString);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testFileTextReadWrite(std::string &reportString)
{
  std::string testName = "rsFileTextReadWrite";
  bool testResult = true;

  //rsFile testTextFile("d:\\TmpData\\testTextFile.txt");  // if you get an error here, make sure that the directory d:/TmpData exists
  //rsFile testTextFile("testTextFile.txt");
  //rsFile testTextFile("..\\..\\GeneratedFiles\\testTextFile.txt");
  rsFile testTextFile(rsGetCurrentApplicationDirectory() + rsString("testTextFile.txt"));

  rsString stringOriginal = createStringWithAllCharacters();
  testTextFile.writeStringToFile(stringOriginal);
  rsString stringReconstructed = testTextFile.readFileAsString();
  testResult &= ( stringOriginal == stringReconstructed );

  int fileSizeInBytes = testTextFile.getSizeInBytes();
  testResult &= ( fileSizeInBytes == 256 );

  //stringOriginal.printToStandardOutput();
  //stringReconstructed.printToStandardOutput();

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testFileStreamReadWrite(std::string &reportString)
{
  std::string testName = "rsFileStream";
  bool testResult = true;

  //rsFileStream fileStream("d:\\TmpData\\testBinaryFile.dat");
  //rsFileStream fileStream("..\\..\\GeneratedFiles\\testBinaryFile.dat");
  rsFileStream fileStream(rsGetCurrentApplicationDirectory() + rsString("testBinaryFile.txt"));

  static const int numElements = 100;

  // create and write binary data:
  double doubleArrayOriginal[numElements];
  rsFillWithIndex(doubleArrayOriginal, numElements);
  fileStream.openForWrite();
  fileStream.appendData(doubleArrayOriginal, numElements);
  fileStream.close();

  // retrieve the data again:
  double doubleArrayReconstructed[numElements];
  fileStream.openForRead();
  fileStream.readData(doubleArrayReconstructed, numElements);
  fileStream.close();

  testResult &= rsAreBuffersEqual(doubleArrayOriginal, doubleArrayReconstructed, numElements);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testFileWaveReadWrite(std::string &reportString)
{
  std::string testName = "rsFileWaveReadWrite";
  bool testResult = true;

  // create test data and write into files:
  static const int writeNumFrames  = 20000;
  int              writeSampleRate = 88200;
  int              writeBitDepth   = 16;

  rsString monoFile(  rsGetCurrentApplicationDirectory() + rsString("testWaveFileMono.wav"));
  rsString stereoFile(rsGetCurrentApplicationDirectory() + rsString("testWaveFileStereo.wav"));
  char *pathMonoFile   = monoFile.getAsZeroTerminatedString();
  char *pathStereoFile = stereoFile.getAsZeroTerminatedString();

  double *waveOriginalL  = new double[writeNumFrames];
  double *waveOriginalR  = new double[writeNumFrames];
  rsFillWithRandomValues(waveOriginalL, writeNumFrames, -1.0, 1.0, 1);
  rsFillWithRandomValues(waveOriginalR, writeNumFrames, -1.0, 1.0, 1);
  writeToMonoWaveFile(  pathMonoFile,   waveOriginalL,                writeNumFrames, writeSampleRate, writeBitDepth);
  writeToStereoWaveFile(pathStereoFile, waveOriginalL, waveOriginalR, writeNumFrames, writeSampleRate, writeBitDepth);

  // retrieve the data again from the stereo files and check if equal to original (within numerical precision limits):
  int readNumChannels, readNumFrames, readSampleRate;
  double tolerance = 1.0 / pow(2.0, (writeBitDepth-1));
  double **waveReconstructed = readFromWaveFile(pathStereoFile, readNumChannels, readNumFrames, readSampleRate);
  testResult &= (readNumChannels == 2);
  testResult &= (readNumFrames == writeNumFrames);
  if( testResult == true )
  {
    //double diff = 0.0; // only needed for inspection in the debugger when test fails
    //diff = rsMaxDeviation(waveOriginalL, waveReconstructed[0], readNumFrames);
    //diff = rsMaxDeviation(waveOriginalR, waveReconstructed[1], readNumFrames);
    testResult &= rsAreBuffersApproximatelyEqual(waveOriginalL, waveReconstructed[0], readNumFrames, tolerance);
    testResult &= rsAreBuffersApproximatelyEqual(waveOriginalR, waveReconstructed[1], readNumFrames, tolerance);
    delete[] waveReconstructed[0];
    delete[] waveReconstructed[1];
    delete[] waveReconstructed;
  }

  // now the same test with the mono-file:
  waveReconstructed = readFromWaveFile(pathMonoFile, readNumChannels, readNumFrames, readSampleRate);
  testResult &= (readNumChannels == 1);
  testResult &= (readNumFrames == writeNumFrames);
  if( testResult == true )
  {
    testResult &= rsAreBuffersApproximatelyEqual(waveOriginalL, waveReconstructed[0], readNumFrames, tolerance);
    delete[] waveReconstructed[0];
    delete[] waveReconstructed;
  }

  // memory cleanup:
  delete[] waveOriginalL;
  delete[] waveOriginalR;
  delete[] pathMonoFile;
  delete[] pathStereoFile;

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}




