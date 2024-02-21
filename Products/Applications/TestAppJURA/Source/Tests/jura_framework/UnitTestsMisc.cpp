#include "UnitTestsMisc.h"
using namespace juce;
using namespace jura;

UnitTestMisc::UnitTestMisc() 
  : juce::UnitTest("Misc", "Misc")
{

}

void UnitTestMisc::runTest()
{
  runTestColor();
  runTestStateFileManager();
}

void UnitTestMisc::runTestColor()
{
  beginTest("Color");

  using HSL = jura::ColourAHSL;
  using RGB = juce::Colour;

  // Tests, if the given hsl values convert correctly to the given target rgb values:
  auto testHSL2RGB = [&](float h, float s, float l, uint8 r, uint8 g, uint8 b)
  {
    HSL hsl(h, s, l);
    RGB rgb = hsl.getAsJuceColour();
    expect(rgb == RGB(r, g, b));
  };

  // Create HSL colors with S=0.5, L=0.5 and H going from 0 to 1 in 0.1 steps and check the result
  // of converting to RGB as juce::Color object. The target RGB colors in the following tests were 
  // obtained by running the code on a windows machine on which the conversion is known to work:
  testHSL2RGB(0.0f, 0.5f, 0.5f, 191,  63,  63);
  testHSL2RGB(0.1f, 0.5f, 0.5f, 191, 140,  63);
  testHSL2RGB(0.2f, 0.5f, 0.5f, 165, 191,  63);
  testHSL2RGB(0.3f, 0.5f, 0.5f,  89, 191,  63);
  testHSL2RGB(0.4f, 0.5f, 0.5f,  63, 191, 114);
  testHSL2RGB(0.5f, 0.5f, 0.5f,  63, 191, 191);
  testHSL2RGB(0.6f, 0.5f, 0.5f,  63, 114, 191);
  testHSL2RGB(0.7f, 0.5f, 0.5f,  89,  63, 191);
  testHSL2RGB(0.8f, 0.5f, 0.5f, 165,  63, 191);
  testHSL2RGB(0.9f, 0.5f, 0.5f, 191,  63, 140);
  testHSL2RGB(1.0f, 0.5f, 0.5f, 191,  63,  63);

  // OK - this test passes on mac too. I guess, the wrong colors are caused by a different 
  // ordering of r,g,b on may (maybe r and b swapped) and in our bilinear gradient rendering 
  // function we assume an r,g,b ordering. ToDo: write a function that figures out what ordering is 
  // used on the given machine by using juce::Color and setting it to full red (or green or blue) 
  // and retrieving the 32-bit integer that represents the RGBA value. Instead of directly using
  // +0, +1, +2 for the indexing of the channels, use +R, +G, +B where R,G,B are some permutation 
  // of 0,1,2 determined by the color ordering of the machine

  //RGB red   = RGB::fromRGBA(255,   0,   0, 255);
  //RGB green = RGB::fromRGBA(  0, 255,   0, 255);
  //RGB blue  = RGB::fromRGBA(  0,   0, 255, 255);


  // Figure out the in-memory ordering of the channels:
  using ARGB = juce::PixelARGB;
  PixelARGB fullAlpha(255,   0,   0,   0);  // full alpha, all colors zero
  PixelARGB fullRed(    0, 255,   0,   0);  // full red, everything else zero (including alpha)
  PixelARGB fullGreen(  0,   0, 255,   0);
  PixelARGB fullBlue(   0,   0,   0, 255);  



  // test setPixelRGB in jura_GraphicsTools.h
  // see PixelARGB::getInARGBMemoryOrder


}


// A class for testing jura::StateFileManager. We cannot instantiate the class directly because it
// is abstract, so we need this dummy subclass:
class JUCE_API TestStateFileManager : public jura::StateFileManager
{
  
public:

  void setStateFromXml(const juce::XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean) override
  {
    // Do nothing
  }

  juce::XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override
  {
    return nullptr;
  }

};
void UnitTestMisc::runTestStateFileManager()
{
  // Create a couple of example .xml files in some special subdirectory of the executable's 
  // directory, if they don't exist already:
  juce::File   exeFile = juce::File::getSpecialLocation(File::currentApplicationFile);
  juce::File   exeDir  = exeFile.getParentDirectory();
  juce::String str     = exeDir.getFullPathName();
  str += "/StateFileManagerTest";
  juce::File tmpDir(str);
  if(!tmpDir.exists())
  {
    bool ok = tmpDir.createDirectory();
    if(!ok)
      showWarningBox("Warning", "Temporary directory for tests of class jura::StateFileManager could not be created");
  }
  int numFiles = 5;
  for(int i = 1; i <= numFiles; i++)
  {
    juce::String fileName = str + "/TestFile" + juce::String(i) + ".xml";
    juce::File   tmpFile  = juce::File(fileName);
    if(!tmpFile.exists())
    {
      bool ok = tmpFile.create();
      if(!ok)
        showWarningBox("Warning", "Temporary file for tests of class jura::StateFileManager could not be created");
    }
  }

  // Now it should be ensured that in tmpDir there are numFiles empty .xml files. These will now be
  // managed by a (dummy-subclass of) jura::StateFileManager:
  TestStateFileManager mngr;
  bool ok = true;
  ok &= mngr.setRootDirectory(exeDir, false);
  ok &= mngr.setActiveDirectory(tmpDir); 
  ok &= mngr.getNumFilesInList() == numFiles;
  ok &= mngr.getActiveFile() == File();        // No active file yet.

  // Get a list of all files in the tmpDir:
  juce::Array<File> fileList;
  juce::File(mngr.getActiveDirectory()).findChildFiles(fileList, File::findFiles, false);
  ok &= fileList.size() == numFiles;
  // The files in this list should be only the .xml files that were created above.



  expect(ok);





  int dummy = 0;

  // ToDo: 
  // -Create also a special kind of FileManagerListener that logs the callbacks and check that they
  //  are called correctly. Mayb do the same for StateWatcher
}