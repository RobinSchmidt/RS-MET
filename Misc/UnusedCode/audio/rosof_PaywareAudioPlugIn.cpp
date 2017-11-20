#include "rosof_PaywareAudioPlugIn.h"
using namespace rosof;

PaywareAudioPlugIn::PaywareAudioPlugIn()
{
  productIndex = 0;
  isInDemoMode = 0.0;   // not a demo version by default
  demoGain     = 1.0;   // unit gain by default

  // setup the smoother for the demo-gain:
  dgSmoother.setMode(rosic::OnePoleFilter::LOWPASS);
  dgSmoother.setCutoff(10.0);

  alertLock.enter();
  alertBox = new RMessageBox();
  alertBox->setAlwaysOnTop(true);
  alertBox->centreWithSize(320, 160);
  alertBox->setOpaque(true);
  alertLock.exit();

  //initializePlugIn();  // call this at the end of the constructor of the subclass
}

PaywareAudioPlugIn::~PaywareAudioPlugIn()
{    
  alertLock.enter();
  alertBox->setVisible(false);
  delete alertBox;  
  alertLock.exit();
}

void PaywareAudioPlugIn::prepareToPlay(double sampleRate, int samplesPerBlock)
{
  // do your pre-playback setup stuff here..
  dgSmoother.setSampleRate(sampleRate);
  AudioPlugIn::prepareToPlay(sampleRate, samplesPerBlock);
}

void PaywareAudioPlugIn::initializePlugIn()
{
  checkForKeyFile();
  if( isInDemoMode != 0.0 )
  {
    plugInName = plugInName; 
    startTimer(1200000); // 1'200'000 milliseconds == 20 minutes
    //startTimer(5000); // 5 seconds for test
    if( underlyingAudioModule != NULL )
      underlyingAudioModule->setModuleNameAppendix( juce::String(T(" - Demo Version")) );
  }
}

void PaywareAudioPlugIn::checkForKeyFile()
{
  // retrieve the directory in which the plugIn resides:
  File         thisPlugInAsFile      = File::getSpecialLocation(File::currentExecutableFile);
  File         thisDirectoryAsFile   = thisPlugInAsFile.getParentDirectory();
  juce::String thisDirectoryAsString = thisDirectoryAsFile.getFullPathName();

  // search for a key file, and switch to demo-mode if no valid keyfile is found:
  juce::Array<File> foundKeyFiles;
  thisDirectoryAsFile.findChildFiles(foundKeyFiles, File::findFiles, false, plugInName+juce::String(T("Key"))+juce::String(T("*.xml")) );
  File keyFile = File::nonexistent;
  if( foundKeyFiles.size() != 0 )
  {
    keyFile      = foundKeyFiles.getFirst();
    isInDemoMode = 7.5 * !checkIfKeyFileIsValid(keyFile);
  }
  else // no keyfile is present
    isInDemoMode = 7.5;
}

bool PaywareAudioPlugIn::checkIfKeyFileIsValid(const File &keyFileToCheck)
{
  XmlDocument keyDocument(keyFileToCheck);
  XmlElement* keyXml = keyDocument.getDocumentElement();

  if( keyXml == NULL )
    return false;
  else
  {
    keyValidator.setProductIndex(productIndex);

    // read out the encoded serial and the key from the xml-element:
    juce::String encodedSerial  = keyXml->getStringAttribute(juce::String(T("Serial")), juce::String::empty);
    juce::String keyToBeChecked = keyXml->getStringAttribute(juce::String(T("Key")),    juce::String::empty);

    if( encodedSerial.length() < 8 || keyToBeChecked.length() < 64 )
    {
      showAlertBox(plugInName + juce::String(T(": Invalid Keyfile!")), 
        juce::String(T("A keyfile was found but it did not contain valid data. PlugIn will run in demo mode.")));
      delete keyXml;
      return false;
    }

    // read out the licensee's data (name and address) and pass it to the validator:
    juce::String licenseeName    = keyXml->getStringAttribute(T("Licensee"), T("No Name"));
    juce::String licenseeAddress = keyXml->getStringAttribute(T("Address"),  juce::String::empty);
    char* licenseeNameAndAddress = 
      toZeroTerminatedString(licenseeName + juce::String(T(", ")) + licenseeAddress);
    keyValidator.setLicenseeNameAndAddress(licenseeNameAndAddress);

    // convert encoded serial to chcaracter array and pass it to the validator:
    char* serialCharArray = toZeroTerminatedString(encodedSerial);
    keyValidator.setEncodedSerialNumber(serialCharArray);

    // generate the correct key-string for comparison::
    int numCharacters = keyToBeChecked.length();
    char* correctKeyCharArray = keyValidator.getKeyString(numCharacters);
    juce::String correctKey = juce::String(correctKeyCharArray);

    // free temporarily allocated memory:
    delete keyXml;
    if( licenseeNameAndAddress != NULL )
      delete[] licenseeNameAndAddress;
    if( serialCharArray != NULL )
      delete[] serialCharArray;
    if( correctKeyCharArray != NULL )
      delete[] correctKeyCharArray;

    // compare the key-string form the xml-file to the correct string:
    if( keyToBeChecked == correctKey )
    {
      return true;
    }
    else
    {
      showAlertBox(plugInName + juce::String(T(": Invalid Keyfile!")), 
        juce::String(T("A keyfile was found but the key is invalid. PlugIn will run in demo mode.")));
      return false;
    }
  }
}
    
void PaywareAudioPlugIn::processBlock(AudioSampleBuffer& buffer, MidiBuffer& midiMessages)
{
  AudioPlugIn::processBlock(buffer, midiMessages);

  if( isInDemoMode != 0.0 )
  {
    float *left  = buffer.getSampleData(0, 0);
    float *right = buffer.getSampleData(1, 0);
    float  g;
    for(int n=0; n<buffer.getNumSamples(); n++)
    {
      g = (float) dgSmoother.getSample(demoGain);
      left[n]  *= g;
      right[n] *= g;
    }
  }
}

//==============================================================================

void PaywareAudioPlugIn::getStateInformation(MemoryBlock& destData)
{
  if( isInDemoMode != 0.0 )
  {
  /*
  // this alert-window crashes reaper when one plugs in straightliner, adds an effect-plugin and 
  // switches back and forth between straightliner and the effect
    AlertWindow::showMessageBox(AlertWindow::WarningIcon, 
      plugInName + juce::String(T(": Total Recall disabled!")), 
      juce::String(T("The host is trying to save the state of the plugin for later total recall but this is disabled in this demo version. The plugin-state will be saved, but only the full version will allow to restore it. If you have just created the uber-killer-sound, you can also save it via the plugin's own preset management.")), 
      juce::String(T("OK")));
  */
  }
  XmlElement* xmlState = getStateAsXml();
  copyXmlToBinary(*xmlState, destData);
  delete xmlState;
}

void PaywareAudioPlugIn::setStateInformation (const void* data, int sizeInBytes)
{
  if( isInDemoMode != 0.0 )
  {
    showAlertBox(plugInName + juce::String(T(": Total Recall disabled!")), 
      juce::String(T("The host is trying recall the state of the plugin but total recall is disabled in this demo version. Only the full version will allow to restore the saved state.")));
  }
  else
  {
    AudioPlugIn::setStateInformation(data, sizeInBytes);
    /*
    XmlElement* const xmlState = getXmlFromBinary (data, sizeInBytes);
    ParameterObserver::globalAutomationSwitch = false; // why this - threading problems?
    setStateFromXml(*xmlState, juce::String(T("recalled by host")));
    ParameterObserver::globalAutomationSwitch = true;
    //underlyingAudioModule->updateCoreObjectAccordingToParameters(); // this messes up total recall in straightliner
    //if( underlyingAudioModule->underlyingRosicInstrument != NULL )
    //  underlyingAudioModule->underlyingRosicInstrument->setPresetName("recalled by host", true);
    delete xmlState;
    */
  }
}

void PaywareAudioPlugIn::timerCallback()
{
  if( isInDemoMode != 0.0 )
  {
    demoGain = 0.0;
    stopTimer();
    showAlertBox(plugInName + juce::String(T(": Demo timed out!")),
      juce::String(T("The 20 minute demo period has expired. The plugIn will be muted now until it is reloaded.")));
  }
}

bool PaywareAudioPlugIn::isDemoVersion()
{
  return (isInDemoMode != 0.0);
}

void PaywareAudioPlugIn::showAlertBox(const juce::String &headlineText, const juce::String &bodyText)
{
  alertLock.enter();
  alertBox->setHeadlineText(headlineText);
  alertBox->setBodyText(bodyText);
  alertBox->addToDesktop(ComponentPeer::windowIsTemporary);
  alertBox->setVisible(true);
  alertLock.exit();
}