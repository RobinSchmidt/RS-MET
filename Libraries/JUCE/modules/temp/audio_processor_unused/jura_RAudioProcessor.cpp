//#include "RAudioProcessor.h"
//#include "RAudioProcessorEditor.h"

RAudioProcessor::RAudioProcessor()
{
  //plugInEngine = NULL; // from old codebase

  productIndex       = 0;
  isKeyFileProtected = 0.0;   // freebie by default
  isInDemoMode       = 0.0;   // not a demo version by default
  demoGain           = 1.0;   // unit gain by default

  // from old codebase:
  //// setup the smoother for the demo-gain:
  //dgSmoother.setMode(rosic::OnePoleFilter::LOWPASS);
  //dgSmoother.setCutoff(10.0);

  //initializePlugIn();
}

RAudioProcessor::~RAudioProcessor()
{
  // from old codebase:
  //if( plugInEngine != NULL )
  //{
  //  delete plugInEngine;
  //  plugInEngine = NULL;
  //}
}

const String RAudioProcessor::getName() const
{
  return plugInName;
}

int RAudioProcessor::getNumParameters()
{
  return 0;
}

float RAudioProcessor::getParameter(int index)
{
  return 0.f;
}

void RAudioProcessor::setParameter(int index, float newValue)
{

}

const String RAudioProcessor::getParameterName(int index)
{
  return String::empty;
}

const String RAudioProcessor::getParameterText(int index)
{
  return String::empty;
}

//const String RAudioProcessor::getInputChannelName(const int channelIndex) const
//{
//  if( channelIndex == 0 )
//    return String("Left Input");
//  else if( channelIndex == 1 )
//    return String("Right Input");
//  else
//    return String("Input Channel ") + String(channelIndex+1);
//}
//
//const String RAudioProcessor::getOutputChannelName(const int channelIndex) const
//{
//  if( channelIndex == 0 )
//    return String("Left Output");
//  else if( channelIndex == 1 )
//    return String("Right Output");
//  else
//    return String("Output Channel ") + String(channelIndex+1);
//}
//
//bool RAudioProcessor::isInputChannelStereoPair(int index) const
//{
//  return true;
//}
//
//bool RAudioProcessor::isOutputChannelStereoPair(int index) const
//{
//  return true;
//}

bool RAudioProcessor::acceptsMidi() const
{
  return false;
}

bool RAudioProcessor::producesMidi() const
{
  return false;
}

void RAudioProcessor::prepareToPlay(double sampleRate, int samplesPerBlock)
{
  // from old codebase:
  //// do your pre-playback setup stuff here..
  //dgSmoother.setSampleRate(sampleRate);
  //if( plugInEngine != NULL )
  //  plugInEngine->setSampleRate(sampleRate);
}

void RAudioProcessor::releaseResources()
{
  // when playback stops, you can use this as an opportunity to free up any
  // spare memory, etc.
}

void RAudioProcessor::processBlock(AudioSampleBuffer& buffer, MidiBuffer& midiMessages)
{
  // make sure that the buffer has the right number of channels:
  if( buffer.getNumChannels() != 2 )
    return;

  // request time-info from host and update the bpm-values for the modulators accordingly (if they 
  // are in sync mode):
  AudioPlayHead::CurrentPositionInfo info;
  if( getPlayHead() != 0  &&  getPlayHead()->getCurrentPosition(info) )
    setBeatsPerMinute(info.bpm); 

  // some stuff for the input midi-buffer
  MidiBuffer::Iterator midiBufferIterator(midiMessages);
  MidiMessage currentMidiMessage(0x80,0,0,0.0);
  int currentMidiMessageOffset;
  bool aMidiMessageWasRetrieved = false;

  // process the buffers sample-by-sample:
  float *left  = buffer.getWritePointer(0);
  float *right = buffer.getWritePointer(1);
  for(int i = 0; i < buffer.getNumSamples(); i++)
  {
    // respond to midi-messages at this instant of time:
    midiBufferIterator.setNextSamplePosition(i);
    currentMidiMessageOffset = -1; 
    aMidiMessageWasRetrieved = midiBufferIterator.getNextEvent(
      currentMidiMessage, currentMidiMessageOffset);
    while( aMidiMessageWasRetrieved && currentMidiMessageOffset == i )
    {
      respondToMidiMessage(currentMidiMessage);
      aMidiMessageWasRetrieved = midiBufferIterator.getNextEvent(
        currentMidiMessage, currentMidiMessageOffset);
    }

    // calculate the output samples:
    double doubleLeft  = (double) (left[i]);
    double doubleRight = (double) (right[i]);
    getSampleFrameStereo(&doubleLeft, &doubleRight);
    left[i]  = (float) doubleLeft;
    right[i] = (float) doubleRight;
  }

  //midiMessages.clear();
}
//void RAudioProcessor::processBlock(AudioSampleBuffer& buffer, MidiBuffer& midiMessages)
//{
//  // we need to make sure that the pointer to the audio-engine is valid because it will be 
//  // de-referenced inside th audio-loop:
//  if( plugInEngine == NULL )
//  {
//    buffer.clear();
//    return;
//  }
//
//  // make sure that the buffer has the right number of channels:
//  if( buffer.getNumChannels() != 2 )
//    return;
//
//  // request time-info from host and update the bpm-values for the modulators accordingly (if they 
//  // are in sync mode):
//  AudioPlayHead::CurrentPositionInfo info;
//  if( getPlayHead() != 0  &&  getPlayHead()->getCurrentPosition(info) )
//    setBeatsPerMinute(info.bpm); 
//
//  // some stuff for the input midi-buffer
//  MidiBuffer::Iterator midiBufferIterator(midiMessages);
//  MidiMessage          currentMidiMessage(0x80,0,0,0.0);
//  int                  currentMidiMessageOffset;
//  bool                 aMidiMessageWasRetrieved = false;
//
//  float *floatSamplePointerLeft, *floatSamplePointerRight;
//  // pointers to the current input- and output-samples
//
//  double doubleSampleLeft, doubleSampleRight;  
//  // the current input- and output-samples, typecasted to double
//
//  // process the buffers sample-by-sample:
//  for(int i=0; i<buffer.getNumSamples(); i++)
//  {
//    floatSamplePointerLeft  = buffer.getSampleData(0, i);
//    floatSamplePointerRight = buffer.getSampleData(1, i);
//
//    // typecast the inputs to double values to make them compatible with the
//    // audio engine:
//    doubleSampleLeft  = (double) (*floatSamplePointerLeft);
//    doubleSampleRight = (double) (*floatSamplePointerRight);
//
//    // check if there are midi-messages at this instant of time:
//    midiBufferIterator.setNextSamplePosition(i);
//    currentMidiMessageOffset = -1; 
//    // initialization - will be set to another value by the function call 
//    // in the next line...
//
//    aMidiMessageWasRetrieved = midiBufferIterator.getNextEvent(
//      currentMidiMessage, currentMidiMessageOffset);
//
//    while( aMidiMessageWasRetrieved && currentMidiMessageOffset == i )
//    {
//      // respond to the midi-message:
//      respondToMidiMessage(currentMidiMessage);
//
//      // get the next message at this sample:
//      aMidiMessageWasRetrieved = midiBufferIterator.getNextEvent(
//        currentMidiMessage, currentMidiMessageOffset);
//    }
//
//    // calculate the output samples:
//    getSampleFrameStereo(&doubleSampleLeft, &doubleSampleRight);
//
//    // typecast the calculated output-samples to float and store them into
//    // their dedicated memory locations:
//    *floatSamplePointerLeft  = (float) doubleSampleLeft;
//    *floatSamplePointerRight = (float) doubleSampleRight;
//  }
//
//  //midiMessages.clear();
//}

void RAudioProcessor::initializePlugIn()
{
  checkForKeyFile();
  if( isInDemoMode != 0.0 )
  {
    //startTimer(10000);      // 10'000 milliseconds == 10 second for test purposes
    startTimer(1200000); // 1'200'000 milliseconds == 20 minutes
  }
}

void RAudioProcessor::checkForKeyFile()
{
  isInDemoMode = 0;

  // code temporarily deactivated

  //// retrieve the directory in which the plugIn resides:
  //File   thisPlugInAsFile      = File::getSpecialLocation(File::currentExecutableFile);
  //File   thisDirectoryAsFile   = thisPlugInAsFile.getParentDirectory();
  //String thisDirectoryAsString = thisDirectoryAsFile.getFullPathName();

  //// search for a key file, and switch to demo-mode if no valid keyfile is found:
  //OwnedArray<File> foundKeyFiles;
  ////thisDirectoryAsFile.findChildFiles(foundKeyFiles, File::findFiles, 
  ////  false, JUCE_T("MagicCarpetKey*.xml") );
  //thisDirectoryAsFile.findChildFiles(foundKeyFiles, File::findFiles, 
  //  false, plugInName + String("*.xml") );
  //File keyFile = File::nonexistent;
  //if( foundKeyFiles.size() != 0 )
  //{
  //  keyFile      = *foundKeyFiles.getFirst();
  //  isInDemoMode = 7.5 * !checkIfKeyFileIsValid(keyFile);
  //}
  //else // no keyfile is present
  //{
  //  isInDemoMode = 7.5;
  //}
}

bool RAudioProcessor::checkIfKeyFileIsValid(const File &keyFileToCheck)
{
  return true;

  // code temporarily deactivated

  //XmlDocument keyDocument(keyFileToCheck);
  //XmlElement* keyXml = keyDocument.getDocumentElement();

  //if( keyXml == NULL )
  //  return false;
  //else
  //{
  //  keyValidator.setProductIndex(productIndex);

  //  // read out the encoded serial and the key from the xml-element:
  //  String encodedSerial  = keyXml->getStringAttribute(String("Serial"), 0);
  //  String keyToBeChecked = keyXml->getStringAttribute(String("Key"),    0);

  //  if( encodedSerial.length() < 8 || keyToBeChecked.length() < 64 )
  //  {
  //    AlertWindow::showMessageBox(AlertWindow::WarningIcon, 
  //      plugInName + String(": Invalid Keyfile!"), 
  //      String("A keyfile was found but it did not contain valid data. PlugIn will run in demo mode."), 
  //      String("OK"));
  //    return false;
  //  }

  //  // read out the licensee's data (name and address) and pass it to the validator:
  //  String licenseeName    = keyXml->getStringAttribute("Licensee", "No Name");
  //  String licenseeAddress = keyXml->getStringAttribute("Address",  String::empty);
  //  char* licenseeNameAndAddress = 
  //    toZeroTerminatedString(licenseeName + String(", ") + licenseeAddress);
  //  keyValidator.setLicenseeNameAndAddress(licenseeNameAndAddress);

  //  // convert encoded serial to chcaracter array and pass it to the validator:
  //  char* serialCharArray = toZeroTerminatedString(encodedSerial);
  //  keyValidator.setEncodedSerialNumber(serialCharArray);

  //  // generate the correct key-string for comparison::
  //  int numCharacters = keyToBeChecked.length();
  //  char* correctKeyCharArray = keyValidator.getKeyString(numCharacters);
  //  String correctKey = String(correctKeyCharArray);

  //  // free temporarily allocated memory:
  //  if( licenseeNameAndAddress != NULL )
  //    delete[] licenseeNameAndAddress;
  //  if( serialCharArray != NULL )
  //    delete[] serialCharArray;
  //  if( correctKeyCharArray != NULL )
  //    delete[] correctKeyCharArray;

  //  // compare the key-string form the xml-file to the correct string:
  //  if( keyToBeChecked == correctKey )
  //  {
  //    return true;
  //  }
  //  else
  //  {
  //    AlertWindow::showMessageBox(AlertWindow::WarningIcon, 
  //      plugInName + String(": Invalid Keyfile!"), 
  //      String("A keyfile was found but the key is invalid. PlugIn will run in demo mode."), 
  //      String("OK"));
  //    return false;
  //  }
  //}
}

void RAudioProcessor::getSampleFrameStereo(double *left, double *right)
{
  plugInEngine->getSampleFrameStereo(left, right, left, right);
  if( isInDemoMode != 0.0 )
  {
    double g = dgSmoother.getSample(demoGain);
    *left  *= g;
    *right *= g;
  }
}

void RAudioProcessor::respondToMidiMessage(MidiMessage theMessage)
{
  if( plugInEngine == NULL )
    return;

  if( theMessage.isController() )
  {
    int controllerNumber = theMessage.getControllerNumber();
    int controllerValue  = theMessage.getControllerValue();
    plugInEngine->setMidiController(controllerNumber, controllerValue);
  }
  //if( theMessage.isNoteOn() )
  //  plugInEngine->noteOn(theMessage.getNoteNumber(), theMessage.getVelocity() );
  //else if( theMessage.isNoteOff() )
  //  synth.noteOn(theMessage.getNoteNumber(), 0);
  //else if( theMessage.isAllNotesOff() )
  //  synth.allNotesOff();
  //else if( theMessage.isPitchWheel() )
  //{
  //  int    wheelValue       = theMessage.getPitchWheelValue();
  //  double wheelValueMapped = (double) (wheelValue-8192) / 8192.0; // check this
  //  synth.setPitchBend(wheelValueMapped);
  //}
}

void RAudioProcessor::setBeatsPerMinute(double newBpm)
{
  if( plugInEngine != NULL )
    plugInEngine->setBeatsPerMinute(newBpm);
}

AudioProcessorEditor* RAudioProcessor::createEditor()
{
  // its not a good idea to update slider during the construction of the GUI, so we temporarily 
  // deactivate automation:
  AutomationListener::guiAutomationSwitch = false;

  AudioProcessorEditor* editor = new RAudioProcessorEditor(this);

  // switch automation for GUI object on again and return the editor:
  AutomationListener::guiAutomationSwitch = true;
  return editor;
}

void RAudioProcessor::getStateInformation(JUCE_NAMESPACE::MemoryBlock& destData)
{
  if( isInDemoMode != 0.0 )
  {
    AlertWindow::showMessageBox(AlertWindow::WarningIcon, 
      plugInName + String(": Total Recall disabled!"), 
      String("The host is trying to save the state of the plugin for later total recall but this is disabled in this demo version. The plugin-state will be saved, but only the full version will allow to restore it. If you have just created the uber-killer-sound, you can also save it via the plugin's own preset management."), 
      String("OK"));
  }
  XmlElement* xmlState = getStateAsXml();
  copyXmlToBinary(*xmlState, destData);
  delete xmlState;
}

void RAudioProcessor::setStateInformation(const void* data, int sizeInBytes)
{
  if( isInDemoMode != 0.0 )
  {
    AlertWindow::showMessageBox(AlertWindow::WarningIcon, 
      plugInName + String(": Total Recall disabled!"), 
      String("The host is trying recall the state of the plugin but total recall is disabled in this demo version. Only the full version will allow to restore the saved state."), 
      String("OK"));
  }
  else
  {
    XmlElement* const xmlState = getXmlFromBinary (data, sizeInBytes);
    AutomationListener::globalAutomationSwitch = false;
    setStateFromXml(*xmlState);
    AutomationListener::globalAutomationSwitch = true;
    plugInEngine->setPresetName("recalled by host", true);
    delete xmlState;
  }
}

/*
XmlElement* RAudioProcessor::getStateAsXml()
{
  return new XmlElement(String(T("RAudioProcessorState")));
  //return straightlinerStateToXml(&synth);
}

void RAudioProcessor::setStateFromXml(const XmlElement &xmlState)
{
  //straightlinerStateFromXml(&synth, xmlState);
  //synth.resetAllVoices();
  //synth.markPresetAsClean();
}
*/

void RAudioProcessor::timerCallback()
{
  if( isInDemoMode != 0.0 )
  {
    demoGain = 0.0;
    AlertWindow::showMessageBox(AlertWindow::WarningIcon, 
      plugInName + String(": Demo timed out!"), 
      String("The 20 minute demo period has expired. The plugIn will be muted now until it is 
        reloaded."), 
      String("OK"));
    stopTimer();
  }
}

bool RAudioProcessor::isDemoVersion()
{
  return (isInDemoMode != 0.0);
}