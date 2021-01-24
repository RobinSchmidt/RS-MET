#ifndef jura_VoiceManager_h
#define jura_VoiceManager_h


/** A class for managing polyphony in instrument modules. */

class rsVoiceManager
{

public:

  rsVoiceManager()
  {
    setMaxNumVoices(16);
  }


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets the maximum number of voices that should be supported. The function is supposed to be 
  called once shortly after construction and then that setting should remain fixed for the lifetime
  of the plugin. If we later want to allow the user to change that setting at runtime, we will need 
  facilities to trigger a re-allocation of all the required resources (buffer, DSP-objects, 
  etc.). We don't have such facilities yet and maybe that feature is not worth the increased 
  complexity anyway...we'll see... */
  void setMaxNumVoices(int newNumber)
  {
    maxNumVoices    = newNumber;
    numVoices       = RAPT::rsMin(numVoices,       maxNumVoices);
    numActiveVoices = RAPT::rsMin(numActiveVoices, maxNumVoices);
    playingVoices.resize(maxNumVoices);
  }

  /** Sets the number of voices that should be available for playing. */
  void setNumVoices(int newNumber)
  {
    numVoices = RAPT::rsMin(newNumber, maxNumVoices);

    //numActiveVoices = RAPT::rsMin(numActiveVoices, numVoices); 
    // hmmm...might be a bit harsh to kill off voices immediately - maybe we should wait for them 
    // to be released naturally
  }


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  int getMaxNumVoices() const { return maxNumVoices; }

  int getNumVoices() const { return numVoices; }

  int getNumActiveVoices() const { return numActiveVoices; } 


protected:

  int maxNumVoices    = 16;  // maximum number of voices
  int numVoices       =  8;  // number of available voices
  int numActiveVoices =  0;  // number of currently playing voices

  std::vector<int> playingVoices;

};








/** A thin wrapper around rosic::rsVoiceManager to provide the glue to the juce framework, 
specifically, the translation of juce::MidiMessage. */

/*
class JUCE_API rsVoiceManager : public rosic::rsVoiceManager
{

public:


protected:


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsVoiceManager)
};
*/

#endif