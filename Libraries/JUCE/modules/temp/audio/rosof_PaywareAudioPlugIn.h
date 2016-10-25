#ifndef rosof_PaywareAudioPlugIn_h
#define rosof_PaywareAudioPlugIn_h

#include "../../rosic/filters/rosic_OnePoleFilter.h"
#include "../../rosic/others/rosic_KeyGenerator.h"
using namespace rosic;

#include "rosof_AudioPlugIn.h"

#include "../../rojue/components/widgets/rojue_RMessageBox.h"

namespace rosof
{

  class PaywareAudioPlugIn : public AudioPlugIn, public Timer
  {

  public:

    PaywareAudioPlugIn();
    virtual ~PaywareAudioPlugIn();

    virtual void prepareToPlay(double sampleRate, int samplesPerBlock);

    virtual void getStateInformation(MemoryBlock& destData);

    virtual void setStateInformation(const void* data, int sizeInBytes);

    /** Overriden to multiply the block with the 'demo-gain' factor. */
    virtual void processBlock(AudioSampleBuffer& buffer, MidiBuffer& midiMessages);

    /** This is the callback which triggers the gain to be set to zero after 30 minutes for demo 
    versions. */
    virtual void timerCallback();

    /** Informs whether this is a demo version of not. */
    virtual bool isDemoVersion();

    //===============================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    /** This method is supposed to be used for initialization which have to be done on loading the 
    plugIn such as checking whether this is a demo version and imposing the required 
    demo-restrictions (initialize the timer, etc.). Make sure to override it in your subclass and 
    call it in the constructor of the subclass. Can be ignored for freeware plugins. */
    virtual void initializePlugIn();

    /** This function checks for the existence and validity of a keyfile and returns true if a valid 
    keyfile is there and false otherwise. */
    virtual void checkForKeyFile();

    /** Check if a keyfile is valid or not. */
    virtual bool checkIfKeyFileIsValid(const juce::File& keyFileToCheck);

    /** This must be overriden in subclasses to calculate one stereo sample frame - the slots serve 
    as inputs and outputs at the same time. */
    //virtual void getSampleFrameStereo(double *left, double *right);

    /** Shows an alert box. */
    virtual void showAlertBox(const juce::String& headlineText, const juce::String& bodyText);


    int    productIndex; // index for the product - needed to check for the keyfile
    double isInDemoMode; // flag to indicate whether this is a demo-version or not - we use a 
                         // double to obfuscate the meaning for crackers

    //double isKeyFileProtected;
    /**< Flag to indicate whether this is a commercial product proteccted by a keyfile or not - we 
    use a double to obfuscate the meaning for crackers. */

    double demoGain;
    /**< 'Demo-Gain' - a gain factor to multiply the signal with. The factor will be 1.0 for plugIns 
    which are either freeware, licensed payware or the 30 minute demo period is still running, 0.0 
    otherwise (demo-versions running longer than 30 minutes). */

    rosic::OnePoleFilter dgSmoother;
    /**< A Smoother for the 'Demo-Gain' to avoid clicks after the 30 minute demo period. */

    rosic::KeyGenerator keyValidator;
    /**< The object which will be used to verify the validity of the key. */

    CriticalSection alertLock;
    RMessageBox*    alertBox;

  };

}

#endif  
