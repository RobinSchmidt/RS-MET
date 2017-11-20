#ifndef rosof_ChannelMatrix2x2AudioModule_h
#define rosof_ChannelMatrix2x2AudioModule_h

#include "../../../rosic/basics/rosic_ChannelMatrix2x2.h"
using namespace rosic;

#include "../rosof_AudioModule.h"

namespace rosof
{

  /**

  This class wraps rosic::ChannelMatrix2x2 into a rosof::AudioModule to facilitate its use as 
  plugIn.

  */

  class ChannelMatrix2x2AudioModule : public AudioModule
  {

    friend class ChannelMatrix2x2ModuleEditor;

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    ChannelMatrix2x2AudioModule(CriticalSection *newPlugInLock, rosic::ChannelMatrix2x2 *channelMatrix2x2ToWrap);   

    //---------------------------------------------------------------------------------------------
    // automation and state management:

    virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
      bool markAsClean);

    virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);

    virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
    { wrappedChannelMatrix2x2->getSampleFrameStereo(inOutL, inOutR); }

    //=============================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    /** Pointer to the underlying rosic object which is wrapped. */
    rosic::ChannelMatrix2x2 *wrappedChannelMatrix2x2;

  };

}

#endif 
