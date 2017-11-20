#include "rosof_ChannelMatrix2x2AudioModule.h"
using namespace rosof;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

ChannelMatrix2x2AudioModule::ChannelMatrix2x2AudioModule(CriticalSection *newPlugInLock, rosic::ChannelMatrix2x2 *channelMatrix2x2ToWrap)
 : AudioModule(newPlugInLock)
{
  jassert(channelMatrix2x2ToWrap != NULL); // you must pass a valid rosic-object to the constructor
  wrappedChannelMatrix2x2 = channelMatrix2x2ToWrap;
  moduleName              = juce::String(T("ChannelMatrix2x2"));
  setActiveDirectory(getApplicationDirectory() + juce::String(T("/ChannelMatrix2x2Presets")) );
  initializeAutomatableParameters();
}

//-------------------------------------------------------------------------------------------------
// state management:

XmlElement* ChannelMatrix2x2AudioModule::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  // store the inherited controller mappings:
  XmlElement *xmlState = AudioModule::getStateAsXml(stateName, markAsClean);

  if( wrappedChannelMatrix2x2 == NULL )
    return xmlState;

  // add attributes for the non-automatable parameters (the automatable ones are already taken care
  // of by AudioModule::getStateAsXml()):
  /*
  xmlState->setAttribute(T("M11"), wrappedChannelMatrix2x2->getM11());
  xmlState->setAttribute(T("M12"), wrappedChannelMatrix2x2->getM12());
  xmlState->setAttribute(T("M21"), wrappedChannelMatrix2x2->getM21());
  xmlState->setAttribute(T("M22"), wrappedChannelMatrix2x2->getM22());
  */
  xmlState->setAttribute(T("LeftToLeft"),   wrappedChannelMatrix2x2->getLeftToLeftGain());
  xmlState->setAttribute(T("RightToLeft"),  wrappedChannelMatrix2x2->getRightToLeftGain());
  xmlState->setAttribute(T("LeftToRight"),  wrappedChannelMatrix2x2->getLeftToRightGain());
  xmlState->setAttribute(T("RightToRight"), wrappedChannelMatrix2x2->getRightToRightGain());

  return xmlState;
}

void ChannelMatrix2x2AudioModule::setStateFromXml(const XmlElement& xmlState,
                                                 const juce::String& stateName, bool markAsClean)
{
  // restore the settings of the inherited AudioModule object:
  AudioModule::setStateFromXml(xmlState, stateName, markAsClean);

  if( wrappedChannelMatrix2x2 == NULL )
    return;

  // restore the values of the non-automatable parameters:
  wrappedChannelMatrix2x2->setLeftToLeftGain(  xmlState.getDoubleAttribute(T("M11"), 1.0));
  wrappedChannelMatrix2x2->setRightToLeftGain( xmlState.getDoubleAttribute(T("M12"), 0.0));
  wrappedChannelMatrix2x2->setLeftToRightGain( xmlState.getDoubleAttribute(T("M21"), 0.0));
  wrappedChannelMatrix2x2->setRightToRightGain(xmlState.getDoubleAttribute(T("M22"), 1.0));

  // these M-attributes are only a legacy from the old version, nowadays we use more intuitive
  // attribute names:
  wrappedChannelMatrix2x2->setLeftToLeftGain(  xmlState.getDoubleAttribute(T("LeftToLeft"), 1.0));
  wrappedChannelMatrix2x2->setRightToLeftGain( xmlState.getDoubleAttribute(T("RightToLeft"), 0.0));
  wrappedChannelMatrix2x2->setLeftToRightGain( xmlState.getDoubleAttribute(T("LeftToRight"), 0.0));
  wrappedChannelMatrix2x2->setRightToRightGain(xmlState.getDoubleAttribute(T("RightToRight"), 1.0));
}