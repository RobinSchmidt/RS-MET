#ifndef jura_MidiMessageFilter_h
#define jura_MidiMessageFilter_h

/** This class can be used to define a filter for MidiMessages by setting it up via a bunch of 
flags that define certain types of messages as suitable or not and then later letting is check 
specific messages on their suitability. */

class MidiMessageFilter
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  MidiMessageFilter();

  /** Destructor. */
  virtual ~MidiMessageFilter();

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Informs if this message is suitable according to the filter settings of this object. */
  virtual bool isMessageSuitable(const MidiMessage& m);

  //-----------------------------------------------------------------------------------------------
  // public data members:

  // Flags to indicate whether the respective types of messages are classified as suitable or not:
  bool noteOn, noteOff, sysEx, programChange, pitchWheel, aftertouch, channelPressure, controller, 
    metaEvent, activeSense, transport, clock, songPosition, machineControl, other;

  juce_UseDebuggingNewOperator;
};


#endif	
