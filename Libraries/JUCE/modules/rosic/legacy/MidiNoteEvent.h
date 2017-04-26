#ifndef MidiNoteEvent_h
#define MidiNoteEvent_h


#include "Time.h"

/**

This is a class for MIDI note-events.

*/

class MidiNoteEvent  
{

public:
 
 //---------------------------------------------------------------------------
 // construction/destruction:

	MidiNoteEvent();          ///< Default constructor.
 MidiNoteEvent(int  initKey, 
               int  initVel, 
               int  initDetune = 0.0, 
               int  initPriority = 0, 
               Time initTimeStamp = Time() );
 /**< Constructor with initializations. */

	virtual ~MidiNoteEvent(); ///< Destructor.

 //---------------------------------------------------------------------------
 // public access-functions:

 int    getKey();
 void   setKey(int newKey);
 int    getVel();
 void   setVel(int newVel);
 double getDetune();
 void   setDetune(double newDetune);
 int    getPriority();
 void   setPriority(int newPriority);
 Time   getTimeStamp();
 void   setTimeStamp(Time newTimeStamp);

 //---------------------------------------------------------------------------
 // overloaded operators:

 bool operator==(const MidiNoteEvent& note2) const  
 {
  if( note2.key == key )
   return true;
  else
   return false;
 }
 /**< Note events are interpreted as equal if the have the same key. */

protected:

 // data members:
 int    key;       // key of the note in the range 0...127
 int    vel;       // velocity of the note in the range 0...127
 double detune;    // detuning in cents (for microtuning) 
 int    priority;  // a priority
 Time   timeStamp; // time of occurence

};

#endif // MidiNoteEvent_h
