/* 

MIDIEventInterpreter.h: interface for the MIDIEventInterpreter class.

© Braindoc 2006 (www.braindoc.de)



*/

#if !defined(MIDIEventInterpreter_h_Included)
#define MIDIEventInterpreter_h_Included

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "VstTools.h"

//enumreate the types of MIDI-events:
enum
{
 kNoteOn,
 kNoteOff,
 kControlChange,
 kPitchWheel,
 kPatchChange,
 kSystemExclusive,

 kOther,
};

class MIDIEventInterpreter  
{
public:

 //construction/destruction
	MIDIEventInterpreter();
	virtual ~MIDIEventInterpreter();

 //evaluate event:
 void interpretEvent(uint8 midiBytes[3], uint8 *type, uint8 *channel, 
                     uint8 *data1, uint8 *data2);

 //convert pitch-wheel data bytes into a value bewtween -1 and +1:
 sample interpretPitchWheel(uint8 HighByte, uint8 LowByte);

};

#endif // !defined(MIDIEventInterpreter_h_Included)
