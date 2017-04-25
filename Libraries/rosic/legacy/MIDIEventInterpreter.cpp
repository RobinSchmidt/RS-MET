// MIDIEventInterpreter.cpp: implementation of the MIDIEventInterpreter class.
//
//////////////////////////////////////////////////////////////////////

#include "MIDIEventInterpreter.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

MIDIEventInterpreter::MIDIEventInterpreter()
{

}

MIDIEventInterpreter::~MIDIEventInterpreter()
{

}

//-----------------------------------------------------------------------------
//event processing
void MIDIEventInterpreter::interpretEvent(uint8 midiBytes[3], uint8 *type, 
                                          uint8 *channel, uint8 *data1, uint8 *data2)
{




}

sample MIDIEventInterpreter::interpretPitchWheel(uint8 HighByte, uint8 LowByte)
{
 static uint16 bendValue;  //value of the wheel (between 0x0000 and 0x3FFF, center: 0x2000(=8192))
 static double pitchValue; //mapped value (into +-1)

 //bendValue   = 0x0000;
 bendValue   = (uint16) LowByte;
 bendValue <<= 7;
 bendValue  |= (uint16) HighByte;
 //bendValue now holds a value between 0x0000 and 0x3FFF 
 //where 0x2000 means: pitch wheel is centered

 pitchValue = (double)bendValue - 8192;   //pitchValue is between -8192 and +8191
 if(pitchValue == -8192)
  pitchValue = -8191;                     //for symmetry
 pitchValue = (pitchValue/8191);          //pitchValue is between -1 and +1

 return pitchValue;
}

