#include "MidiNoteEvent.h"

MidiNoteEvent::MidiNoteEvent()
{
 key       = 64;
 vel       = 64;
 detune    = 0.0;
 priority  = 0;
 timeStamp = Time();
}

MidiNoteEvent::MidiNoteEvent(int initKey, 
                             int initVel, 
                             int initDetune, 
                             int initPriority, 
                             Time initTimeStamp)
{
 if( initKey >=0 && initKey <= 127)
  key = initKey;
 else 
 {
  // throw exception...
 }
 if( initVel >=0 && initVel <= 127)
  vel = initVel;
 else 
 {
  // throw exception...
 }
 detune = initDetune;
 if( initPriority >=0 )
  priority = initPriority;
 else 
 {
  // throw exception...
 }
 timeStamp = initTimeStamp;
}

MidiNoteEvent::~MidiNoteEvent()
{

}

int MidiNoteEvent::getKey()
{
 return key;
}

void MidiNoteEvent::setKey(int newKey)
{
 if( newKey >=0 && newKey <= 127)
  key = newKey;
 else 
 {
  // throw exception...
 }
}

int MidiNoteEvent::getVel()
{
 return vel;
}

void MidiNoteEvent::setVel(int newVel)
{
 if( newVel >=0 && newVel <= 127)
  vel = newVel;
 else 
 {
  // throw exception...
 }
}

double MidiNoteEvent::getDetune()
{
 return detune;
}

void MidiNoteEvent::setDetune(double newDetune)
{
 detune = newDetune;
}

int MidiNoteEvent::getPriority()
{
 return priority;
}

void MidiNoteEvent::setPriority(int newPriority)
{
 if( newPriority >=0 )
  priority = newPriority;
 else 
 {
  // throw exception...
 }
}

Time MidiNoteEvent::getTimeStamp()
{
 return timeStamp;
}

void MidiNoteEvent::setTimeStamp(Time newTimeStamp)
{
 timeStamp = newTimeStamp;
}