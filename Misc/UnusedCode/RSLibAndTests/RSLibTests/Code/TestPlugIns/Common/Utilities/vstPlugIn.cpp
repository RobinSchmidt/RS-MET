#include "vstPlugIn.h"

vstPlugIn::vstPlugIn(audioMasterCallback audioMaster, VstInt32 numPrograms, VstInt32 numParams) 
  : AudioEffectX(audioMaster, numPrograms, numParams)
{
  setNumInputs(2);
  setNumOutputs(2);
  canProcessReplacing();
  canDoubleReplacing();

  programsAreChunks(false); // check, if this is needed
  vst_strncpy(programName, "Default", kVstMaxProgNameLen);
  params  = new float[numParams];

  bufSize = 0;
  buf[0]  = nullptr;
  buf[1]  = nullptr;
}

vstPlugIn::~vstPlugIn()
{
  delete[] params;
  freeBuffers();
}

void vstPlugIn::processReplacing(float** inputs, float** outputs, VstInt32 sampleFrames)
{
  //mutex.lock();

  // enlarge buffer, if necessary:
  if(sampleFrames > bufSize)
    allocateBuffers(sampleFrames);

  // convert input to double precision, store in internal buffers:
  for(int n = 0; n < sampleFrames; n++)
  {
    buf[0][n] = (double) inputs[0][n];
    buf[1][n] = (double) inputs[1][n];
  }

  // process the internal double precision buffer:
  processDoubleReplacing(buf, buf, sampleFrames);

  // convert back to single precision:
  for(int n = 0; n < sampleFrames; n++)
  {
    outputs[0][n] = (float) buf[0][n];
    outputs[1][n] = (float) buf[1][n];
  }

  //mutex.unlock();
}
void vstPlugIn::processDoubleReplacing(double** inputs, double** outputs, VstInt32 sampleFrames)
{
  //mutex.lock();
  for(int n = 0; n < sampleFrames; n++)
    processStereoFrame(&inputs[0][n], &inputs[1][n], &outputs[0][n], &outputs[1][n]);
  //mutex.unlock();
}

void vstPlugIn::setParameter(VstInt32 index, float value)
{
  params[index] = value;
  //mutex.lock();
  updateCoreParameter(index, value);
  //mutex.unlock();
}

float vstPlugIn::getParameter(VstInt32 index)
{
  return params[index];
}

void vstPlugIn::setProgramName (char* name)
{
  vst_strncpy(programName, name, kVstMaxProgNameLen);
}
void vstPlugIn::getProgramName (char* name)
{
  vst_strncpy(name, programName, kVstMaxProgNameLen);
}

bool vstPlugIn::getVendorString (char* text)
{
  vst_strncpy(text, "RS-MET", kVstMaxVendorStrLen);
  return true;
}

bool vstPlugIn::getProductString(char* text)
{
  return getEffectName(text);
}

VstInt32 vstPlugIn::getVendorVersion ()
{ 
  return 1000; 
}

void vstPlugIn::allocateBuffers(int size)
{
  if(size == bufSize)
    return;  // nothing to do
  freeBuffers();
  buf[0]  = new double[size];
  buf[1]  = new double[size];
  bufSize = size;
}

void vstPlugIn::freeBuffers()
{
  delete[] buf[0];
  delete[] buf[1];
  buf[0]  = nullptr;
  buf[1]  = nullptr;
  bufSize = 0;
}

//=================================================================================================

vstInstrument::vstInstrument(audioMasterCallback audioMaster, VstInt32 numPrograms, 
  VstInt32 numParams)
  : vstPlugIn(audioMaster, numPrograms, numParams)
{
  //setNumInputs(0);  // maybe not necessary - we can process input, too
  isSynth();
}

vstInstrument::~vstInstrument()
{

}

void vstInstrument::handleEvent(VstMidiEvent midiEvent)
{
  char* midiData = midiEvent.midiData;

  // converts all messages to corresponding messages on MIDI channel 1:
  long status = midiData[0] & 0xf0;

  // respond to note-on and note-off on channel 1:
  if( status == 0x90 || status == 0x80 )
  {
    // bitwise AND with 0x7f=01111111 sets the first bit of an 8-bit word to zero, allowing
    // only numbers from 0 to 127 to pass unchanged, higher numbers (128-255) will be mapped
    // into this range by subtracting 128:
    long note     = midiData[1] & 0x7f;
    long velocity = midiData[2] & 0x7f;
    long detune   = midiEvent.detune;

    // respond to note-off on channel 1 (status: 0x80)
    if( status == 0x80 )
      velocity = 0;	// note off by note-off message

    // respond to note-on (note-offs are handled there too by recognizing zero velocity):
    onNoteOn(note, velocity, detune);
  }

  // respond to all notes off:
  // (status=0xb0: controller on ch 1, midiData[1]=0x7b: control 123: all notes off):
  else if( (status) == 0xb0 && (midiData[1] == 0x7b) )
  {
    for(int i=0; i <= 127; i++)
      onNoteOn(i, 0, 0);
  }
  else if( status == 0xb0 )  // respond to controllers
    onControlChange((int)midiData[1], (int)midiData[2]);
  else if (status == 0xe0)  // respond to pitchbend
  {
    unsigned short bendValue;
    bendValue   = (unsigned short) midiData[2];
    bendValue <<= 7;
    bendValue  |= (unsigned short) midiData[1]; // between 0x0000 and 0x3FFF, center: 0x2000(=8192)
    onPitchWheel(bendValue - 8192);             // between -8192 and +8191
  }
}

inline void copyVstEvent(VstEvent *src, VstEvent *dst)
{
  dst->byteSize    = src->byteSize;
  dst->deltaFrames = src->deltaFrames;
  dst->flags       = src->flags;
  dst->type        = src->type;
  for(int j = 0; j < 16; j++)
    dst->data[j]   = src->data[j];
}
VstInt32 vstInstrument::processEvents(VstEvents* ev)
{
  if(ev->numEvents == 0)
    return 1;

  // collect the events into our internal buffer:
  events.resize(ev->numEvents);
  for(int i = 0; i < ev->numEvents; i++)
    copyVstEvent(ev->events[i], &events[i]);

  return 1;	 
  // return value 1 indicates that we want to continue to receive events
}

void vstInstrument::processDoubleReplacing(double** inputs, double** outputs, 
  VstInt32 sampleFrames)
{
  // call baseclass method, if there are no events to handle:
  if(events.size() == 0)
    vstPlugIn::processDoubleReplacing(inputs, outputs, sampleFrames);
  else
  {
    // loop for the case when there are events to be considered:
    int eventCounter = 0;
    for(int n = 0; n < sampleFrames; n++)
    {
      // check, if at this sample one or more events have occurred:
      while( (eventCounter < events.size()) && (n == events[eventCounter].deltaFrames) )
      {
        // yes, an event has occurred now -> process the event, if it's a MIDI message:
        if( events[eventCounter].type == kVstMidiType )
        {
          // cast VstEvent to VstMidiEvent with a pointer trick:
          VstEvent     *pVstEvent  = &(events[eventCounter]);
          VstMidiEvent *pMidiEvent =  (VstMidiEvent*)(pVstEvent);
          VstMidiEvent midiEvent   =  *pMidiEvent;
          handleEvent(midiEvent);
          eventCounter++;
        }
        else // the event was not of kVstMidiType and is ignored:
          eventCounter++;
      }

      // all events at this sample instant have been handled - now render the audio sample:
      processStereoFrame(&inputs[0][n], &inputs[1][n], &outputs[0][n], &outputs[1][n]);
    }
  }

  events.clear(); // buffer was used up and is invalid now
}

VstInt32 vstInstrument::canDo(char* text)
{
  if (!strcmp (text, "receiveVstEvents"))
    return 1;
  if (!strcmp (text, "receiveVstMidiEvent"))
    return 1;
  return -1;	// explicitly can't do; 0 => don't know
}
