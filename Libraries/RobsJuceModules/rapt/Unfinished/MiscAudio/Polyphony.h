#pragma once


// move to rosic

/** Baseclass for all DSP classes that need to implement polyphony handling... */


class rsVoice
{

public:

  virtual void noteOn(int key, int vel) = 0;

  virtual void noteOff(int key, int vel) = 0;

  //virtual void setPitchBend(double amount) {}

  //virtual bool hasFinished();

protected:

  // ...it would probably be best, if the baseclass does not have any data members at all

};

//=================================================================================================

/** For recursively composing more complex voice objects from simpler ones - for example, a 
synthesizer voice could be composed of two oscillator voices, a filter voice, envelope generator 
voices etc. */

class rsComplexVoice : public rsVoice
{

public:

  // addChildVoice, removeChildVoice, ...

protected:

  std::vector<rsVoice> childVoices;

};

//=================================================================================================

class rsVoiceManager
{

public:



};