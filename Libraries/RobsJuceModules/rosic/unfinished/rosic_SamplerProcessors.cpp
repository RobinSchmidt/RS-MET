
void rsSamplerFilter::setup(rsSamplerFilter::Type type, float cutoff, float resonance)
{


  int dummy = 0;
}

void rsSamplerFilter::processFrame(float& L, float& R)
{


  int dummy = 0;
}


/*

rsSamplerFilter:
-Maybe make also a struct for very basic 1-pole/1-zero filter. It can be realized by the biquad 
 (and by the other structures, too), but maybe it's more efficient to do it like that. I expect 
 that a simple 1st order lowpass is a quite common thing to use.
-Maybe
-Use the filter also for the equalizer opcode. No need to define a different class for that. Maybe
 extend sfz to support 4 instead of 3 bands when we later can realize 2 bands per filter...

*/


//=================================================================================================

/** Class where all the boilerplate for making DSP processors available in the sample goes. */

class rsSamplerProcessors
{

public:

  class Filter : public rsSamplerEngine::SignalProcessor
  {

  public:

    void processFrame(rsFloat64x2& inOut) override
    {

    }

    void processBlock(rsFloat64x2* inOut, int N) override
    {

    }

    void resetState() override
    {

    }

    void resetSettings() override
    {

    }

  protected:

    rsSamplerFilter core;

  };

};