
void rsSamplerFilter::setup(rsSamplerFilter::Type type, float cutoff, float resonance)
{

}

void rsSamplerFilter::initCoeffs()
{
  /*
  switch(type) // use switch type.topology
  {

  }
  */
}

void rsSamplerFilter::processFrame(float& L, float& R)
{

}

void rsSamplerFilter::resetState()
{

}

//=================================================================================================

/** Class where all the boilerplate for making DSP processors available in the sample goes. */

class rsSamplerProcessors
{

public:

  class Filter : public rsSamplerEngine::SignalProcessor
  {

  public:

    void processFrame(rsFloat64x2& inOut) override {}
    void processBlock(rsFloat64x2* inOut, int N) override {}
    void resetState() override { core.resetState(); }
    void resetSettings() override { core.initCoeffs(); }

  protected:

    rsSamplerFilter core;

  };


  class WaveShaper : public rsSamplerEngine::SignalProcessor
  {

  public:

    void processFrame(rsFloat64x2& inOut) override {}
    void processBlock(rsFloat64x2* inOut, int N) override {}
    void resetState() override {}
    void resetSettings() override {}

  protected:

    rsSamplerWaveShaper core;

  };

  // before writing too much boilerplate, fix the API:
  // -processFrame should work either with 2 floats or rsFloat32x4...or maybe make a class
  //  rsFloat32x2...could be just a synonym for rsFloat32x4 (because that can be simdified), but 
  //  the name makes clear that we use only 2 of the 4 available slots



};



//=================================================================================================
/*

Notes:

rsSamplerFilter:
 -Maybe make also a struct for very basic 1-pole/1-zero filter. It can be realized by the biquad 
 (and by the other structures, too), but maybe it's more efficient to do it like that. I expect 
 that a simple 1st order lowpass is a quite common thing to use.
 -Maybe
 -Use the filter also for the equalizer opcode. No need to define a different class for that. Maybe
 extend sfz to support 4 instead of 3 bands when we later can realize 2 bands per filter...



*/
