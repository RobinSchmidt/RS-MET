#include "../../../../RS-MET/Misc/APE Scripts/rapt_for_ape.cpp"
// relative path from "Audio Programming Environment/includes"

#include <effect.h> // or maybe we should use generator.h?


//=================================================================================================
// The core DSP object

template<class T>
class rsShepardToneGenerator
{
  
public:

  //-----------------------------------------------------------------------------------------------
  // \name Setup

   /** Sets the sample rate at which this object runs. This determines the time increment for our
   phasor */
   inline void setSampleRate(T newRate) { dt = T(1) / newRate; }

   /** Sets the lower und upper cutoff frequency between which all the action happens. */
   inline void setCutoffFrequencies(T lowCutoff, T highCutoff)
   {
     wLo = T(2) * PI * lowCutoff  * dt;
     wHi = T(2) * PI * highCutoff * dt;
   }

  //-----------------------------------------------------------------------------------------------
  // \name Processing

   /** Computes the desired gain factor for a given radian frequency w. */
   inline T getGainForOmega(T w)
   {
     return T(1); // preliminary - todo: implement a bell curve here
   }

  /** Produces one output sample at a time. */
  inline T getSample()
  {
    T y = T(0);   // output

    // Create all frequencies from wRef up to the upper cutoff:
    T w = wRef;
    while(w <= wHi)
    {
      y += getGainForOmega(w) * sin(w*t);
      w *= T(2);
    }

    // Create all frequencies from 0.5*wRef down to the lower cutoff:
    w = T(0.5)*wRef;
    while(w >= wLo)
    {
      y += getGainForOmega(w) * sin(w*t);
      w *= T(0.5);
    }

    // Increment time and return output:
    t += dt;
    return y;
  }

  inline void reset()
  {
    t = T(0);
  }
  
protected:

  T t;         // current time, normalized to 0..1. our phasor
  T dt;        // time increment per sample for t. equal to 1/sampleRate
  T wRef;      // current reference omega
  T wLo, wHi;  // high and low cutoff omegas

};


//=================================================================================================
// The APE effect


GlobalData(ShepardTone, "Endless glissando generator");

class ShepardTone : public ape::Effect
{
public:


  ShepardTone() {}
  
  // Shorthands for convenience to reduce boilerplate:
  using Par = ape::Param<float>;
  using Rng = ape::Range;
  using Map = Rng::Mapping;

  Par parGain{   "Gain",     Rng(-48,  12)              }; // in dB
  Par parFreqLo{ "FreqLo",   Rng(20,   20000, Map::Exp) }; // lower cutoff in Hz
  Par parFreqHi{ "FreqHi",   Rng(20,   20000, Map::Exp) }; // upper cutoff in Hz

  Par parFreq{   "Freq",     Rng(20,   20000, Map::Exp) }; // current reference freq
  // this is preliminary - later this is the freq that should automatically increase or decrease
  // and wrap around when it reaches to times its original value



private:   

  rsShepardToneGenerator<float> core;

  /** Resets the internal state. */
  void reset()
  {
    core.reset();
  }

  //-----------------------------------------------------------------------------------------------
  // \name Overriden callbacks (why are they private?) 

  void start(const ape::IOConfig& cfg) override
  {
    reset();
  }

  void process(ape::umatrix<const float> ins, ape::umatrix<float> outs, size_t numFrames) override
  {
    // Pull the values of the parameters for this block and update the DSP object:
    const float gain = RAPT::rsDbToAmp((float)parGain);
    float s  = 2*PI / config().sampleRate;
    float wl = s * (float)parFreqLo;
    float wh = s * (float)parFreqHi;
    //core.setOmegas(wl, wh);

    // Loop over the sample frames:
    const auto numChannels = sharedChannels();     // shared? with whom?
    for(size_t n = 0; n < numFrames; ++n)
    {
      float y = gain * core.getSample();      // compute output
      for(size_t c = 0; c < numChannels; ++c) // write generated sample into all channels
        outs[c][n] = y;                 
    }

    clear(outs, numChannels); // what does this do? clear unused channels?
  }
};

