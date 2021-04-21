#include "../../../../RS-MET/Misc/APE Scripts/rapt_for_ape.cpp"
// relative path from "Audio Programming Environment/includes"

#include <effect.h>

using namespace ape;

GlobalData(LadderFilter, "Moog-style Ladder Filter");

class LadderFilter : public ape::Effect
{
  
public:

  using LDR = RAPT::rsLadderFilter<float, float>;

  LadderFilter() {}
  

  static constexpr typename Param<LDR::Mode>::Names ModeNames = 
  { 
    "Bypass",
    "Lowpass_6",
    "Lowpass_12",
    "Lowpass_18",
    "Lowpass_24",
    "Highpass_6",
    "Highpass_12",
    "Highpass_18",
    "Highpass_24",
    "Bandpass_6_6",
    "Bandpass_6_12",
    "Bandpass_6_18",
    "Bandpass_12_6",
    "Bandpass_12_12",
    "Bandpass_18_6",
  };
  

  // Shorthands for convenience to reduce boilerplate:
  using Par = ape::Param<float>;
  using Rng = ape::Range;
  using Map = Rng::Mapping;

  Par parFreq { "Freq",  Rng(20,  20000, Map::Exp) }; // in Hz
  Par parReso { "Reso",  Rng(0,   400)             }; // in %
  Par parDrive{ "Drive", Rng(-12, 48)              }; // in dB
  Par parGain { "Gain",  Rng(-48, 12)              }; // in dB
  Par parB1   { "B1",    Rng(0,   0.5)             }; // 0: allpole, 0.5: bilinear
  
  // ToDo: Mode, Drive (pre-gain), Gain (post-gain), DC/Asym, SatShape
  
  
  Param<LDR::Mode> parMode { "Mode", ModeNames };
  
  //Param<bool> parBilinear { "Bilinear" };
  
  
  // How can we set the parameters to reasonable default values?


private:

  LDR flt;
  float sampleRate;

  
  /** Resets the internal state. */
  void reset()
  {
    flt.reset();
  }

  //-----------------------------------------------------------------------------------------------
  // \name Overriden callbacks (why are they private?) 

  void start(const ape::IOConfig& cfg) override
  {
    sampleRate = (float)cfg.sampleRate;
    flt.setSampleRate(sampleRate);
    reset();
  }

  void process(ape::umatrix<const float> ins, ape::umatrix<float> outs, size_t numFrames) override
  {
    // Pull the values of the parameters for this block and update the DSP objects:
    const float freq  = (float)parFreq;
    const float reso  = (float)parReso * 0.01f;             // convert from percent to raw
    const float drive = RAPT::rsDbToAmp((float)parDrive);
    const float gain  = RAPT::rsDbToAmp((float)parGain);
    flt.setCutoff(freq);
    flt.setResonance(reso);
    flt.setMode(parMode);
    flt.setB1(parB1);
    //flt.setBilinear(parBilinear);
    // todo: set them per sample
    
    // Loop over the sample frames:
    const auto numChannels = sharedChannels();
    for(size_t n = 0; n < numFrames; ++n)
    {
      float x = drive * ins[0][n];            // read input - we just use the 1st channel
      float y = gain  * flt.getSample(x);     // compute output
      for(size_t c = 0; c < numChannels; ++c) // write generated sample into all channels
        outs[c][n] = y;                 
    }

    clear(outs, numChannels); // what does this do? clear unused channels?
  }
};

// Notes:
// -with high cutoff (> 12 kHz) and 100% resonance, it starts behaving chatocially
//  -> oversample by factor 2 (factor out a DSP class that wraps that)
// -with even more resonance (> 100%), the artifacts begin to show already at lower frequencies
// -When a filter should support resonances > 100, we need a nonlinear parameter mapping, 
//  maybe 100% should occur at 2/3 or 3/4 of the full range. actually, that would be good anyway
//  because from 0..50%, not much seems to be happening and before 100, a lot is happening. be 
//  sure to test this with settings where the drive does not diminish the resonance
// -it seems like the artifacts can be combatted by using smoother nonlinearities
// -when the resonance is above 100, the drive can be used to diminish the resonance - the signal
//  gets smashed aginst the limits...however, in a sawtooth, the resonace recovers in the middle
//  where the signal values are low
// -the bilinear version behaves differently with respect to these artifacts. it seems to tend to
//  alias more. regular seems to produce more chaos/noise. but that could be because the 
//  resonance gets less diminished?
// -by adjusting the drive just right (during self-oscillations), growl can be achieved
// -introduce an envelope that can be used cut off the self-oscillation when the input is quiet
//  -> use an env-follower with instantaneous attack and adjustable decay on the input
//  -> 2 parameters: env-amount, env-decay
//  -> maybe it should have a hold-time, too
//  ...because it can be annyoing when the filter constantly produces sound even without input








