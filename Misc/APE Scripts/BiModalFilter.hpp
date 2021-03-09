#include "../../../../RS-MET/Misc/APE Scripts/rapt_for_ape.cpp"
// relative path from "Audio Programming Environment/includes"

#include <effect.h>


GlobalData(BiModalFilter, "2 modal filters with nonlinear feedback");

class BiModalFilter : public ape::Effect
{
  
public:

  BiModalFilter() {}

  // Shorthands for convenience to reduce boilerplate:
  using Par = ape::Param<float>;
  using Rng = ape::Range;
  using Map = Rng::Mapping;

  // Why are the parameters public? is there a compelling reason for this? perhaps the gui widgets 
  // need to access them? -> figure out
  Par parGain{     "Gain",     Rng(-48,  12)              }; // in dB
  Par parFeedback{ "Feedback", Rng(-1,   +1)              }; // as raw factor
  Par parFreq1{    "Freq1",    Rng(20,   20000, Map::Exp) }; // in Hz
  Par parFreq2{    "Freq2",    Rng(20,   20000, Map::Exp) };
  Par parAmp1{     "Amp1",     Rng(-1,   +1)              }; // as raw factor
  Par parAmp2{     "Amp2",     Rng(-1,   +1)              };
  Par parDecay1{   "Decay1",   Rng(0.1,  1000,  Map::Exp) }; // in milliseconds
  Par parDecay2{   "Decay2",   Rng(0.1,  1000,  Map::Exp) };
  Par parPhase1{   "Phase1",   Rng(-180, +180)            }; // in degrees
  Par parPhase2{   "Phase2",   Rng(-180, +180)            };
  // How can we set the parameters to reasonable default values?


private:

  RAPT::rsModalFilter<float, float> mf1, mf2;
  float sampleRate;
  float yOld;        // output signal from previous sample
  float ya, yd;      // average and difference between current and previous output
  
  /** Resets the internal state. */
  void reset()
  {
    mf1.reset();
    mf2.reset();
    yOld = ya = yd = 0.f;
  }

  /** Feedback nonlinearity that produces upward spikes at all (upward and downward) 
  zero-crossings. How spikey the are is controlled by the local parameter t in the code. 
  ToDo:
  -make that t-parameter a user parameter "feedback shape" or "curve" or something
  -maybe introduce other feedback functions (see Experiments.cpp in the Research codebase), let the 
   user select them via an int parameter  */
  inline float spikesAtAllZeros(float y, float y_t) // y_t not used, but may be in other functions
  { 
    float t  = 0.99f;                  // shape parameter - values close to 1 make the feedback spikier
    float v  = 1.f-y*y;
    float tv = t*v;
    float z  = (tv-v)/(2.f*tv-t-1.f);  // rationally mapped v
    return z;
  };

  /** Produces one output sample at a time, taking an input sample x as input and a feedback 
  parameter fb. */
  inline float getSample(float x, float fb)
  {
    // Compute output:
    float g  = fb * spikesAtAllZeros(ya, yd);  // nonlinear feedback signal
    float y1 = mf1.getSample(x + g);           // output of filter 1
    float y2 = mf2.getSample(x + g);           // output of filter 2
    float y  = y1 + y2;                        // overall output

    // Update state and return output:
    ya   = 0.5f*(y + yOld);                    // average
    yd   = 0.5f*(y - yOld);                    // approximate derivative
    yOld = y;
    return y;
  }

  //-----------------------------------------------------------------------------------------------
  // \name Overriden callbacks (why are they private?) 

  void start(const ape::IOConfig& cfg) override
  {
    sampleRate = (float)cfg.sampleRate;
    reset();
  }

  void process(ape::umatrix<const float> ins, ape::umatrix<float> outs, size_t numFrames) override
  {
    // Pull the values of the parameters for this block and update the DSP objects:
    const float gain = RAPT::rsDbToAmp((float)parGain);
    const float fb   = (float)parFeedback;
    mf1.setModalParameters(parFreq1, parAmp1, 0.001f*parDecay1, parPhase1, sampleRate);
    mf2.setModalParameters(parFreq2, parAmp2, 0.001f*parDecay2, parPhase2, sampleRate);

    // Loop over the sample frames:
    const auto numChannels = sharedChannels();     // shared? with whom?
    for(size_t n = 0; n < numFrames; ++n)
    {
      float x = ins[0][n];                    // read input - we just use the 1st channel
      float y = gain * getSample(x, fb);      // compute output
      for(size_t c = 0; c < numChannels; ++c) // write generated sample into all channels
        outs[c][n] = y;                 
    }

    clear(outs, numChannels); // what does this do? clear unused channels?
  }
};

// ToDo: factor out a class for the pure DSP process