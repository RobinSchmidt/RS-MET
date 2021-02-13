#include "../../../../RS-MET/Libraries/RobsJuceModules/rapt/rapt_for_ape.cpp"
// relative path from "Audio Programming Environment/includes"

#include <effect.h>


GlobalData(BiModalFilter, "2 modal filters with nonlinear feedback");

class BiModalFilter : public ape::Effect
{
public:


  BiModalFilter() {}
  
  using Par = ape::Param<float>;

  
  // why are the parameter public? is there a compelling reason fro this? -> figure out
  Par parGain{ "Gain", ape::Range(-48, 12) }; // in dB
  
  Par parFeedback{ "Feedback", ape::Range(-1, +1) };  // as raw factor

  Par parFreq1{ "Freq1", ape::Range(20, 20000, ape::Range::Mapping::Exp) }; // in Hz
  Par parFreq2{ "Freq2", ape::Range(20, 20000, ape::Range::Mapping::Exp) };
  
  Par parAmp1{ "Amp1", ape::Range(-1, +1) };  // as raw factor
  Par parAmp2{ "Amp2", ape::Range(-1, +1) };
  
  Par parDecay1{ "Decay1", ape::Range(0.1, 1000, ape::Range::Mapping::Exp) }; // in ms
  Par parDecay2{ "Decay2", ape::Range(0.1, 1000, ape::Range::Mapping::Exp) };
  
  Par parPhase1{ "Phase1", ape::Range(-180, +180) };  // in degrees
  Par parPhase2{ "Phase2", ape::Range(-180, +180) };


private:   

  RAPT::rsModalFilter<float, float> mf1, mf2;
  float sampleRate;
  float yOld;           // output signal for previous sample
  
  void reset()
  {
    mf1.reset();
    mf2.reset();
    yOld = 0;
  }
  
  
  // why is start and process private? 

  void start(const ape::IOConfig& cfg) override
  {
    sampleRate = (float)cfg.sampleRate;
    reset();
  }

  void process(ape::umatrix<const float> inputs, ape::umatrix<float> outputs, size_t frames) override
  {
    // Pull the values of the parameters for this block and update the DSP objects:
    const float gain = RAPT::rsDbToAmp((float)parGain);
    const float fb   = (float)parFeedback;
    mf1.setModalParameters(parFreq1, parAmp1, parDecay1, parPhase1, sampleRate);
    mf2.setModalParameters(parFreq2, parAmp2, parDecay2, parPhase2, sampleRate);

    // Loop over the sample frames:
    const auto numChannels = sharedChannels();    
    for(std::size_t n = 0; n < frames; ++n)
    {
      float x  = inputs[0][n];
      float y1 = mf1.getSample(x);
      float y2 = mf2.getSample(x);
      
      // to do:
      // -include nonlinear feedback into the inputs to the filter
      
      
      float y = gain * (y1 + y2);

      // Write generated sample into all channels:
      for(std::size_t c = 0; c < numChannels; ++c)
        outputs[c][n] = y;
    }

    clear(outputs, numChannels); // what does this do? clear unused channels?
  }
};