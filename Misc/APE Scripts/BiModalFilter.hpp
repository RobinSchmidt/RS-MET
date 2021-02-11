#include "../../../../RS-MET/Libraries/RobsJuceModules/rapt/rapt_for_ape.cpp"
// relative path from "Audio Programming Environment/includes"

#include <effect.h>


GlobalData(BiModalFilter, "2 modal filters with nonlinear feedback");

class BiModalFilter : public ape::Effect
{
public:


  BiModalFilter() {}

  ape::Param<float> param{ "Gain", ape::Range(0, 4) };





private:

  RAPT::rsModalFilter<float, float> mf1, mf2;
  float sampleRate;
  

  void start(const ape::IOConfig& cfg) override
  {
    sampleRate = (float)cfg.sampleRate;
  }

  void process(ape::umatrix<const float> inputs, ape::umatrix<float> outputs, size_t frames) override
  {
    const auto numChannels = sharedChannels();
    const float gain = param;

    for(std::size_t c = 0; c < numChannels; ++c)
    {
      for(std::size_t n = 0; n < frames; ++n)
        outputs[c][n] = gain * inputs[c][n];
    }

    clear(outputs, numChannels); // what does this do? clear unused channels?
  }
};