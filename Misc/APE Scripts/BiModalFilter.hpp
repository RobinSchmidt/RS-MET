#include "../../../../RS-MET/Libraries/RobsJuceModules/rapt/rapt_for_ape.cpp"
// relative path from "Audio Programming Environment/includes"

#include <effect.h>


GlobalData(BiModalFilter, "2 modal filters with nonlinear feedback");

class BiModalFilter : public ape::Effect
{
public:


  BiModalFilter() {}

  ape::Param<float> param{ "Gain", ape::Range(0, 4) };  // why is the parameter public?





private:   

  RAPT::rsModalFilter<float, float> mf1, mf2;
  float sampleRate;
  
  
  // why is start and process private? 

  void start(const ape::IOConfig& cfg) override
  {
    sampleRate = (float)cfg.sampleRate;
  }

  void process(ape::umatrix<const float> inputs, ape::umatrix<float> outputs, size_t frames) override
  {
    const auto numChannels = sharedChannels();
    const float gain = param;


    for(std::size_t n = 0; n < frames; ++n)
    {
      // to do:
      // -set up the filter's accorind to the user parameters
      // -update the filter coefficients
    
    
      float y1 = mf1.getSample(inputs[0][n]);
      
      // to do:
      // -include nonlinear feedback into the inputs to the filter
      
      float y  = y1;

      for(std::size_t c = 0; c < numChannels; ++c)
        outputs[c][n] = gain * y;
    }

    clear(outputs, numChannels); // what does this do? clear unused channels?
  }
};