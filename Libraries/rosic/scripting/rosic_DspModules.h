#ifndef rosic_DspModules_h
#define rosic_DspModules_h

/**

This file defines custom functions for the ExpressionEvaluator, which are used for digital 
signal processing tasks - these 'functions' often are actually DSP block with their own memory as
opposed to the instantaneous functions defined in ExpressionEvaluatorFunctions.h.

*/

// rosic-indcludes:
#include "rosic_ExpressionEvaluatorFunctions.h"
#include "../delaylines/rosic_IntegerDelayLine.h"
#include "../generators/rosic_SineOscillatorStereo.h"
#include "../filters/rosic_OnePoleFilterStereo.h"

//using namespace rosic;


namespace rosic
{




  //-----------------------------------------------------------------------------------------------
  // lowpass6:

  class FunctionNodeLowpass6 : public FunctionNode
  {
  public:
    FunctionNodeLowpass6(Expression *expression) : FunctionNode(expression)
    { 
      SetArgumentCount(2, 2, 0, 0, 0, 0); 
      filter.setMode(OnePoleFilterStereo::LOWPASS);
    }

    double DoEvaluate()
    { 
      // interpret inputs:
      setFrequency( m_nodes[1]->Evaluate() );

      // generate output:
      return filter.getSample( m_nodes[0]->Evaluate() );
    }

  protected:

    void setFrequency(double newFrequency)
    {
      if( newFrequency != filter.cutoff )
        filter.setCutoff(newFrequency);
    }

    OnePoleFilterStereo filter;
  };

  class FunctionFactoryLowpass6 : public FunctionFactory
  {
  public:
    std::string GetName() const { return "lowpass6"; }
    FunctionNode* DoCreate(Expression *expression) 
    { return new FunctionNodeLowpass6(expression);  }
  };



  //-----------------------------------------------------------------------------------------------
  // sineOscillator:

  class FunctionNodeSineOscillator : public FunctionNode
  {
  public:
    FunctionNodeSineOscillator(Expression *expression) : FunctionNode(expression)
    { 
      SetArgumentCount(1, 2, 0, 0, 0, 0);  
      osc.trigger();
    }

    double DoEvaluate()
    { 
      // interpret inputs:
      if( m_nodes.size() < 2 )  
        setStartPhase(0.0);
      else
        setStartPhase( m_nodes[1]->Evaluate() );
      setFrequency( m_nodes[0]->Evaluate() );

      // generate output:
      double outL, outR;
      osc.getSampleFrameStereo(&outL, &outR); 
      return outR;
    }

  protected:

    void setFrequency(double newFrequency)
    {
      if( newFrequency != osc.frequency )
        osc.setFrequency(newFrequency); 
    }
    void setStartPhase(double newStartPhase)
    {   
      osc.startPhase = newStartPhase; // just an assignment -> cheap -> no conditional
    }

    SineOscillatorStereo osc;
  };

  class FunctionFactorySineOscillator : public FunctionFactory
  {
  public:
    std::string GetName() const                    { return "sineOscillator"; }
    FunctionNode* DoCreate(Expression *expression) 
    { return new FunctionNodeSineOscillator(expression);  }
  };





  //-----------------------------------------------------------------------------------------------
  // sampleDelay:

  class FunctionNodeSampleDelay : public FunctionNode
  {
  public:
    FunctionNodeSampleDelay(Expression *expression) : FunctionNode(expression)
    { SetArgumentCount(1, 2, 0, 0, 0, 0);  }

    double DoEvaluate()
    { 
      if( m_nodes.size() < 2 )  
        delayLine.setDelayInSamples(1);
      else
        delayLine.setDelayInSamples((int) m_nodes[1]->Evaluate());

      double out = delayLine.getSample( m_nodes[0]->Evaluate() );
      return out; 
    }

  protected:

    IntegerDelayLine delayLine;

  };

  class FunctionFactorySampleDelay : public FunctionFactory
  {
  public:
    std::string GetName() const                    { return "sampleDelay"; }
    FunctionNode* DoCreate(Expression *expression) 
    { return new FunctionNodeSampleDelay(expression);  }
  };

} // end of namespace rosic



#endif  