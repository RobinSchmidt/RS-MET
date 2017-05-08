#ifndef romos_ConcretePerformanceTests_h
#define romos_ConcretePerformanceTests_h

#include "romos_PerformanceTest.h"


namespace romos
{

  /**

  This file contains concrete subclasses of PerformanceTest. 
  maybe split into AtomicPerformanceTests and FrameworkPerformanceTests

  */


  /** A biquad made from basic modules (unit-delay, multiply, etc.).  */
  class BiquadMacroPerformanceTest : public PerformanceTest
  {
  public:
    BiquadMacroPerformanceTest();
  };

  /** A chain made from 20 interconnected identity modules. */
  class IdentityChainPerformanceTest : public PerformanceTest
  {
  public:
    IdentityChainPerformanceTest();
  };

  /** Same as IdentityChainPerformanceTest but with an additional (feedback) connection between the last and the first identity module.
  This additional connection should cause the block-processing functions to fall back to per-frame functions. */
  class IdentityChainWithFeedbackPerformanceTest : public PerformanceTest
  {
  public:
    IdentityChainWithFeedbackPerformanceTest();
  };

  /** Similar to IdentityChainPerformanceTest but with adder modules. */
  class AdderChainPerformanceTest : public PerformanceTest
  {
  public:
    AdderChainPerformanceTest();
  };

  /** Similar to AdderChainWithFeedbackPerformanceTest but with adder modules. */
  class AdderChainWithFeedbackPerformanceTest : public PerformanceTest
  {
  public:
    AdderChainWithFeedbackPerformanceTest();
  };




  // maybe use a macro (like CREATE_PERFORMANCE_TEST(Biquad) ) to create the boilerplate code here, too - maybe use a namespace to make
  // the "::" available as separator for the ClassName :

  /** Tests the performance of the FirstOrderFilter module. */
  class FirstOrderFilterPerformanceTest : public PerformanceTest
  {
  public:
    FirstOrderFilterPerformanceTest();
  };

  /** Tests the performance of the Biquad module. */
  class BiquadPerformanceTest : public PerformanceTest
  {
  public:
    BiquadPerformanceTest();
  };

  class BiquadDesignerPerformanceTest : public PerformanceTest
  {
  public:
    BiquadDesignerPerformanceTest();
  };

  class BandlimitedImpulseTrainPerformanceTest : public PerformanceTest
  {
  public:
    BandlimitedImpulseTrainPerformanceTest();
  };

  class SawOscillatorPerformanceTest : public PerformanceTest
  {
  public:
    SawOscillatorPerformanceTest();
  };


} // end namespace romos

#endif 
