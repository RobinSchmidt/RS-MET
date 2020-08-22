#ifndef romos_ConcretePerformanceTests_h
#define romos_ConcretePerformanceTests_h

//#include "romos_PerformanceTest.h"

namespace rsTestRomos
{

/** This file contains concrete subclasses of PerformanceTest. 
maybe split into AtomicPerformanceTests and FrameworkPerformanceTests  */


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

class Formula11PerformanceTest : public PerformanceTest
{
public:
  Formula11PerformanceTest() : PerformanceTest("Formula_1_1")
  { 
    //moduleFactory.registerModuleType(new FormulaModule_1_1TypeInfo); 
    // may have to be re-activated if the unit-tests are not run before the performance tests
    // ...find a better solution..

    moduleToTest = moduleFactory.createModule("Formula_1_1"); 

  }
};

class FormulaN1PerformanceTest : public PerformanceTest
{
public:
  FormulaN1PerformanceTest() : PerformanceTest("Formula_N_1")
  { 
    //moduleFactory.registerModuleType(new FormulaModule_N_1TypeInfo);
    // may have to be re-activated if the unit-tests are not run before the performance tests
    // ...find a better solution..

    moduleToTest = moduleFactory.createModule("Formula_N_1"); 
  }
};

class FormulaNMPerformanceTest : public PerformanceTest
{
public:
  FormulaNMPerformanceTest() : PerformanceTest("Formula_N_M")
  { moduleToTest = moduleFactory.createModule("Formula"); }
};

// The performance penalty for using FormulaNM instead of Formula11 is about 10% (for a trivial
// identity formula y=x)...hmm...i guess that justifies removing the _1_1 and _N_1 versions from
// Liberty and include only the _N_M case



} // end namespace romos

#endif 
