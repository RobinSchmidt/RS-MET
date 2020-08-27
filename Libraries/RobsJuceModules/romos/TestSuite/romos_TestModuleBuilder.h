#ifndef romos_TestModuleBuilder_h
#define romos_TestModuleBuilder_h

//#include "../framework/romos_ContainerModule.h"

namespace rsTestRomos
{

/** This class contains a bunch of static functions that create and return a pointer to modules 
that are programmatically created. They are intended to be used for automated unit test and 
demonstrate the usage of the framework. The caller must take care to eventually delete the returned 
modules again. */

class TestModuleBuilder
{

public:

  // todo: pass the string as 1st argument, update the CodeGenerator accordingly:

  static romos::Module* createWrappedAdder(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createGain(const rosic::rsString &name, int x, int y, bool polyphonic);  // simple gain
  static romos::Module* createSumDiff(const rosic::rsString &name, int x, int y, bool polyphonic);  // sum and difference
  static romos::Module* createWrappedSumDiff(const rosic::rsString &name, int x, int y, bool polyphonic);  // sum and difference, wrapped into container
  static romos::Module* createSumDiffProd(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createWrappedSumDiffProd(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createWrappedAdderN(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createDifferences(const rosic::rsString &name, int x, int y, bool polyphonic);  // some differences of inputs
  static romos::Module* createSummedDiffs(const rosic::rsString &name, int x, int y, bool polyphonic);  // some sums of differences
  static romos::Module* createMovingAverage(const rosic::rsString &name, int x, int y, bool polyphonic);  // simple, two-point moving average filter
  static romos::Module* createDelayedConnection(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createLeakyIntegrator(const rosic::rsString &name, int x, int y, bool polyphonic);  // leaky integrator filter
  static romos::Module* createTestFilter1(const rosic::rsString &name, int x, int y, bool polyphonic);  // some arithmetic combinations of MA and LI
  static romos::Module* createBiquadMacro(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createPinSortingInner(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createPinSortTest(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createDifferencer(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createImpulse(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createBlip(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createMonoToPoly(const rosic::rsString &name, int x, int y, bool polyphonic);  // monophonic constant into polyphonic minus
  static romos::Module* createVoiceCombiner(const rosic::rsString &name, int x, int y, bool polyphonic);  // polyphonic constant into voice combiner
  static romos::Module* createIn3Out5(const rosic::rsString &name, int x, int y, bool polyphonic);  // 3 inputs, 5 outputs, computes some diffs
  static romos::Module* createGateAndKill(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createTriggerAndKill(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createGatedNoteFrequency(const rosic::rsString &name, int x, int y, bool polyphonic);


  // test modules for system tests:
  static romos::Module* createVoiceKill(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createIn1Out2(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createPolyBlipStereo(const rosic::rsString &name, int x, int y, bool polyphonic);

  static romos::Module* createNoteFilter(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createGatedNoise(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createNoiseFlute(const rosic::rsString &name, int x, int y, bool polyphonic);



  // test modules for container manipulations:

  static romos::Module* createContainerize01(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createContainerize02(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createOneTwoThree(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createOutputDeletion(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createAddedConstants(const rosic::rsString &name, int x, int y, bool polyphonic);  // adds constants to yield the result of 900



  // test modules for performance tests:
  static romos::Module* createIdentityChain(const rosic::rsString &name, int x, int y, bool polyphonic);  // 20 chained identity modules in a container
  static romos::Module* createIdentityChainWithFeedback(const rosic::rsString &name, int x, int y, bool polyphonic);
  static romos::Module* createAdderChain(const rosic::rsString &name, int x, int y, bool polyphonic);  // 20 chained adder modules in a container
  static romos::Module* createAdderChainWithFeedback(const rosic::rsString &name, int x, int y, bool polyphonic);


  // createNestedContainers20 - 
  // container 1 in, 1 identity, 1 out - nested a 20-fold, outermost conatiner feeds the whole thing with a constant







  /** Creates a biquad filter that computes: y[n] = b0[n]*x[n] + b1[n]*x[n-1] + b2[n]*x[n-2] - a1[n]*y[n-1] - a2[n]*y[n-2]. */
  //romos::Module* createBiquadModule(ModuleProperties *propertiesToUse);

  /** Creates a Moog style ladder filter using 4 leaky integrators. */
  /*
  static romos::Module* createMoogFilterModule(ModuleProperties *propertiesToUse);

  static romos::Module* createContainerizeTestModule(ModuleProperties *propertiesToUse);

  // for testing uncontainerize:
  static romos::Module* createUnContainerizeInnerModule(ModuleProperties *propertiesToUse);
  static romos::Module* createUnContainerizeTestModule(ModuleProperties *propertiesToUse);

  // for testing I/O minimization:
  static romos::Module* createMinimizeInputsInnerModule1(ModuleProperties *propertiesToUse);
  static romos::Module* createMinimizeInputsTestModule1(ModuleProperties *propertiesToUse);
  */

};

}

#endif 
