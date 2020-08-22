#ifndef romos_GlobalFrameworkTests_h
#define romos_GlobalFrameworkTests_h

//#include "romos_UnitTest.h"
//#include "romos_TestEventGenerator.h"

//#include "romos_ConcreteProcessingTests.h"
//#include "../framework/romos_ModuleFactory.h"
//#include "../romos.h"

namespace rsTestRomos
{

  /**

  This file contains test-classes that test some global stuff, such as the voice-allocator, etc.

  */


  class VoiceAllocatorTest : public UnitTest
  {
  public:
    VoiceAllocatorTest();
    virtual ~VoiceAllocatorTest();
  protected:
    virtual bool runTest();
    virtual bool testStealOldestWithoutRetrigger();
    virtual bool testStealOldestWithRetrigger();

    virtual bool areAllNoteOnTriggerFlagsUnchecked(const VoiceAllocator &voiceAllocator);
    virtual bool isNoteOnTriggerFlagCheckedExclusively( const VoiceAllocator &voiceAllocator, int voiceIndexThatShouldHaveFlagSet);
    virtual bool isNoteOffTriggerFlagCheckedExclusively(const VoiceAllocator &voiceAllocator, int voiceIndexThatShouldHaveFlagSet);
  };


  class TriggerAndKillTest : public ProcessingTest
  {
  public:
    TriggerAndKillTest();
    virtual ~TriggerAndKillTest();
  protected:
    virtual bool runTest();
    romos::Module         *triggerAndKillModule;
    //romos::TopLevelModule *top
    //romos::ModularSynth *theSynth;  
  };



  class TopLevelModuleTest : public UnitTest
  {
  public:
    TopLevelModuleTest();
    virtual ~TopLevelModuleTest();
  protected:
    virtual bool runTest();
    romos::TopLevelModule *moduleToTest;
  };








} // end namespace romos

#endif 
