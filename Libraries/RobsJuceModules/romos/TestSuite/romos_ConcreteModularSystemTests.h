#ifndef romos_ModularSystemTests_h
#define romos_ModularSystemTests_h

//#include "romos_ModularSystemTest.h"
//#include "romos_GenerateDesiredOutput.h"

namespace rsTestRomos
{

  /**

  This file contains concrete subclasses of ModularSystemTest. 

  */



  /** A patch that just passes the inputs through to the outputs. */
  class BypassTest : public ModularSystemTest
  {
  public:
    BypassTest(const char *testName = "BypassTest");
  protected: 
    virtual void createAndConnectTestChildModules();
    virtual void fillDesiredOutputSignalArrays();
  };


  /** Like BypassTest, but adds a child-module to the TopLevelModule after establishing the connections between inputs and outputs. There
  was a bug that turned the output silent when adding a child module. It was due to attempting to re-connect outside connections of a 
  TopLevelModule in ContainerModule::sortChildModuleArray - the method has now been overriden in TopLevelModule to fix it. */
  class BypassWithChildTest : public BypassTest
  {
  public:
    BypassWithChildTest();
  protected: 
    virtual void createAndConnectTestChildModules();
  };


  /** A polyphonic blip (but only 1 note is played). This test also sorts the childmodule  array of the embedded In1Out2 object - this is 
  done to expose a bug which occurred in updateHasDelayedConnectionsFlag (and was fixed since then). */
  class PolyBlipStereoTest : public ModularSystemTest
  {
  public:
    PolyBlipStereoTest();
  protected: 
    virtual void createAndConnectTestChildModules();
    virtual void fillDesiredOutputSignalArrays();
    ContainerModule *polyBlipStereo;
  };


  /** A stereophonic flute-like instrument using 2 noise sources feeding 2 bandpasses. */
  class NoiseFluteTest : public ModularSystemTest
  {
  public:
    NoiseFluteTest();
  protected: 
    virtual void createAndConnectTestChildModules();
    virtual void fillDesiredOutputSignalArrays();
    ContainerModule *noiseFlute;
  };


} // end namespace romos

#endif 
