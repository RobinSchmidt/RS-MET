#ifndef rosic_Module_h
#define rosic_Module_h

// rosic-indcludes:
#include "../datastructures/rosic_Array.h"
#include "rosic_MutexLock.h"

namespace rosic
{

  // this file contains the classes that establish the infrastructure for handling modules in a
  // (semi) modular framework, including modulation routing

  //===============================================================================================
  // class ModulationSource:

  /**

  This class serves as baseclass for various (routable) modulation sources such as envelopes, LFOs,
  etc. It is meant to be used in conjunction with classes ModulatableParameter and 
  ModulationRouter.

  */

  class ModulationSource
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    ModulationSource() {}   

    /** Destructor. */
    virtual ~ModulationSource() {} 

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Override this to set up the sample-rate. */
    virtual void setSampleRate(double newSampleRate) = 0;

    /** Override this to set up the tempo, if your effect needs this. */
    virtual void setTempoInBPM(double newTempo) 
    { 
      newTempo = 0.f;  // to avoid unreferenced parameter warning
    }

    //---------------------------------------------------------------------------------------------
    // inquiry:


    //---------------------------------------------------------------------------------------------
    // processing:

    /** Generates one sample of the modulation signal at a time. */
    virtual double getSample();

    //---------------------------------------------------------------------------------------------
    // others:

    /** Override this to re-trigger the modulation generator. */
    virtual void trigger() {}

    //=============================================================================================

  protected:

  };


  //===============================================================================================
  // class ModulatableParameter:

  /**

  This class serves as baseclass for (routable) parameters. It is meant to be used in conjunction 
  with classes ModulationSource and ModulationRouter.

  */

  class ModulatableParameter
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    ModulatableParameter(const char* name = "unknown parameter", double initialValue = 0.0);

    /** Destructor. */
    //~ModulatableParameter();

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets the nominal value of the parameter (without any modulators applied). */
    void setNominalValue(double newValue) { nominalValue = newValue; }

    /** Initializes the instantaneous value with the nominal value. It is supposed to be called 
    before subsequent calls to addModulationSignal. */
    void initInstantaneousValue() { instantaneousValue = nominalValue; }

    /** Adds a modulation signal to the (instantaneous) value of this parameter. This function is 
    supposed to be used to accumulate all the modulation signals after the instantaneous value has
    been initialized via to its nominal value initInstantaneousValue. */
    void addModulationSignal(double valueToAdd) { instantaneousValue += valueToAdd; }

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the nominal value of this parameter. */
    double getNominalValue() const { return nominalValue; }

    /** Returns the instantaneous value of this parameter. */
    double getInstantaneousValue() const { return instantaneousValue; }

    /** Returns the name of this parameter as zero-terminated c-string. */
    char* getName() const { return name; }

    //=============================================================================================

  protected:

    double nominalValue, instantaneousValue;
    char   *name;

  };


  //===============================================================================================
  // class ModulationConnection:

  /**

  This class establishes a connection with adjustable strength between an ModulationSource object 
  and a ModulatableParameter object.

  */

  class ModulationConnection
  {

    friend class ModulationRouter;
  
  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    ModulationConnection(ModulationSource *sourceToConnect = NULL, 
      ModulatableParameter *parameterToConnect = NULL, double connectionStrength = 0.0)
    {
      source    = sourceToConnect;
      parameter = parameterToConnect;
      strength  = connectionStrength;
    }

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets the modulation-source for this connection. */
    void setSource(ModulationSource *newSourceToConnect) { source = newSourceToConnect; }

    /** Sets the modulation-target parameter for this connection. */
    void setSource(ModulatableParameter *newParameterToConnect) 
    { parameter = newParameterToConnect; }

    /** Sets the strength for this modulation connection. */
    void setStrength(double newStrength) { strength = newStrength; }

  protected:

    ModulationSource     *source;
    ModulatableParameter *parameter;
    double                strength;

  };


  //===============================================================================================
  // class ModulationRouter:

  /**

  This class is to be used to connect an arbitrary number of ModulationSource objects to an 
  arbirary number of ModulatableParameter objects.

  */

  class ModulationRouter
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    ModulationRouter();

    /** Destructor. */
    ~ModulationRouter();

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Establishes a new connection between a ModualtionSource and a ModulatableParameter with the given 
    (initial) strength. */
    void establishNewConnection(ModulationSource *sourceToConnect, 
      ModulatableParameter *parameterToConnect, double connectionStrength = 0.0);

    /** Removes the modulation connection with the given index. */
    void removeConnection(int index);


    //---------------------------------------------------------------------------------------------
    // inquiry:



    //---------------------------------------------------------------------------------------------
    // others:

    /** Iterates through all the modulation connections and applies the modulation-signals to their 
    respective target parameters, optionally initializing the parameter with its nominal value 
    before (this is usually what you want, unless you are using more than one ModulationRouter that
    can address the parameter in question). When the function returns, all parameters will have 
    their instantantaneous values set up. */
    void applyModulations(bool initWithNominalValue = true);

    //=============================================================================================

  protected:

    MutexLock mutex;

    //std::vector<ModulationConnection> connections;
    Array<ModulationConnection> connections;

  };


  //===============================================================================================
  // class Module:


  /**

  This class defines the interface (I/O functions, etc.) for audio-processing modules that 
  take a stereo-pair of input signals and produce a stereo-pair of output signals. Each concrete 
  effect module should be subclassed from this Module baseclass and implement these I/O 
  functions in a suitable manner. The class facilitates the use of the modules in a 
  (semi) modular framework such as Quadrifex.

  */

  class Module
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    Module() {}   

    /** Destructor. */
    virtual ~Module() {} 

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Override this to set up the sample-rate. */
    virtual void setSampleRate(double newSampleRate) = 0;

    /** Override this to set up the tempo, if your effect needs this. */
    virtual void setTempoInBPM(double newTempo) 
    {
      newTempo = 0.f;
    }

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the number of modulatable parameters. Override this in your subclass if you want
    to allow some parameter to be modulated over a modulation matrix (or similar routing 
    facility). */
    virtual int getNumModulatableParameters() { return modulatableParameters.getNumElements(); }

    /** Returns the name of one of the modulatable parameters as a c-string. */
    virtual char* getModulatableParameterName(int index);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Override this to render one output sample-frame from and input sample-frame at a time. The
    function is supposed to be wrapped into calls to acquireLock/releaseLock to make it 
    thread-safe.*/
    virtual void processSampleFrame(double *inOutL, double *inOutR) = 0;

    //---------------------------------------------------------------------------------------------
    // others:

    /** Override this to reset the internal state (buffers, oscillator-phases and such) of the 
    effect. */
    virtual void reset() = 0;

    /** Override this to re-trigger oscillators. */
    virtual void trigger() {}

    //=============================================================================================

  protected:

    Array<ModulatableParameter> modulatableParameters;

  };

} // end namespace rosic

#endif // rosic_Module_h
