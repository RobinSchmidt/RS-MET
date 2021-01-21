#ifndef rosic_ModulationRouting_h
#define rosic_ModulationRouting_h

namespace rosic
{

  //===============================================================================================
  // class ModulationSource:

  /**

  This class serves as baseclass for various (routable) modulation sources such as envelopes, LFOs,
  etc. It is meant to be used in conjunction with classes ModulatableParameter and
  ModulationRouter.

  -the design doesn't seem quite right yet  -this is under cosntruction - maybe move to unfinished
   also the class Module

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
    virtual void setTempoInBPM(double /*newTempo*/) {}

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

    /** Returns the name of this parameter */
    std::string getName() const { return name; }

    ///** Returns the name of this parameter as zero-terminated c-string. */
    //char* getName() const { return name; }

    //=============================================================================================

  protected:

    double nominalValue, instantaneousValue;
    //char   *name;
    std::string name;

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
    rsDynamicArray<ModulationConnection> connections;

  };


} // end namespace rosic

#endif // rosic_ModulationRotuing_h
