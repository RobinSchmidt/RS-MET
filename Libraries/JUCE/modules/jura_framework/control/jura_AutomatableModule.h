#ifndef jura_AutomatableModule_h
#define jura_AutomatableModule_h



/** This class serves as base class for effect- or instrument plugins that provide support for  
parameter automation. It derives from ParameterObserver and extends this class by a vector of 
pointers to Parameter objects - this vector is the place where we keep track of the parameters we 
are interested in. The class is still abstract in that it does not provide an implementation for 
the inherited (purely virtual) parameterChanged() method.

\todo factor out a baseclass ParameterManager - can be used for all classes that need to maintain 
a set of parameters, these classes need not to be audio related - this should perhaps also be a 
subclass of StateManager - OR rename this class ParametrizedObject/Module and get rid of the MIDI 
stuff -> move it to AudioModule


*/

class JUCE_API AutomatableModule : public ParameterObserver
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. You need to pass a CriticalSection that will be used for accesses to the 
  parameters. */
  AutomatableModule(CriticalSection *lockToUse);

  /** Destructor. */
  virtual ~AutomatableModule();

  //-----------------------------------------------------------------------------------------------
  // MIDI controller stuff:

  /** Assigns a MIDI controller to one of the observed parameters. */
  virtual void assignMidiController(const juce::String& nameOfParameter, int controllerNumber);
  // move to AudioModule

  /** Receives MIDI controller messages and dispatches them to the appropriate to one of the
  Parameter object. */
  virtual void setMidiController(int controllerNumber, float controllerValue);
  // move to AudioModule

  /** Reverts all observed parameters to their default settings. */
  virtual void revertToDefaultMapping();
  // move to AudioModule

  //-----------------------------------------------------------------------------------------------
  // retrieve pointers to the observed parameters:

  /** Retrieves a pointer to an Parameter object which has a given name - if no
  parameter with the given name exists in the vector of observed parameters, NULL will be
  returned. */
  virtual Parameter* getParameterByName(const juce::String& nameOfParameter) const;

  /** Retrieves a pointer to an Parameter object which has a given index - if no
  parameter with the given index exists in the vector of observed parameters, NULL will be
  returned. */
  virtual Parameter* getParameterByIndex(int indexOfParameter) const;

  /** Retrieves the index of a parameter in the array of observed parameters, returns -1 if the
  parameter was not found. */
  virtual int getIndexOfParameter(Parameter* parameterToRetrieveIndexOf) const;

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the number of automatable parameters in the vector. */
  virtual int getNumParameters() const;

  //-----------------------------------------------------------------------------------------------
  // add/remove observed parameters:

  /** Adds a pointer to an Parameter object to the array of observed parameters and registers this 
  instance as listener to the passed parameter and sets the Parameter up to use the same mutex 
  lock as this object. */
  virtual void addObservedParameter(Parameter *parameterToAdd);
  // rename to addParameter

  /** Removes a pointer to an Parameter object from the array of observed parameters and optionally 
  deletes the object itself. */
  virtual void removeObservedParameter(Parameter *parameterToRemove, bool deleteObject);
  // rename to removeParameter

  /** Removes all the pointers to the observed parameters and optionally deletes the objects 
  themselves. */
  virtual void removeAllObservedParameters(bool deleteObjects);
  // rename to removeAllParameters

  /** Overrides the inherited parameterChanged function in order to pass the new value of the 
  parameter to the appropriate handler
  function. */
  virtual void parameterChanged(Parameter *parameterThatHasChanged);

  /** Overrides the inherited parameterIsGoingToBeDeleted function in order to delete the 
  to-be-deleted parameter from our array of observed parameters (if it is in there, otherwise does 
  nothing). */
  virtual void parameterIsGoingToBeDeleted(Parameter* parameterThatWillBeDeleted);


protected:

  /** An array of our observed parameters. */
  juce::Array<Parameter*, CriticalSection> observedParameters;

  /** Returns the index of the parameter in the array or -1 if the parameter was not found .*/
  int getParameterIndex(Parameter *parameterToLookFor);

  /** This is for retrieving the number of inherited parameters - each sub(sub...etc.) class
  should write its number of parameters (including all the inherited ones) into a new slot of
  this array (at the end of initializeAutomatableParameters) - this facilitates later extensions
  of a baseclass's parameter set. */
  //juce::Array<int> numInheritedParameters;

  /** Assuming that subclasses correctly append their number of parameters to this array after
  adding all their parameters
  virtual int getNumInheritedParameters()
  { return numInheritedParameters[numInheritedParameters.size()-1]; }

  /** A mutex-lock for accesses to the vector of observed parameters. */
  CriticalSection *lock = nullptr;     // mutex to access the wrapped core dsp object

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(AutomatableModule)
};

#endif 
