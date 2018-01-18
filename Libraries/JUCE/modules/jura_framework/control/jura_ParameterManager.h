#ifndef jura_ParameterManager_h
#define jura_ParameterManager_h

/** This class serves as baseclass fo objects that maintain an array of Parameters

\todo: maybe do state-management here (save/load settings)
\todo; maybe maintain an array of sub-ParameterManagers */

class JUCE_API ParameterManager : public ParameterObserver
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. You need to pass a CriticalSection that will be used for accesses to the 
  parameters. */
  ParameterManager(CriticalSection *lockToUse);

  /** Destructor. */
  virtual ~ParameterManager();

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

  /** Returns the pointer to the mutex lock that is used by the parameters. */
  CriticalSection* getCriticalSection() const { return lock; }

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
  virtual void parameterWillBeDeleted(Parameter* parameterThatWillBeDeleted);


protected:

  /** An array of our observed parameters. */
  std::vector<Parameter*> parameters;

  /** Returns the index of the parameter in the array or -1 if the parameter was not found .*/
  int getParameterIndex(Parameter *parameterToLookFor);

  /** A mutex-lock for accesses to the vector of observed parameters. */
  CriticalSection* lock = nullptr;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ParameterManager)
};

#endif 
