#ifndef jura_StateFileManager_h
#define jura_StateFileManager_h

/** This class manages files that presumably contain application settings, plugin presets, etc.
(that is: states or mementos) stored as xml-files. */

class JUCE_API StateFileManager : virtual public FileManager, public StateManager
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  StateFileManager();

  /** Destructor. */
  virtual ~StateFileManager();

  //-----------------------------------------------------------------------------------------------
  // FileManager overrides:

  /** Loads a new state-file and informs whether or not this was operation was successful. */
  virtual bool loadFile(const juce::File& fileToLoad);

  /** Saves the current state into a file and informs whether or not this was operation was
  successful.  */
  virtual bool saveToFile(const juce::File& fileToSaveTo);



  /** Tries to load the xml-file passed in the argument and if this was successful, it calls
  setStateFromXml with the XmlElement obtained from the file. */
  virtual bool loadStateFromXmlFile(const juce::File& fileToLoadFrom);

  /** Retrieves the current state by calling getStateAsXml() and saves this state into an
  xml-file as specified by the argument. */
  virtual bool saveStateToXmlFile(const juce::File& fileToSaveTo);

protected:
  juce_UseDebuggingNewOperator;
};

#endif  