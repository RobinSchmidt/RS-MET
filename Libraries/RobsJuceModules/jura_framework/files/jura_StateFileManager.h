#ifndef jura_StateFileManager_h
#define jura_StateFileManager_h

// ToDo: rename file to FileManagers or ConcreteFileManagers and move the code for the 
// AudioFileManager in here, too. We don't want so many small source files.

/** This class manages files that presumably contain application settings, plugin presets, etc.
(that is: states or mementos) stored as xml-files. 

\todo 
Why do we need virtual inheritance from FileManager? Try to get rid of virtual or document, 
why it's needed! */

class JUCE_API StateFileManager : virtual public FileManager, public StateManager
{


public:

  //-----------------------------------------------------------------------------------------------
  // \name Lifetime

  /** Constructor. */
  StateFileManager();

  /** Destructor. */
  virtual ~StateFileManager();

  //-----------------------------------------------------------------------------------------------
  // \name Overrides:

  /** Loads a new state-file and informs whether or not this was operation was successful. */
  bool loadFile(const juce::File& fileToLoad) override;

  /** Saves the current state into a file and informs whether or not this was operation was
  successful.  */
  bool saveToFile(const juce::File& fileToSaveTo) override;

  //-----------------------------------------------------------------------------------------------
  // \name Xml load/save:

  /** Tries to load the xml-file passed in the argument and if this was successful, it calls
  setStateFromXml with the XmlElement obtained from the file. */
  virtual bool loadStateFromXmlFile(const juce::File& fileToLoadFrom);

  /** Retrieves the current state by calling getStateAsXml() and saves this state into an
  xml-file as specified by the argument. */
  virtual bool saveStateToXmlFile(const juce::File& fileToSaveTo);


protected:

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(StateFileManager)
};









#endif