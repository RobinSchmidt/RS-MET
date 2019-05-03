#ifndef jura_TuningFileManager_h
#define jura_TuningFileManager_h

/** This class manages tuning files in the .tun format and in an proprietary .xml format.

\todo:
-introduce a factor and a divisor for all frequencies
-introduce extrpolation by octave-equivalence (or some other equivalent interval) with a
 tag ExtrapolationFactor
-allow expressions to define tuneing frequencies (introduce a child-element with a list of
 constants)

...420 is good as base-freq (divisible by 2,3,5,7 and also by 4,6,10,14)
->Factor=440, Divisor=420 -> 440/420 = 1.0476... as scaler to bring it back to 440 Hz basefreq */

class TuningFileManager : public FileManager
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  TuningFileManager();

  /** Destructor. */
  virtual ~TuningFileManager();

  //---------------------------------------------------------------------------------------------
  // others:

  /** Assigns the rosic::TuningTable object which will be affected when a new tuning file is
  loaded. */
  virtual void assignTuningTable(rosic::TuningTable* newTable);

  /** Triggers the loading of a new tuning table from an xml file. */
  virtual bool loadTuningFromFile(const juce::File& fileToLoadFrom);

  //---------------------------------------------------------------------------------------------
  // FileManager overrides:

  /** Loads a new tuning-file and informs whether or not this was operation was successful. */
  virtual bool loadFile(const juce::File& fileToLoad);

  /** Saves the current tuning-data into a file. !!!NOT YET FUNCTIONAL!!! */
  virtual bool saveToFile(const juce::File& fileToSaveTo);

protected:

  /** Pointer to a rosic::TuningTable object which will be affected when a new tuning was
  loaded. */
  rosic::TuningTable* theTable;

  juce_UseDebuggingNewOperator;
};

#endif  