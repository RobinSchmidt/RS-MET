#ifndef rosic_TuningTable_h
#define rosic_TuningTable_h

//#include "../_third_party/MarkHenning/TuningMap.h"
//#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{
  /**

  This class implements a note-to-frequency mapping via a table which maps the 128 (0...127)
  MIDI notes to their corresponding frequency. The class can be used to realize different
  tunings/scales. By default, it realizes the equal tempered scale.

  */

  class TuningTable
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    TuningTable();
    /**< Constructor. */

    ~TuningTable();
    /**< Destructor. */

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    void setName(const char* newName);
    /**< Sets the name (i.e. the relative path from the standard tuning file directory) of the
    current tuning table - this function is only for conviniently handling loading/saving and
    preset management in a plugIn-context. The char-array is supposed to have 'newLength'+1
    elements where the last element should be the terminating zero.*/

    void resetToDefaults();
    /**< Sets this tuning table into its default state (equal temperament, A4 = 440.0 Hz,
    name = "EqualTemperament"). */

    void assignFrequency(int note, double newFrequency);
    /**< Assigns a frequency to a particular note-number. */

    void assignFrequencies(double* newFrequencies);
    /**< Assigns all frequencies at once - the passed array is assumed to have 128 entries
    corresponding to the 128 MIDI notes. */

    void setMasterTuneA4(double newTuneA4);
    /**< Sets the tuning frequency for the A4 reference tone (usually 440 Hz). */

    bool loadFromTunFile(char* path);
    /**< Loads a tuning table from a .tun file. */

    //---------------------------------------------------------------------------------------------
    // note to frequency conversion:

    char* getName();
    /**< Returns the name of the tuning (i.e. the relative path from the standard tuning file
    directory) as a zero-terminated string. */

    bool isInDefaultState();
    /**< Informs whether or not this object is in its default state. */

    double getFrequency(int note);
    /**< Converts a MIDI note number to a corresponding frequency by reading out the table at the
    position given by 'note' */

    double getFrequency(double note);
    /**< Converts a MIDI note number to a corresponding frequency by reading out the table at the
    position given by a possibly non-integer note-number by linearly interpolating  between the two
    adjacent pitches (not the frequencies). */

    double getMasterTuneA4();
    /**< Returns the tuning frequency for the A4 reference tone (usually 440 Hz). */

  protected:

    double table[128];  // ToDo: use a symbolic constant tableLength
    /**< The table which contains the mapping from the note numbers to frequencies. */



    double masterTuneA4;
    /**< The master tuning frequency in Hz. */

    double detuneFactor;
    /**< A global detuning factor given by the ratio between the master tune for A4 and 440 Hz. */

    //bool defaultState;
    /**< A flag to indicate whether or not this object is in default state. */

    char* name;
    /**< A name assigned to the current tuning - this is supposed to be a relative path from the
    standard tuning file directory to the xml-file which represents the tuning. */

  };

} // end namespace rosic

#endif //  rosic_TuningTable_h
