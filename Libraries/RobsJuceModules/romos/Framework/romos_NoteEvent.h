#ifndef romos_NoteEvent_h
#define romos_NoteEvent_h


  /**

  This class is used to represents a note-events within the RoMoS library. Note events contain a key a velocity and an offset in 
  sample-frames which are to be interpreted with respect to the signal-block that is currently processed.

  \todo generalize this to MusicalEvent (can also be a control-change, pitch-bend, etc.)
  -> needs an eventType field, key and velocity fileds should be renamed to data1, data2 and interpreted differently according to the 
  event-type

  */

  class NoteEvent
  {

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. Creates a note event with the given initial values for the data fields. */
    NoteEvent(unsigned int deltaFrames = 0, unsigned int key = 0, unsigned int velocity = 0);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns true when this is a note-off event, false otherswise. Note-off events are internally represented by velocity == 0 
    (following the MIDI convention). */
    inline bool isNoteOff() const
    {
      return velocity == 0;
    }

    /** Returns the number of sample-frames that this event is offset with respect to the start-frame of the current block. */
    inline unsigned int getDeltaFrames() const
    {
      return deltaFrames;
    }

    /** Returns the note key as a value between 0 and 127 according to the MIDI convention. */
    inline unsigned int getKey() const
    {
      return key;
    }

    /** Returns the note velocity as a value between 0 and 127 according to the MIDI convention. */
    inline unsigned int getVelocity() const
    {
      return velocity;
    }

    //-------------------------------------------------------------------------------------------------------------------------------------
    // operators:

    /** Two note events are considered to be equal, if all data fields match. */
    bool operator==(const NoteEvent& other) const  
    {
      if( deltaFrames == other.deltaFrames && key == other.key && velocity == other.velocity )
        return true;
      else
        return false;
    }

  protected:

    unsigned int deltaFrames;
    unsigned int key;
    unsigned int velocity; 

  };

  /** Defines the less-than relation between note-events by considering the time of occurence (i.e. the deltaFrames field) as primary 
  indicator. When the deltaFrames fileds are equal, that is, the events are simultaneous, they will be ordered by their velocity where 
  lower velocities come first to make sure that note-offs always precede note-ons (this is important when a note ends and is immediatetly 
  re-triggered). If velocities are also equal, the key will be considered, where lower keys precede higher ones. */
  bool noteEventLessByDeltaFrames(NoteEvent left, NoteEvent right); 


#endif
