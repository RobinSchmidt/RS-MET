#ifndef romos_TestEventGenerator_h
#define romos_TestEventGenerator_h

//#include "../framework/romos_NoteEvent.h"

namespace rsTestRomos
{

  /** 

  This class generates certain sequences of events that are used to drive the tests.

  */

  class TestEventGenerator
  {

  public:

    /** Returns a vector containing a note-on and a corresponding note-off (indicated by velocity == 0) - if zero is passed as duration, it 
    will return an empty vector. */
    static std::vector<romos::NoteEvent> generateNoteOnOffPair(unsigned int key, unsigned int velocity, 
      unsigned int deltaFramesForNoteOn, unsigned int durationInFrames);

    /** Generates a bunch of simultaneous notes with equal spacing between the individual notes, given by noteSpacing. */
    static std::vector<romos::NoteEvent> generateSimultaneousNotes(unsigned int key, unsigned int velocity, 
      unsigned int deltaFramesForNoteOn, unsigned int durationInFrames,
      unsigned int numNotes, unsigned int noteSpacing);

    /** Merges two vectors of events. The result will be sorted by the time of occurence of the events. Simlutaneous events may occur in any
    order wihtin the array (\todo maybe have well defined ordering criterions in these cases, too). */
    static std::vector<romos::NoteEvent> mergeEvents(const std::vector<romos::NoteEvent> &eventArray1, const std::vector<romos::NoteEvent> &eventArray2);

    /** Converts an array of note-on events with corresponding note-off events into an array that contains only the note-ons and another array 
    that contains the corresponding durations (in sample-frames). */
    static void convertNoteEventsToStartsAndDurations(const std::vector<romos::NoteEvent> &events, std::vector<romos::NoteEvent> &noteOns, 
      std::vector<int> &durations);

    /** Given a vector of note-events and a particular note-on event, this function returns the index of the note-off event that corresponds to
    the given note-on event. If none is found, -1 will be returned. */
    static int findIndexOfMatchingNoteOff(const std::vector<romos::NoteEvent> &events, romos::NoteEvent noteOnEvent);

    /** Sorts the passed sequence of events by time of occurence. */
    static void sortEventsByTimeOfOccurence(const std::vector<romos::NoteEvent> &eventArray);

  };
  
} 

#endif 
