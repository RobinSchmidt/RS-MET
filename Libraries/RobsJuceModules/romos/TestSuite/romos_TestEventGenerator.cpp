#include "romos_TestEventGenerator.h"
//using namespace rsTestRomos;

namespace rsTestRomos
{

std::vector<NoteEvent> TestEventGenerator::generateNoteOnOffPair(unsigned int key, unsigned int velocity,
  unsigned int deltaFramesForNoteOn, unsigned int durationInFrames)
{
  std::vector<NoteEvent> result;
  if(durationInFrames == 0)
    return result; // return an empty vector - note-ons with simultaneous note-offs are discarded
  result.push_back(NoteEvent(deltaFramesForNoteOn, key, velocity));
  result.push_back(NoteEvent(deltaFramesForNoteOn + durationInFrames, key, 0));
  return result;
}

std::vector<NoteEvent> TestEventGenerator::generateSimultaneousNotes(unsigned int key, unsigned int velocity,
  unsigned int deltaFramesForNoteOn, unsigned int durationInFrames,
  unsigned int numNotes, unsigned int noteSpacing)
{
  std::vector<NoteEvent> result;
  unsigned int currentKey = key;
  for(unsigned int i = 1; i <= numNotes; i++)
  {
    result = mergeEvents(result, generateNoteOnOffPair(currentKey, velocity, deltaFramesForNoteOn, durationInFrames));
    currentKey += noteSpacing;
  }
  return result;
}

std::vector<NoteEvent> TestEventGenerator::mergeEvents(const std::vector<NoteEvent>& eventArray1,
  const std::vector<NoteEvent>& eventArray2)
{
  std::vector<NoteEvent> result;
  result.reserve(eventArray1.size() + eventArray2.size());
  unsigned int i;
  for(i = 0; i < eventArray1.size(); i++)
    result.push_back(eventArray1[i]);
  for(i = 0; i < eventArray2.size(); i++)
    result.push_back(eventArray2[i]);
  std::sort(result.begin(), result.end(), noteEventLessByDeltaFrames);
  return result;
}

void TestEventGenerator::convertNoteEventsToStartsAndDurations(const std::vector<NoteEvent>& events,
  std::vector<NoteEvent>& noteOns, std::vector<int>& durations)
{
  unsigned int i;
  noteOns.clear();
  for(i = 0; i < events.size(); i++)
  {
    if(events[i].getVelocity() != 0)
      noteOns.push_back(events[i]);
  }

  durations.reserve(noteOns.size());
  for(i = 0; i < noteOns.size(); i++)
  {
    int noteOffIndex = findIndexOfMatchingNoteOff(events, noteOns[i]);
    if(noteOffIndex > -1)
      durations.push_back(events[noteOffIndex].getDeltaFrames() - noteOns[i].getDeltaFrames());
    else
      durations.push_back(INT_MAX);  // no matching note-off - potentially infinite duration, but INT_MAX is the largest we have
  }
}

int TestEventGenerator::findIndexOfMatchingNoteOff(const std::vector<NoteEvent>& events, NoteEvent noteOnEvent)
{
  unsigned int i;
  unsigned int startIndex = rosic::findElement(events, noteOnEvent) + 1;
  for(i = startIndex; i < events.size(); i++)
  {
    if(events[i].isNoteOff() && events[i].getKey() == noteOnEvent.getKey())
      return i;
  }
  return -1;
}

}