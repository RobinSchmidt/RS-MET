#include "UnitTestsToolChain.h"
using namespace juce;
using namespace jura;


void UnitTestToolChain::runTest()
{
  runTestVoiceManager();
}


void UnitTestToolChain::runTestVoiceManager()
{
  beginTest("rsVoiceManager");

  rsVoiceManager voiceMan;
  std::vector<double> voiceBuffer(2*voiceMan.getMaxNumVoices());
  voiceMan.setVoiceSignalBuffer(&voiceBuffer[0]);

  // Test voice killing:
  using KM = rsVoiceManager::KillMode;
  double killThresh = 0.2;
  double killTime   = 0.1;
  double sampleRate = 100;
  voiceMan.setKillMode(KM::afterSilence);
  voiceMan.setSampleRate(sampleRate);
  voiceMan.setKillThreshold(killThresh);
  voiceMan.setKillTime(killTime);
  int killSamples = voiceMan.getKillTimeSamples();

  using Msg = juce::MidiMessage;
  int    key1  = 10;
  float  vel1  = 1.0f;
  double vel1Q = voiceMan.quantize7BitUnsigned(vel1);
  int    active, releasing;

  // Trigger a voice, then simulate a generated signal for the voice above the kill-threshhold by
  // assigning that value to the first two slots of the voiceBuffer:
  voiceMan.handleMidiMessage(Msg::noteOn(1, key1, vel1));
  active    = voiceMan.getNumActiveVoices();
  releasing = voiceMan.getNumReleasingVoices();
  expectEquals(active,    1);
  expectEquals(releasing, 0);
  voiceBuffer[0] = voiceBuffer[1] = 1.1*killThresh;

  // Release the voice just triggered. It should nevertheless remain active (in release state),
  // because the signals in the corresponding voiceBuffer slots are above the kill threshold:
  voiceMan.handleMidiMessage(Msg::noteOn(1, key1, 0.f));
  active    = voiceMan.getNumActiveVoices();
  releasing = (int)voiceMan.getNumReleasingVoices();
  expectEquals(active,    1);
  expectEquals(releasing, 1);
  for(int n = 0; n < 2*killSamples; n++) {
    voiceMan.findAndKillFinishedVoices();               // should find no voice to kill...
    active    = voiceMan.getNumActiveVoices();          // ...so this should remain at 1
    releasing = voiceMan.getNumReleasingVoices();       // ...and this should be 1, too
    expectEquals(active,    1);
    expectEquals(releasing, 1); }

  // Now set the content of the voiceBuffer equal to the kill threshold. The voice should still 
  // remain active for "killSamples" samples and then turn off:
  voiceBuffer[0] = voiceBuffer[1] = killThresh;
  for(int n = 1; n < killSamples; n++) {
    voiceMan.findAndKillFinishedVoices();               // should find no voice to kill...
    active    = voiceMan.getNumActiveVoices();          // ...so this should remain at 1
    releasing = voiceMan.getNumReleasingVoices();       // ...and this should be 1, too
    expectEquals(active,    1);
    expectEquals(releasing, 1); }
  voiceMan.findAndKillFinishedVoices();               // should find and kill the voice
  active    = voiceMan.getNumActiveVoices();          // ...so this should go to 0
  releasing = voiceMan.getNumReleasingVoices();       // ...and this too
  expectEquals(active,    0);
  expectEquals(releasing, 0);

  // Maybe later test other kill modes. Currently the only other mode is the trivial 
  // "immediately" mode, but if we later come up with more modes, test for them should go here. But
  // i really can't think of any other meaningful modes at the moment.

  // Test voice stealing:
  using SM = rsVoiceManager::StealMode;
  voiceMan.setNumVoices(4);
  voiceMan.setStealMode(SM::oldest);
  voiceMan.setKillMode(KM::immediately);
  voiceMan.reset();
  int key2 = 20, key3 = 30, key4 = 40, key5 = 50, key6 = 60;
  int voice;
  voice = voiceMan.noteOnReturnVoice(key1, 100); expectEquals(voice, 0);
  voice = voiceMan.noteOnReturnVoice(key2, 100); expectEquals(voice, 1);
  voice = voiceMan.noteOnReturnVoice(key3, 100); expectEquals(voice, 2);
  voice = voiceMan.noteOnReturnVoice(key4, 100); expectEquals(voice, 3);

  // Now all 4 voices are used up and stealing should take place:
  voice = voiceMan.noteOnReturnVoice(key5, 100); expectEquals(voice, 0);
  voice = voiceMan.noteOnReturnVoice(key6, 100); expectEquals(voice, 1); // fails!

  // Releae voice 3 and trigger another note - for the new note, the juts freed voice should be 
  // used again:
  voice = voiceMan.noteOffReturnVoice(key4);      expectEquals(voice, 3);
  voice = voiceMan.noteOnReturnVoice( key1, 100); expectEquals(voice, 3);


  voiceMan.setNumVoices(2);
  voiceMan.setStealMode(SM::oldest);
  voiceMan.setKillMode(KM::afterSilence);
  voiceMan.reset();

  // noteOn for key1:
  voice = voiceMan.noteOnReturnVoice(key1, 100); 
  expectEquals(voice, 0);
  expectEquals(voiceMan.getNumActiveVoices(),    1);
  expectEquals(voiceMan.getNumReleasingVoices(), 0);

  // noteOff for key1 (releases it, it remains active):
  voice = voiceMan.noteOnReturnVoice(key1,   0); 
  expectEquals(voice, 0);
  expectEquals(voiceMan.getNumActiveVoices(),    1);
  expectEquals(voiceMan.getNumReleasingVoices(), 1);

  // noteOn for key2:
  voice = voiceMan.noteOnReturnVoice(key2, 100); 
  expectEquals(voice, 1);
  expectEquals(voiceMan.getNumActiveVoices(),    2);
  expectEquals(voiceMan.getNumReleasingVoices(), 1);

  // noteOff for key2 (releases it, it remains active):
  voice = voiceMan.noteOnReturnVoice(key2,   0); 
  expectEquals(voice, 1);
  expectEquals(voiceMan.getNumActiveVoices(),    2);
  expectEquals(voiceMan.getNumReleasingVoices(), 2);

  // noteOn for key3, should steal voice 0, so it leaves release mode:
  voice = voiceMan.noteOnReturnVoice(key3, 100); 
  expectEquals(voice, 0);
  expectEquals(voiceMan.getNumActiveVoices(),    2);
  expectEquals(voiceMan.getNumReleasingVoices(), 1);

  // noteOff for key3 (releases it, it remains active):
  voice = voiceMan.noteOnReturnVoice(key3,   0); 
  expectEquals(voice, 0);
  expectEquals(voiceMan.getNumActiveVoices(),    2);
  expectEquals(voiceMan.getNumReleasingVoices(), 2);

  // noteOn for key1 (again) - should steal voice 1:
  voice = voiceMan.noteOnReturnVoice(key1, 100); 
  expectEquals(voice, 1);
  expectEquals(voiceMan.getNumActiveVoices(),    2);
  expectEquals(voiceMan.getNumReleasingVoices(), 1);

  // noteOff for key1 (releases it, it remains active):
  voice = voiceMan.noteOnReturnVoice(key1,   0); 
  expectEquals(voice, 1);
  expectEquals(voiceMan.getNumActiveVoices(),    2);
  expectEquals(voiceMan.getNumReleasingVoices(), 2);

  // noteOn for key1 (again) - should steal voice 0:
  voice = voiceMan.noteOnReturnVoice(key1, 100); 
  expectEquals(voice, 0);
  expectEquals(voiceMan.getNumActiveVoices(),    2);
  expectEquals(voiceMan.getNumReleasingVoices(), 1);

  // noteOff for key1 - should release voice 0 and 1
  voice = voiceMan.noteOnReturnVoice(key1,   0); 
  expectEquals(voice, 0);
  //expectEquals(voice, voiceMan.manyVoices);  // fails - returns the last voice that was released
  expectEquals(voiceMan.getNumActiveVoices(),    2);
  expectEquals(voiceMan.getNumReleasingVoices(), 2);

  // todo: check that both voices are released in a single noteOff in cases where they both play 
  // the same note

  // Trigger 3 different notes, then release them and loop through the killSamples, feeding a value
  // below the threshold, so it should eventually kill the voices:
  voiceMan.reset();
  voice = voiceMan.noteOnReturnVoice(key1, 100);
  voice = voiceMan.noteOnReturnVoice(key2, 100);
  voice = voiceMan.noteOnReturnVoice(key3, 100);
  voice = voiceMan.noteOnReturnVoice(key1,   0);
  voice = voiceMan.noteOnReturnVoice(key2,   0);
  voice = voiceMan.noteOnReturnVoice(key3,   0);

  double* vb = &voiceBuffer[0];
  vb[0] = vb[1] = vb[2] = vb[3] = killThresh;
  bool ok = true;
  for(int n = 1; n < killSamples; n++) {
    voiceMan.findAndKillFinishedVoices();
    active    = voiceMan.getNumActiveVoices();
    releasing = voiceMan.getNumReleasingVoices();
    ok &= active == 2 && releasing == 2; 
    jassert(ok);  }
  expect(ok);
  voiceMan.findAndKillFinishedVoices();         // should find and kill the 2 voices
  active    = voiceMan.getNumActiveVoices();    // ...so this should go to 0
  releasing = voiceMan.getNumReleasingVoices(); // ...and this too
  expectEquals(active,    0);
  expectEquals(releasing, 0);

  int dummy = 0;

  /*
  voiceMan.reset();
  voice = voiceMan.noteOnReturnVoice(key1, 100);
  voice = voiceMan.noteOnReturnVoice(key1, 100);
  voice = voiceMan.noteOnReturnVoice(key1, 100);
  voice = voiceMan.noteOnReturnVoice(key1,   0);
  voice = voiceMan.noteOnReturnVoice(key1,   0);
  voice = voiceMan.noteOnReturnVoice(key1,   0);
  */


  // now the releasingVoices array is overfull (size == 3 where numVoices is only 2) - todo:
  // handle releasingVoices in the same way as activeVoices (resize to maxNumVoices and keep a
  // numReleasingVoices variable






  // ToDo: test voice stealing in the various modes, voice retriggering, etc.



}
// actually, that test belongs into the framework tests - maybe make files UnitTestsAudio/Midi
// and put this test there