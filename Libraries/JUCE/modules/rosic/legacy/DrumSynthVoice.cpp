#include "DrumSynthVoice.h"

//----------------------------------------------------------------------------
// construction/destruction:

DrumSynthVoice::DrumSynthVoice()
{
 sampleRate        = 44100.0;

 // init member variables:
 currentNote       = -1; 
 currentNoteVel    = 0; // velocity = 0 (-> note is off)
 currentNoteDetune = 0;
 currentNoteAge    = 0;
 isPlaying         = true;

 osc1Amplitude     = 1.0;
 osc1TuneByKey     = 100.0;
 osc1LevelByVel    = 0.0;

 //osc2Amplitude     = 1.0;
 //osc2TuneByKey     = 100.0;
 //osc2ByVel         = 0.0;

 fltFreq           = 1000.0;
 fltFreqByKey      = 0.0;
 fltFreqByVel      = 0.0;
 fltFreqScaler     = 1.0;


 // initialize the amplitude envelope:
 ampEnv.setSampleRate(sampleRate);
 ampEnv.insertBreakpoint(0.2, 1.0, 3, 1.0); // peak-value, index==1
 ampEnv.insertBreakpoint(0.4, 0.5, 3, 1.0); // loop-start, index==2
 ampEnv.insertBreakpoint(0.6, 0.6, 2, 1.0); // loop-in-between, 3
 ampEnv.insertBreakpoint(1.0, 0.5, 2, 1.0); // loop-end, 4

 ampEnv.modifyBreakpoint(ampEnv.lastBreakpointIndex(), 0.1, 0.0, 3, 1.0);

 //ampEnv.setLoopEndIndex(ampEnv.lastBreakpointIndex()-1);
 ampEnv.setLoopEndIndex(4);
 ampEnv.setLoopStartIndex(2);
 ampEnv.setLoopMode(true);

 // insert a breakpoint inside the loop:
 //ampEnv.insertBreakpoint(0.75, 0.2, 3, 1.0); // loop-start, index==2

}

DrumSynthVoice::~DrumSynthVoice()
{
 
}

//----------------------------------------------------------------------------
// parameter settings: 

void DrumSynthVoice::setSampleRate(double newSampleRate)
{
 if(newSampleRate > 0.1)
  sampleRate = newSampleRate;

 osc1.setSampleRate(sampleRate);
 osc2.setSampleRate(sampleRate);
 ampEnv.setSampleRate(sampleRate);
}

//----------------------------------------------------------------------------
// event processing:

void DrumSynthVoice::noteOn(int newNoteNumber, int newVelocity)
{
 currentNote       = newNoteNumber;
 currentNoteVel    = newVelocity;
 // currentNoteDetune = Detune;

 // handle note-off events:
 if( currentNoteVel==0 )
 {
  ampEnv.noteOff();   // trigger release on note-off events (vel=0)
  //fltEnv.noteOff();  
  return;
 }

 //else:

 // reset the note-age counter:
 currentNoteAge = 0;

 // calculate and set the amplitude-factors of the individual slots according
 // to the key and velocity:
 /*
 double ampScaler;
 ampScaler  = dB2amp(topLeftLevelByKey*(newNoteNumber-64.0)/63.0);
 ampScaler *= dB2amp(topLeftLevelByVel*(newVelocity-64.0)/63.0);
 topLeftSource.setAmpScaler(ampScaler);
 ampScaler  = dB2amp(topRightLevelByKey*(newNoteNumber-64.0)/63.0);
 ampScaler *= dB2amp(topRightLevelByVel*(newVelocity-64.0)/63.0);
 topRightSource.setAmpScaler(ampScaler);
 ampScaler  = dB2amp(bottomLeftLevelByKey*(newNoteNumber-64.0)/63.0);
 ampScaler *= dB2amp(bottomLeftLevelByVel*(newVelocity-64.0)/63.0);
 bottomLeftSource.setAmpScaler(ampScaler);
 ampScaler  = dB2amp(bottomRightLevelByKey*(newNoteNumber-64.0)/63.0);
 ampScaler *= dB2amp(bottomRightLevelByVel*(newVelocity-64.0)/63.0);
 bottomRightSource.setAmpScaler(ampScaler);
 */

 // calculate frequencies (not including detuning) for the sources from the
 // note and set up the sources acordingly, thereby take into account the
 // keytracking of the tune-parameter (100% means normal chromatic playing, 
 // 0% means that the sample is not to be transposed at all)
 /*
 double rootKey, pitch, freq;
 rootKey = freqToPitch(topLeftSource.fundamentalFreq);
 pitch   = rootKey + (0.01*topLeftTuneByKey)*(newNoteNumber-rootKey);
 freq    = pitchToFreq(pitch);
 topLeftSource.setFrequency(freq);

 rootKey = freqToPitch(topRightSource.fundamentalFreq);
 pitch   = rootKey + (0.01*topRightTuneByKey)*(newNoteNumber-rootKey);
 freq    = pitchToFreq(pitch);
 topRightSource.setFrequency(freq);

 rootKey = freqToPitch(bottomLeftSource.fundamentalFreq);
 pitch   = rootKey + (0.01*bottomLeftTuneByKey)*(newNoteNumber-rootKey);
 freq    = pitchToFreq(pitch);
 bottomLeftSource.setFrequency(freq);

 rootKey = freqToPitch(bottomRightSource.fundamentalFreq);
 pitch   = rootKey + (0.01*bottomRightTuneByKey)*(newNoteNumber-rootKey);
 freq    = pitchToFreq(pitch);
 bottomRightSource.setFrequency(freq);
 */

 // calculate and set the time-scale-factors for the envelopes:
 /*
 double timeScaler;
 timeScaler  = pow(2.0, 0.01*ampEnvDurByKey*(newNoteNumber-60.0)/12.0);
 timeScaler *= pow(2.0, 0.01*ampEnvDurByVel*(newVelocity-64.0)/63.0);
 ampEnv.setTimeScale(timeScaler);
 */

 /*
 timeScaler  = pow(2.0, 0.01*fltEnvDurByKey*(newNoteNumber-60.0)/12.0);
 timeScaler *= pow(2.0, 0.01*fltEnvDurByVel*(newVelocity-64.0)/63.0);
 fltEnv.setTimeScale(timeScaler);
 */

 // calculate and set the peak-scale-factors for the envelopes:
 /*
 double peakScaler;
 peakScaler = DB2AMP(ampEnvPeakByVel*(newVelocity-64.0)/63.0);
 ampEnv.setPeakScale(peakScaler);
 */

 /*
 peakScaler = pitchOffsetToFreqFactor(fltEnvPeakByVel*(newVelocity-64.0)/63.0);
 fltEnv.setPeakScale(peakScaler);
 */

 // calculate a the key- and velocity-sclaed filter frequency:
 /*
 fltFreqScaled  = fltFreq;
 fltFreqScaler  = pow(2.0, (0.01*fltFreqByKey*(newNoteNumber-60.0)/12.0));
 fltFreqScaler *= pow(2.0, (0.01*fltFreqByVel*(newVelocity-64.0)/63.0));
 */

 if(ampEnv.endIsReached)
 {
  // when we do a reset of the wavetable, we need to calculate some 
  // dummy source-samples (about 2000) in order to let the sources 
  // lowpass/highpass filters setlle to a steady state (warming up):
  double dummyL, dummyR;

  // reset the phase-pointers of the sources (with offset of -2000):
  osc1.reset(-2000.0);
  osc2.reset(-2000.0);

  for(int i=0; i<2000; i++)
  {
   osc1.getSampleFrameStereo(&dummyL, &dummyR);
   osc2.getSampleFrameStereo(&dummyL, &dummyR);
  }

  // O.K. now the source-filters are warmed up, we are now ready to calculate
  // actual output-signals....
 }
 else // the voice has not died out, so we don't retrigger the oscillators
 {
  osc1.decrementPhase();
  osc2.decrementPhase();

  double dummyL, dummyR;
  osc1.getSampleFrameStereo(&dummyL, &dummyR);
  osc2.getSampleFrameStereo(&dummyL, &dummyR);
   // When the new note has a velocity different from the old one and a source 
   // volume is velocity dependent, a step in the source-signal would result. 
   // This step in the source signal would be converted into an impulse by the
   // differentiator. By calling the getSourceSample() once, we can avoid that
   // impulse at the expense of dropping one sample from the output-stream 
   // (which is much less audible). We could compensate this by decrementing
   // the phase-pointer of all sources first, but this would screw up their
   // internal filters states a bit.
 }

 // reset the phase-pointers of the sources, which are not in loop-mode 
 // independently from the envelope:
 if( !osc1.isInLoopMode() )
  osc1.reset();
 if( !osc2.isInLoopMode() )
  osc2.reset();

 // trigger the envelopes:
 ampEnv.noteOn(true);
 //fltEnv.trigger();
}


//----------------------------------------------------------------------------
// others:



