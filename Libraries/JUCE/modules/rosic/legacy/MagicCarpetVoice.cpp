#include "MagicCarpetVoice.h"

//----------------------------------------------------------------------------
// construction/destruction:

MagicCarpetVoice::MagicCarpetVoice()
{
 // init member variables:
 currentNote       = -1; 
 currentNoteVel    = 0; // velocity = 0 (-> note is off)
 currentNoteDetune = 0;
 currentNoteAge    = 0;
 isPlaying         = true;

 //topLeftDetune     = 0.0;
 //topRightDetune    = 0.0;
 //bottomLeftDetune  = 0.0;
 //bottomRightDetune = 0.0;

 topLeftAmplitude       = 0.5;
 topRightAmplitude      = 0.5;
 bottomLeftAmplitude    = 0.5;
 bottomRightAmplitude   = 0.5;

 topLeftLevelByKey      = 0.0;
 topRightLevelByKey     = 0.0;
 bottomLeftLevelByKey   = 0.0;
 bottomRightLevelByKey  = 0.0;
 topLeftTuneByKey       = 100.0;
 topRightTuneByKey      = 100.0;
 bottomLeftTuneByKey    = 100.0;
 bottomRightTuneByKey   = 100.0;

 fltFreq                = 1000.0;
 fltFreqByKey           = 0.0;
 fltFreqByVel           = 0.0;
 fltFreqScaler          = 1.0;


 ampEnvDurByKey         = 0.0;
 fltEnvDurByKey         = 0.0;

 topLeftLevelByVel      = 0.0;
 topRightLevelByVel     = 0.0;
 bottomLeftLevelByVel   = 0.0;
 bottomRightLevelByVel  = 0.0;
 ampByVel               = 0.0;

 ampEnvPeakByVel        = 0.0;
 fltEnvPeakByVel        = 0.0;

 sampleRate        = 44100.0;    

 x1L               = 0.0;
 x1R               = 0.0;

 // init embedded modules:
 // ...
 ampEnv.setStart(-96.0);
 ampEnv.setAttack(0.5);
 ampEnv.setPeak(0.0);
 ampEnv.setDecay(2.0);
 ampEnv.setSustain(0.0);
 ampEnv.setRelease(2.5);
 ampEnv.setTauScale(0.2);

 //filter.setTwoStages(false);
}

MagicCarpetVoice::~MagicCarpetVoice()
{
 
}

//----------------------------------------------------------------------------
// parameter settings: 

void MagicCarpetVoice::setSampleRate(double newSampleRate)
{
 if(newSampleRate > 0.1)
  sampleRate = newSampleRate;

 topLeftSource.setSampleRate(sampleRate);
 topRightSource.setSampleRate(sampleRate);
 bottomLeftSource.setSampleRate(sampleRate);
 bottomRightSource.setSampleRate(sampleRate);

 // tell the LFOs the sample-rate.....

 fltEnv.setSampleRate(sampleRate);
 filter.setSampleRate(sampleRate);

 ampEnv.setSampleRate(sampleRate);
}

/*
void MagicCarpetVoice::setAmpEnvAttack(double newAmpEnvAttack)
{
 ampEnv.setAttack(newEnvAttack);
}


//....
*/

//----------------------------------------------------------------------------
// event processing:

void MagicCarpetVoice::noteOn(int newNoteNumber, int newVelocity)
{
 currentNote       = newNoteNumber;
 currentNoteVel    = newVelocity;
 // currentNoteDetune = Detune;

 // handle note-off events:
 if( currentNoteVel==0 )
 {
  ampEnv.noteOff();   // trigger release on note-off events (vel=0)
  fltEnv.noteOff();  
  return;
 }

 //else:

 // reset the note-age counter:
 currentNoteAge = 0;

 // calculate and set the amplitude-factors of the individual slots according
 // to the key and velocity:
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

 // calculate frequencies (not including detuning) for the sources from the
 // note and set up the sources acordingly, thereby take into account the
 // keytracking of the tune-parameter (100% means normal chromatic playing, 
 // 0% means that the sample is not to be transposed at all)
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

 /*
 // this was the simple version wihtout the keytracking-parameter:
 freq = pitchToFreq(currentNote);
 topLeftSource.setFrequency(freq);
 topRightSource.setFrequency(freq);
 bottomLeftSource.setFrequency(freq);
 bottomRightSource.setFrequency(freq);
 */

 // calculate and set the time-scale-factors for the envelopes:
 double timeScaler;
 timeScaler  = pow(2.0, 0.01*ampEnvDurByKey*(newNoteNumber-60.0)/12.0);
 timeScaler *= pow(2.0, 0.01*ampEnvDurByVel*(newVelocity-64.0)/63.0);
 ampEnv.setTimeScale(timeScaler);

 timeScaler  = pow(2.0, 0.01*fltEnvDurByKey*(newNoteNumber-60.0)/12.0);
 timeScaler *= pow(2.0, 0.01*fltEnvDurByVel*(newVelocity-64.0)/63.0);
 fltEnv.setTimeScale(timeScaler);

 // calculate and set the peak-scale-factors for the envelopes:
 double peakScaler;
 peakScaler = DB2AMP(ampEnvPeakByVel*(newVelocity-64.0)/63.0);
 ampEnv.setPeakScale(peakScaler);

 peakScaler = pitchOffsetToFreqFactor(fltEnvPeakByVel*(newVelocity-64.0)/63.0);
 fltEnv.setPeakScale(peakScaler);

 // calculate a the key- and velocity-sclaed filter frequency:
 //fltFreqScaled  = fltFreq;
 fltFreqScaler  = pow(2.0, (0.01*fltFreqByKey*(newNoteNumber-60.0)/12.0));
 fltFreqScaler *= pow(2.0, (0.01*fltFreqByVel*(newVelocity-64.0)/63.0));

 if(ampEnv.endIsReached())
 {
  resetDifferentiator();
  filter.resetBuffers();

  // when we do a reset of the wavetable, we need to calculate some 
  // dummy source-samples (about 2000) in order to let the sources 
  // lowpass/highpass filters setlle to a steady state (warming up):
  double dummyL, dummyR;

  // reset the phase-pointers of the sources (with offset of -2000):
  topLeftSource.reset(-2000.0);
  topRightSource.reset(-2000.0);
  bottomLeftSource.reset(-2000.0);
  bottomRightSource.reset(-2000.0);

  for(int i=0; i<2000; i++)
   getSourceSampleFrameStereo(&dummyL, &dummyR);

  // O.K. now the source-filters are warmed up, we are now ready to calculate
  // actual output-signals....
 }
 else // the voice has not died out, so we don't retrigger the oscillators
 {
  topLeftSource.decrementPhase();
  topRightSource.decrementPhase();
  bottomLeftSource.decrementPhase();
  bottomRightSource.decrementPhase();

  double dummyL, dummyR;
  getSourceSampleFrameStereo(&dummyL, &dummyR); 
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
 if( !topLeftSource.isInLoopMode() )
  topLeftSource.reset();
 if( !topRightSource.isInLoopMode() )
  topRightSource.reset();
 if( !bottomLeftSource.isInLoopMode() )
  bottomLeftSource.reset();
 if( !bottomRightSource.isInLoopMode() )
  bottomRightSource.reset();

 // trigger the envelopes:
 ampEnv.trigger(true);
 fltEnv.trigger();
}


//----------------------------------------------------------------------------
// others:

void MagicCarpetVoice::resetDifferentiator()
{
 x1L = 0.0;
 x1R = 0.0;
}


