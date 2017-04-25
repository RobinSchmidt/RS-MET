/*

AggressorOscSection.h: interface for the AggressorOscSection class.

© Braindoc 2002 (www.braindoc.de)

*/

#if !defined(AggressorOscSection_h_Included)
#define AggressorOscSection_h_Included

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Oscillator.h"
#include "SuperOscillator.h"
#include "LpfHpfApf.h"
#include "PitchEnvRC.h"
#include "ExponentialRamp.h"

class AggressorOscSection
{
public:

 // construction/destruction:
 AggressorOscSection();
 ~AggressorOscSection();

 //parameter settings:
	//general:
 void setSampleRate(double SampleRate);
 void setMonoPoly  (long   MonoPoly);
 void setSlideTime (double SlideTime); 
 void setAccent    (double Accent);
 void setEnvSlope  (double EnvSlope); 
 void setVibDepth  (double VibDepth);
 void setVibRate   (double VibRate);
 void setVibWave   (long   VibWave);
 void setPwmDepth  (double PwmDepth);
 void setPwmRate   (double PwmRate);
 void setPwmWave   (long   PwmWave);

 //oscillator 1:
 void setWaveForm1        (long   WaveForm1);
 void setStartPhase1      (double StartPhase1);
 void setVol1Start        (double Vol1Start);
 void setVol1StartByVel   (double Vol1StartByVel);
 void setVol1Time         (double Vol1Time);
 void setVol1End          (double Vol1End);   
 void setPitch1Start      (double Pitch1Start);
 void setPitch1StartByVel (double Pitch1StartByVel);
 void setPitch1Time       (double Pitch1Time);
 void setPitch1End        (double Pitch1End);  
 void setOsc1Lpf          (double Osc1Lpf);  
 void setOsc1Hpf          (double Osc1Hpf); 
 void setOsc1Apf          (double Osc1Apf); 
 void setPulseWidth1      (double PulseWidth1); 
 void setDensity          (double Density); 
 void setDetuneStart      (double DetuneStart); 
 void setDetuneStartByVel (double DetuneStartByVel); 
 void setDetuneTime       (double DetuneTime); 
 void setDetuneEnd        (double DetuneEnd);
 void setDetuneSpacing    (int    DetuneSpacing);
 void setDetuneRatio      (double DetuneRatio);
 void setDetunePhaseSpread(double DetunePhaseSpread);

 //modulations:
 void setAmBy2Start       (double AmBy2Start);
 void setAmBy2StartByVel  (double AmBy2StartByVel);
 void setAmBy2Time        (double AmBy2Time);
 void setAmBy2End         (double AmBy2End);
 void setAmBy3Start       (double AmBy3Start);
 void setAmBy3StartByVel  (double AmBy3StartByVel);
 void setAmBy3Time        (double AmBy3Time);
 void setAmBy3End         (double AmBy3End);
 void setFmBy1Start       (double FmBy1Start);
 void setFmBy1StartByVel  (double FmBy1StartByVel);
 void setFmBy1Time        (double FmBy1Time);
 void setFmBy1End         (double FmBy1End);
 void setFmBy2Start       (double FmBy2Start);
 void setFmBy2StartByVel  (double FmBy2StartByVel);
 void setFmBy2Time        (double FmBy2Time);
 void setFmBy2End         (double FmBy2End);
 void setFmBy3Start       (double FmBy3Start);
 void setFmBy3StartByVel  (double FmBy3StartByVel);
 void setFmBy3Time        (double FmBy3Time);
 void setFmBy3End         (double FmBy3End);
 void setModulationMode   (long   ModulationMode);
 void setSyncMode         (long   SyncMode);

 //oscillator 2:
 void setWaveForm2        (long   WaveForm2);
 void setStartPhase2      (double StartPhase2);
 void setVol2Start        (double Vol2Start);
 void setVol2StartByVel   (double Vol2StartByVel);
 void setVol2Time         (double Vol2Time);
 void setVol2End          (double Vol2End);   
 void setPitch2Start      (double Pitch2Start);
 void setPitch2StartByVel (double Pitch2StartByVel);
 void setPitch2Time       (double Pitch2Time);
 void setPitch2End        (double Pitch2End);  
 void setOsc2Lpf          (double Osc2Lpf);  
 void setOsc2Hpf          (double Osc2Hpf); 
 void setOsc2Apf          (double Osc2Apf); 
 void setPulseWidth2      (double PulseWidth2); 
 void setTuneCoarse2      (double TuneCoarse2); 
 void setTuneFine2        (double TuneFine2); 

 //oscillator 3:
 void setWaveForm3        (long   WaveForm3);
 void setStartPhase3      (double StartPhase3);
 void setVol3Start        (double Vol3Start);
 void setVol3StartByVel   (double Vol3StartByVel);
 void setVol3Time         (double Vol3Time);
 void setVol3End          (double Vol3End);   
 void setPitch3Start      (double Pitch3Start);
 void setPitch3StartByVel (double Pitch3StartByVel);
 void setPitch3Time       (double Pitch3Time);
 void setPitch3End        (double Pitch3End);  
 void setOsc3Lpf          (double Osc3Lpf);  
 void setOsc3Hpf          (double Osc3Hpf); 
 void setOsc3Apf          (double Osc3Apf); 
 void setPulseWidth3      (double PulseWidth3); 
 void setTuneCoarse3      (double TuneCoarse3); 
 void setTuneFine3        (double TuneFine3);

 //ringmodulation output volumes:
 void setRm12Start        (double Rm12Start);
 void setRm12StartByVel   (double Rm12StartByVel);
 void setRm12Time         (double Rm12Time);
 void setRm12End          (double Rm12End);
 void setRm13Start        (double Rm13Start);
 void setRm13StartByVel   (double Rm13StartByVel);
 void setRm13Time         (double Rm13Time);
 void setRm13End          (double Rm13End);
 void setRm23Start        (double Rm23Start);
 void setRm23StartByVel   (double Rm23StartByVel);
 void setRm23Time         (double Rm23Time);
 void setRm23End          (double Rm23End);

	//event processing:
 void triggerNote  (long   NoteNumber, long Velocity, long Detune);
 void slideToNote  (long   NoteNumber, long Velocity, long Detune);
 void noteOff      (long   NoteNumber); 
	void setPitchBend (double PitchBend);

 // audio processing:
 __forceinline double getSample();

	// others:
	void   resetOscillators();  //resets oscillators to their start-phases
                     
protected:

 // some switch variables:
 bool   osc1IsPlaying;     // indicate, if the osc is playing
 bool   osc2IsPlaying;
 bool   osc3IsPlaying;
 bool   vibIsActive;
 bool   pwmIsActive;
 bool   modulationsOff;    // true, if modulationMode==0
 bool   syncIsActive;

	// tuning of oscillator 2 and 3:
 doubleA tuneCoarse2, tuneFine2, tuneCoarse3, tuneFine3;
 doubleA tuneFactor2, tuneFactor3;

	// pulsewidth of oscillator 2 and 3:
	doubleA pulseWidth1, pulseWidth2, pulseWidth3;

 // parameter variables:
	doubleA sampleRate;
 doubleA slideTime;       // slide-time in seconds
 doubleA slideSamples;    // slide-time in samples
 doubleA pitchIncPerSamp; // increment for the pitch each sample
 doubleA freqFacPerSamp;  // factor, with which the osc's frequency has to
                          // be multiplied each sample until the target
                          // frequency is reached (which happens after
                          // slideSamples Samples)
	doubleA accent;          // scales the velocity dependence of the ramp
                          // start values
 doubleA vibDepth;        // depth of vibrato (amplitude of the lfo)
 doubleA pwmDepth;        // depth of pulse-width modulation

	doubleA vol1Start;       // start volume of osc1 (has to be scaled here
                          // according to velocity vol1StartByVel and accent
	doubleA vol1StartByVel;  // the velocity dependence of osc1's start volume
	doubleA vol1StartByVelAc;// the same, but scaled by the accent parameter

	// analogous values for osc1's pitch envelope:
	doubleA phase1, pitch1Start, pitch1StartByVel, pitch1StartByVelAc;

	// modulation envelope parameters:
	doubleA amBy2Start, amBy2StartByVel, amBy2StartByVelAc,
		       amBy3Start, amBy3StartByVel, amBy3StartByVelAc,
		       fmBy1Start, fmBy1StartByVel, fmBy1StartByVelAc,
		       fmBy2Start, fmBy2StartByVel, fmBy2StartByVelAc,
		       fmBy3Start, fmBy3StartByVel, fmBy3StartByVelAc;

	// the same for oscillator 2 and 3:
	doubleA phase2, vol2Start, vol2StartByVel, vol2StartByVelAc,
         pitch2Start, pitch2StartByVel, pitch2StartByVelAc;

	doubleA phase3, vol3Start, vol3StartByVel, vol3StartByVelAc,
         pitch3Start, pitch3StartByVel, pitch3StartByVelAc;

	// sync- and modulation-mode:
	int syncMode;
	int modulationMode;

	 // ringmod:
 doubleA rm12Start, rm12StartByVel, rm12StartByVelAc,
         rm13Start, rm13StartByVel, rm13StartByVelAc,
         rm23Start, rm23StartByVel, rm23StartByVelAc;

	// not needed anymore
 //doubleA pitchEnvPeak;      //needs to be stored here in order to scale it by note velocity
 //doubleA pitchEnvPeakByVel;

 int     currentNote;       //most recently triggered note which is on (-1 if none)
 int     currentDetune;     //detune of the current note
 doubleA currentPitch;      //actual pitch of the osc at this moment (even when in slide)
 doubleA targetPitch;       //pitch, which we want to slide to
	doubleA freqBend;          //scale factor for the frequency in response to setPitchBend
 doubleA currentFreq;       //actual frequency of the osc at this moment (even when in slide)
 doubleA targetFreq;        //frequency, which we want to slide to
 int     slideSampCount;    //counter for the samples in which slide already took place

 //embedded audio-module objects:
	//the oscillators (LFO's included):
 SuperOscillator osc1;
	Oscillator      osc2, osc3, vibLfo, pwmLfo;

	// the per-oscillator filters:
 LpfHpfApf       flt1, flt2, flt3;

 // the ramp-envelopes:
	ExponentialRamp ampRamp1, ampRamp2, ampRamp3,    //amplitudes of the oscillators
                 freqRamp1, freqRamp2, freqRamp3, //frequenceies of the oscillators
                 rm12Ramp, rm13Ramp, rm23Ramp;    //ringmodulation volume ramps

	ExponentialRamp detuneRamp;

	ExponentialRamp amBy2Ramp, amBy3Ramp,
		               fmBy1Ramp, fmBy2Ramp, fmBy3Ramp;



 //sample          prevOutputOsc1; //remember previous output of osc 1 for feedback-FM
                                   //...is realzied via a stsic variable in getSample() function

	//doubleA          fm12Index; //for testing

};

//-----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):
__forceinline double AggressorOscSection::getSample()
{
 static doubleA instFreq1, instFreq2, instFreq3; //frequency at this instant of time (for the 3 oscs)

	static doubleA ampFactor1, ampFactor2, ampFactor3; 
	//factor for the amplitudes (comes from envelope and amplitude modulation)

	static doubleA vibValue;

	static doubleA instPulseWidth1, instPulseWidth2, instPulseWidth3, pulseWidthOffset;

	static doubleA freqDeviation; //for FM in osc1

 static doubleA masterPhase, masterIncrement, timeStamp; 
  // for correctly setting up the phases in the sync-slaves

	 //the audio signals (at different satges of the signal chain):
 static doubleA osc1Out,              osc2Out,              osc3Out,
                osc1WithFilter,       osc2WithFilter,       osc3WithFilter,
                osc1WithFilterAndEnv, osc2WithFilterAndEnv, osc3WithFilterAndEnv;

 static doubleA rm12Out, rm13Out, rm23Out;
 rm12Out = 0.0;
 rm13Out = 0.0;
 rm23Out = 0.0;

	static doubleA outputSample;

 //update the frequency (slide):
 if(currentFreq!=targetFreq)
 { 
  currentPitch += pitchIncPerSamp;
  currentFreq  *= freqFacPerSamp;
  slideSampCount++;
 }
 //compensate for roundoff-error in frequency after slide:
 if(slideSampCount>=slideSamples)
 {
  currentFreq  = targetFreq;
  currentPitch = targetPitch;
 }

 //calculate freq-factor for vibrato:
 if(vibIsActive)
	 vibValue = vibLfo.getSample()*vibDepth;  //linear vibrato would be more efficient
 else
  vibValue = 0.0;

 // calculate an offset for the pulseWidths from the PWM-LFO:
 if(pwmIsActive)
	 pulseWidthOffset = pwmDepth*pwmLfo.getSample();  
 else
  pulseWidthOffset = 0.0;

 if(osc2IsPlaying)
 {
  // calculate instantaneous frequency of osc2:
  instFreq2  = freqRamp2.getSample() * currentFreq * tuneFactor2 * freqBend;
  instFreq2 += vibValue*instFreq2;

  /*
	 instFreq2  = currentFreq*freqRamp2.getSample(); // applies by the pitch ramp
  instFreq2 *= tuneFactor2;                       // applies tuning-factor
	 instFreq2 *= freqBend;                          // applies pitchbend
  instFreq2 += vibValue*instFreq2;                // applies the vibrato
  */

  // calculate instantaneous pulse-width of osc2:
	 instPulseWidth2 = pulseWidth2 + pulseWidthOffset;

  // calculate instantaneous amplitude of osc2:
	 ampFactor2 = ampRamp2.getSample();      

  // set up the oscillator:
  osc2.setFreq(instFreq2);
	 osc2.setPulseWidth(instPulseWidth2);
  osc2.calcIncrements();

  // calculate the oscillators output:
  osc2Out = osc2.getSample();

  // apply the static filter:
	 osc2WithFilter = flt2.getSample(osc2Out);

  // apply the amp-ramp:
  osc2WithFilterAndEnv = osc2WithFilter*ampFactor2;
 }
 else
 {
  osc2Out              = 0.0;
  osc2WithFilter       = 0.0;
  osc2WithFilterAndEnv = 0.0;
 }

 if(osc3IsPlaying)
 {
  // calculate instantaneous frequency of osc2:
  instFreq3  = freqRamp3.getSample() * currentFreq * tuneFactor3 * freqBend;
  instFreq3 += vibValue*instFreq3;

  /*
	 instFreq3  = currentFreq*freqRamp3.getSample();   //applies by the pitch ramp:
  instFreq3 *= tuneFactor3;                         //applies tuning-factor:
	 instFreq3 *= freqBend;                            //applies pitchbend
  instFreq3 += vibValue*instFreq3;                  //applies the vibrato
  */

  // calculate instantaneous pulse-width of osc3:
	 instPulseWidth3 = pulseWidth3 + pulseWidthOffset;

  // calculate instantaneous amplitude of osc3:
	 ampFactor3 = ampRamp3.getSample();   

  // set up the oscillator:
  osc3.setFreq(instFreq3);
	 osc3.setPulseWidth(instPulseWidth3);
  osc3.calcIncrements();

  // calculate the oscillators output:
  osc3Out = osc3.getSample();

  // apply the static filter:
	 osc3WithFilter = flt3.getSample(osc3Out);

  // apply the amp-ramp:
  osc3WithFilterAndEnv = osc3WithFilter*ampFactor3;
 }
 else
 {
  osc3Out              = 0.0;
  osc3WithFilter       = 0.0;
  osc3WithFilterAndEnv = 0.0;
 }


 if(osc1IsPlaying)
 {
	 instFreq1  = currentFreq*freqRamp1.getSample();   //applies by the pitch ramp
	 instFreq1 *= freqBend;                            //applies pitchbend
  instFreq1 += vibValue*instFreq1;                  //applies the vibrato

	 instPulseWidth1 = pulseWidth1 + pulseWidthOffset; //applies pwm to the static pulseWidth

	 ampFactor1 = ampRamp1.getSample();                //calculates the output amplitude

  //apply the modulations:
  if(modulationsOff)
  {
   //do nothing
  }
 	else if(modulationMode==1)  //the oscillator outputs directly modulate osc1
 	{
	 	//feedback-FM doesn't work yet
   //calculate frequency deviation of osc1 caused by osc1 (feedback-FM):
	  //freqDeviation = fmBy1Ramp.getSample()*instFreq1; 
	  //use the previous output of osc1 (it keeps its value from one call to another because it is
   //declared as a static variable) to modulate osc's frequency:
	  //instFreq1 += fmBy1Ramp.getSample() * osc1Out * sampleRate * osc1.tableLengthRec; 

   //calculate frequency deviation of osc1 caused by osc2 :
	  freqDeviation = fmBy2Ramp.getSample()*instFreq2; 
	  //use also the previous output of osc2 :
	  instFreq1 += freqDeviation * osc2Out;

	  //calculate frequency deviation of osc1 caused by osc3 :
	  freqDeviation = fmBy3Ramp.getSample()*instFreq3; 
	  //use also the previous output of osc3 :
	  instFreq1 += freqDeviation * osc3Out;

   //calculate the freq-deviation caused by feedback-fm:
   instFreq1 = instFreq1 + fmBy1Ramp.getSample()*instFreq1*osc1Out; //ranges from 0 to 2*instFreq1
   //does not work this way -> pitch decreases when incresing amount
   //osc1.modulateIncrement(fmBy1Ramp.getSample()*osc1Out);

		 //calculate amplitude factor for osc1 from its envelope and the amplitude modulation ramps:
		 ampFactor1 = ampFactor1 * (1 + amBy2Ramp.getSample()*osc2Out);
		 ampFactor1 = ampFactor1 * (1 + amBy3Ramp.getSample()*osc3Out);

   //calculate the ringmodulations output values:
   rm12Out = osc1Out*osc2Out*rm12Ramp.getSample();
   rm13Out = osc1Out*osc3Out*rm13Ramp.getSample();
   rm23Out = osc2Out*osc3Out*rm23Ramp.getSample();
	 }
	 else if(modulationMode==2)
	 {
   //calculate frequency deviation of osc1 caused by osc2 :
	  freqDeviation = fmBy2Ramp.getSample()*instFreq2; 
	  //use also the previous output of osc2 :
	  instFreq1 += freqDeviation * osc2WithFilter;

	  //calculate frequency deviation of osc1 caused by osc3 :
	  freqDeviation = fmBy3Ramp.getSample()*instFreq3; 
	  //use also the previous output of osc3 :
	  instFreq1 += freqDeviation * osc3WithFilter;

   //calculate the freq-deviation caused by feedback-fm:
   instFreq1 = instFreq1 + fmBy1Ramp.getSample()*instFreq1*osc1WithFilter; //ranges from 0 to 2*instFreq1

		 //calculate amplitude factor for osc1 from its envelope and the amplitude modulation ramps:
		 ampFactor1 = ampFactor1 * (1 + amBy2Ramp.getSample()*osc2WithFilter);
		 ampFactor1 = ampFactor1 * (1 + amBy3Ramp.getSample()*osc3WithFilter);

   //calculate the ringmodulations output values:
   rm12Out = osc1WithFilter*osc2WithFilter*rm12Ramp.getSample();
   rm13Out = osc1WithFilter*osc3WithFilter*rm13Ramp.getSample();
   rm23Out = osc2WithFilter*osc3WithFilter*rm23Ramp.getSample();
	 }
	 else if(modulationMode==3)
	 {
   //calculate frequency deviation of osc1 caused by osc2 :
	  freqDeviation = fmBy2Ramp.getSample()*instFreq2; 
	  //use also the previous output of osc2 :
	  instFreq1 += freqDeviation * osc2WithFilterAndEnv;

	  //calculate frequency deviation of osc1 caused by osc3 :
	  freqDeviation = fmBy3Ramp.getSample()*instFreq3; 
	  //use also the previous output of osc3 :
	  instFreq1 += freqDeviation * osc3WithFilterAndEnv;

		 //calculate amplitude factor for osc1 from its envelope and the amplitude modulation ramps:
		 ampFactor1 = ampFactor1 * (1 + amBy2Ramp.getSample()*osc2WithFilterAndEnv);
		 ampFactor1 = ampFactor1 * (1 + amBy3Ramp.getSample()*osc3WithFilterAndEnv);

   //calculate the ringmodulations output values:
   rm12Out = osc1WithFilterAndEnv*osc2WithFilterAndEnv*rm12Ramp.getSample();
   rm13Out = osc1WithFilterAndEnv*osc3WithFilterAndEnv*rm13Ramp.getSample();
   rm23Out = osc2WithFilterAndEnv*osc3WithFilterAndEnv*rm23Ramp.getSample();
	 }

  //set up the oscillator:
	 osc1.setDetune(detuneRamp.getSample());
	 osc1.setPulseWidth(instPulseWidth1);   //the function-call makes sure that pw is between 0.01 and 0.99
  osc1.setFreq(instFreq1);               //optimize here!
  //osc1.calcIncrements();

  //calculate the oscillators output:
  //osc1Out = osc1.getSample();
  osc1Out = osc1.getSample();

  //apply the static filter:
	 osc1WithFilter = flt1.getSample(osc1Out);

  //apply the amp-ramp:
  osc1WithFilterAndEnv = osc1WithFilter*ampFactor1;
 }
 else
 {
  osc1Out              = 0.0;
  osc1WithFilter       = 0.0;
  osc1WithFilterAndEnv = 0.0;
 }

	// sync if sync is active:
 if(syncIsActive)
 {
  masterPhase     = osc1.phaseIndices[0];
  masterIncrement = osc1.getIncrement();
  timeStamp       = masterPhase/masterIncrement;
	 if( osc1.wraparoundOccurred && (syncMode==1 || syncMode==3) )
  {
	 	//osc2.setPhase(osc1.phaseIndices[0]); // for some old tracks
   osc2.triggerSync(timeStamp);
  }
	 if( osc1.wraparoundOccurred && (syncMode==2 || syncMode==3) )
  {
   osc3.triggerSync(timeStamp);
  }
 }

 //add the signals:
	outputSample =   osc1WithFilterAndEnv + osc2WithFilterAndEnv + osc3WithFilterAndEnv 
		              + rm12Out + rm13Out + rm23Out;

 return outputSample;
}


#endif // !defined(AggressorOscSection_h_Included)
