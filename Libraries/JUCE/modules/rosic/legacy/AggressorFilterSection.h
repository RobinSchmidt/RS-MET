/*

AggressorFilterSection.h: interface for the AggressorFilterSection class.

© Braindoc 2002 (www.braindoc.de)

This is a monophonic filter section with an exponential-segment envelope. 

*/

#if !defined(AggressorFilterSection_h_Included)
#define AggressorFilterSection_h_Included

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "CookbookFilter.h"
#include "MoogFilter.h"
#include "PitchEnvRc.h"

const double refCutoff    = 200.0;          //reference cutoff frequency for resonance-scaling by cutoff
const double recRefCutoff = 1.0/refCutoff;  //reciprocal

class AggressorFilterSection
{
public:
 //construction/destruction:
 AggressorFilterSection();
 virtual ~AggressorFilterSection();
 
 //parameter settings:
 void triggerNote        (long   NoteNumber, long Velocity, long Detune);
 void slideToNote        (long   NoteNumber, long Velocity, long Detune);
 void noteOff            (long   NoteNumber); 
	void setPitchBend       (double PitchBend);
 void setSampleRate      (double SampleRate);
 void setSlideTime       (double SlideTime);
 void setEnvSlope        (double EnvSlope);

 void setMode            (long   Mode);
 void setStages          (long   Stages);

 void setDryWet          (double DryWet);      //value is assumed to be in % wet
 void setDryWetByVel     (double DryWetByVel);

 void setCutoff          (double Cutoff);
 void setCutoffByVel     (double CutoffByVel);
	void setCutoffByInput   (double CutoffByInput);
	void setCutoffByOutput  (double CutoffByOutput);

 void setReso            (double Reso);
 void setResoByVel       (double ResoByVel);   //in "dB"
	void setResoByInput     (double ResoByInput);
	void setResoByOutput    (double ResoByOutput);

 void setGain            (double Gain);        //gain for peaking and shelfing filters
 void setResoByCutoff    (double ResoByCut);
 void setKeyTrack        (double KeyTrack);
	void setDrive           (double Drive);
	void setDcOffset        (double DcOffset);
 void setRefKey          (long   RefKey);      //reference key at which keytracking is neutral
 void setOverSamp        (long   OverSamp);

 //filter envelope settings:
 void setFiltEnvStart    (double Start);
 void setFiltEnvAttack   (double Attack);
 void setFiltEnvPeak     (double Peak);
 void setFiltEnvPeakByVel(double PeakByVel);
 void setFiltEnvHold     (double Hold);
 void setFiltEnvDecay    (double Decay);
 void setFiltEnvRelease  (double Release);
 void setFiltEnvEnd      (double End);
 void setFiltEnvDurByVel (double EnvDurByVel);

 //audio processing:
 __forceinline double getSample(double In);

 //others:
	void resetBuffers();  //set the memorized samples in the filter opcodes to zero
 
protected:

 // embedded audio-modules:
 CookbookFilter biQuadCasc;
 MoogFilter     moogFilter,fatFilter;
 PitchEnvRc     filtEnv;

 // parameter variables:
	doubleA sampleRate;      // the sample-rate
 doubleA slideTime;       // slide-time in seconds
 doubleA slideSamples;    // slide-time in samples
 intA    slideSampCount;  // counter for the samples in which slide already
                          // took place
 doubleA pitchIncPerSamp; // increment for the pitch each sample
 doubleA freqFacPerSamp;  // factor, with which the osc's frequency has to be
                          // multiplied each sample until the target frequency
                          // is reached (which happens after slideSamples
                          // samples)

 intA    mode;            // mode has to be stored here to switch between the
                          // filter classes

 doubleA dryWet,
         dryWetByVel,
         dryVolume, 
         wetVolume;

 doubleA cutoff;
 doubleA cutoffByVel;       // velocity scales cutoff by +-cutoffByVel 
                            // semitones (+-0 @ velo=64)
 doubleA cutoffByInput;
 doubleA cutoffByOutput;

 doubleA reso;
 doubleA resoByVel;         // is in "dB"
	doubleA resoByInput;
	doubleA resoByOutput;
 doubleA resoByCutoff;
 doubleA resoScaleExponent; // a value calculated from resoByCutoff to scale
                            // the resonance
	doubleA drive;
 doubleA dcOffset;

 doubleA keyTrack;
 intA    refKey;

 doubleA filtEnvPeak;
 doubleA filtEnvPeakByVel, filtEnvDurByVel;

 // status variables:
 intA    currentNote;    // most recently triggered note which is on (-1 if none)
 intA    currentVelo;   
 intA    currentDetune;
 doubleA currentPitch;   // actual pitch of the filter freq at this moment 
                         // (even when in slide)
 doubleA targetPitch;    // pitch, which we want to slide to
	doubleA pitchBend;      // position of the pitchwheel in semitones
	doubleA freqBend;       // factor for cutoff frequency determined by
                         // pitchbend and scaled by keytrack
 doubleA currentCutoff;  // actual frequency of the filter at this moment
                         // (even when in slide) without being scaled by the
                         // envelope
 doubleA targetCutoff;   // frequency, which we want to slide to
 doubleA currentReso;    // current resonance (determined by reso, resoByVel
                         // and Velocity) without being scaled by the cutoff
                         // frequency
 doubleA currentVel;


	doubleA prevOutput;     // previous output sample (is used for feedback 
                         // filter-FM)
};

//----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):
__forceinline double AggressorFilterSection::getSample(double In)
{
 static doubleA instCutoff, instReso;  //instantaneous cutoff frequency and resonance
 static doubleA out;                   //output signal

 //apply the slide to the cutoff frequency:
 if(currentCutoff!=targetCutoff)
 { 
  //currentPitch  += pitchIncPerSamp;  //not really needed (?)
  currentCutoff *= freqFacPerSamp;
  slideSampCount++;
 }
 //compensate for roundoff-error in frequency after slide:
 if(slideSampCount>=slideSamples)
 {
  currentCutoff = targetCutoff;
  //currentPitch  = targetPitch;  //not really needed (?)
 }

 //apply filter envelope:
 instCutoff = currentCutoff*filtEnv.getSample();

	//apply the pitchwheel:
	instCutoff *= freqBend;

	//modulate cutoff-freq by filters input and output signals:
 instCutoff = instCutoff + cutoffByInput*In*instCutoff;
 instCutoff = instCutoff + cutoffByOutput*prevOutput*instCutoff;

 //calculate the resonance at this instant from the values of currentReso (determined by
 //reso, resoByVel and velocity) and the instantaneous cutoff frequency:
 instReso = pow( (instCutoff*recRefCutoff), resoScaleExponent) * currentReso;

	//modulate resonance by filters input and output signals:
	instReso = instReso + In*resoByInput*instReso;
	instReso = instReso + prevOutput*resoByOutput*instReso;

	//set up the filters and calculate output sample...
 if(mode<=9)
 {
  biQuadCasc.setFreq(instCutoff);
  biQuadCasc.setQ(instReso);
  biQuadCasc.calcCoeffs();
  //biQuadCasc.convertDirectToLadder();
  //out = biQuadCasc.getSampleLadder1(In); 
  out = biQuadCasc.getSampleDirect1(In); 
  //out = biQuadCasc.getSampleDirect2(In); 
  //out = biQuadCasc.getSampleAutoChoose(In); 
 }
 else if(mode<=11)
 {
  moogFilter.setCutoff(instCutoff);
  moogFilter.setFeedback(instReso);
  moogFilter.calcCoeffs();
  out = moogFilter.getSample(In+dcOffset); //+ fatness*fatFilter.getSample(inSamp);
 }
 else
 {
  out = 0.0; //+ fatness*fatFilter.getSample(inSamp);
 }


 //...and store it in prevOutput (will be used in the next call for modulation):
 prevOutput = out;
 CLIP(prevOutput,1); //make sure that it is not larger than one

	//return the sample:
	return wetVolume*out + dryVolume*In;
}


#endif // !defined(AggressorFilterSection_h_Included)
