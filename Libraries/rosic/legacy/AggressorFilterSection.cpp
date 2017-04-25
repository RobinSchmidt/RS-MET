#include "AggressorFilterSection.h"

//construction/destruction
AggressorFilterSection::AggressorFilterSection()
{
	sampleRate        = 44100.0;
 slideTime         = 0.5;   // init slideTime to 0.5 seconds
 slideSamples      = slideTime*sampleRate;
 slideSampCount    = 0;
 pitchIncPerSamp   = 0.0;
 freqFacPerSamp    = 1.0;

 mode              = 1;

 dryWet            = 1.0;
 dryWetByVel       = 0.0;
 dryVolume         = 0.0;
 wetVolume         = 1.0;

 cutoff            = 1000.0;
 cutoffByVel       = 0.0;
 cutoffByInput     = 0.0;
 cutoffByOutput    = 0.0;

 reso              = 1.0;
 resoByVel         = 0.0;
 resoByInput       = 0.0;
 resoByOutput      = 0.0;
 resoByCutoff      = 0.0;
 resoScaleExponent = 1.0;

	drive             = 0.0;
 dcOffset          = 0.0;

 filtEnvPeak       = 0.0; // in semitones
 filtEnvPeakByVel  = 0.0;
 filtEnvDurByVel   = 1.0;

 keyTrack          = 0.0;
 refKey            = 33;
 //currentPitch      = refKey;

 currentNote       = -1;
 currentVelo       =  0;
 currentDetune     =  0;
 currentPitch      =  0.0;
 targetPitch       =  0.0;
 pitchBend         =  0.0;
 freqBend          =  1.0;
 currentCutoff     = 1000.0;
 targetCutoff      = 1000.0;
 currentReso       = 1.0;

 currentVel        = 0.0; // ambiguity with int currentVelo - bäh!

 slideSampCount    = 0;
 prevOutput        = 0.0;

	fatFilter.setMode(1);
	fatFilter.setFeedback(0.0);
}

AggressorFilterSection::~AggressorFilterSection()
{
}
//------------------------------------------------------------------------------------------------------------
//parameter settings:
void AggressorFilterSection::setPitchBend(double PitchBend)
{
 pitchBend = PitchBend;
	freqBend  = PITCHOFFSET2FREQFACTOR(keyTrack * pitchBend);
}
void AggressorFilterSection::setSampleRate(double SampleRate)
{
	if(SampleRate>0)
	{
		sampleRate = SampleRate;
  biQuadCasc.setSampleRate(SampleRate);
	 moogFilter.setSampleRate(SampleRate);
	 fatFilter.setSampleRate(SampleRate);
  filtEnv.setSampleRate(SampleRate);
	}
}
void AggressorFilterSection::setSlideTime(double SlideTime)
{
 if(SlideTime>=0)
  slideTime = SlideTime;
 slideSamples = slideTime*sampleRate;
}
void AggressorFilterSection::setEnvSlope(double EnvSlope)
{
 filtEnv.setTauScale(EnvSlope);
}
void AggressorFilterSection::setMode(long Mode)
{
 mode = Mode;
 if(mode<=9)
  biQuadCasc.setMode(mode);  //use the biquad cascade
 else
  moogFilter.setMode(mode-9);  //use the moog filter
}
void AggressorFilterSection::setStages(long Stages)
{
 biQuadCasc.setNumStages(Stages);
 moogFilter.setOutputTap(Stages);
}

void AggressorFilterSection::setDryWet(double DryWet)
{
 dryWet = 0.01*DryWet;

 wetVolume = dryWet; 
 dryVolume = 1 - dryWet;
}

void AggressorFilterSection::setDryWetByVel(double DryWetByVel)
{
 dryWetByVel = 0.01*DryWetByVel;
}

void AggressorFilterSection::setCutoff(double Cutoff)
{
 cutoff        = Cutoff;

	//calculate the current cutoff frequency at it is determined by the value of cutoff,
 //keyTrack,NoteNumber,refKey,cutoffByVel and velocity:
 //scale cutoff by NoteNumber (according to keyTrack and refKey):
 currentCutoff = PITCHOFFSET2FREQFACTOR((currentPitch-refKey)*keyTrack)*cutoff;

 //scale cutoff frequency by velocity (unscaled @ velocity = 64)
 //maximum: +- cutoffByVelo semitones (@ velo = 127 and velo = 1)
 currentCutoff = PITCHOFFSET2FREQFACTOR((currentVelo-64)*cutoffByVel/63)*currentCutoff;

 //if cutoff is changed during a slide process, we have to update the slide parameters:
 if(slideSampCount<slideSamples)
 {
	 //calculate the target cutoff frequency at it is determined by the value of cutoff,
  //keyTrack,NoteNumber,refKey,cutoffByVel and velocity:
  //scale cutoff by NoteNumber (according to keyTrack and refKey):
  targetCutoff = PITCHOFFSET2FREQFACTOR((targetPitch-refKey)*keyTrack)*cutoff;

  //scale cutoff frequency by velocity (unscaled @ velocity = 64)
  //maximum: +- cutoffByVelo semitones (@ velo = 127 and velo = 1)
  targetCutoff = PITCHOFFSET2FREQFACTOR((currentVelo-64)*cutoffByVel/63)*targetCutoff;
 }
 else
  targetCutoff = currentCutoff;
}

void AggressorFilterSection::setCutoffByVel(double CutoffByVel)
{
	cutoffByVel = CutoffByVel;

 //see setCutoff for comments on the following lines:
 currentCutoff = PITCHOFFSET2FREQFACTOR((currentPitch-refKey)*keyTrack)*cutoff;
 currentCutoff = PITCHOFFSET2FREQFACTOR((currentVelo-64)*cutoffByVel/63)*currentCutoff;

 if(slideSampCount<slideSamples)
 {
  targetCutoff = PITCHOFFSET2FREQFACTOR((targetPitch-refKey)*keyTrack)*cutoff;
  targetCutoff = PITCHOFFSET2FREQFACTOR((currentVelo-64)*cutoffByVel/63)*targetCutoff;
 }
 else
  targetCutoff = currentCutoff;
}

void AggressorFilterSection::setCutoffByInput(double CutoffByInput)
{
 if( CutoffByInput>= -100 && CutoffByInput<=100)
  cutoffByInput = 0.01*CutoffByInput;
}

void AggressorFilterSection::setCutoffByOutput(double CutoffByOutput)
{
 if( CutoffByOutput>= -100 && CutoffByOutput<=100)
  cutoffByOutput = 0.01*CutoffByOutput;
}

void AggressorFilterSection::setReso(double Reso)
{
	reso        = Reso;
	currentReso = DB2AMP((currentVelo-64)*resoByVel/63) * reso;
}
void AggressorFilterSection::setResoByVel(double ResoByVel)
{
	resoByVel   = ResoByVel;
 currentReso = DB2AMP((currentVelo-64)*resoByVel/63) * reso;
}

void AggressorFilterSection::setResoByInput(double ResoByInput)
{
 if(ResoByInput>= -100 && ResoByInput<=100)
	 resoByInput = 0.01*ResoByInput;
}
void AggressorFilterSection::setResoByOutput(double ResoByOutput)
{
 if(ResoByOutput>= -100 && ResoByOutput<=100)
	 resoByOutput = 0.01*ResoByOutput;
}

void AggressorFilterSection::setGain(double Gain)
{biQuadCasc.setGain(Gain);}
void AggressorFilterSection::setResoByCutoff(double ResoByCut)
{
 resoByCutoff      = ResoByCut;
 resoScaleExponent = resoByCutoff/(20*log10(2.0));
}
void AggressorFilterSection::setKeyTrack(double KeyTrack)
{
	keyTrack  = KeyTrack;
	freqBend  = PITCHOFFSET2FREQFACTOR(keyTrack * pitchBend);

	//see setCutoff for comments on the following lines:
 currentCutoff = PITCHOFFSET2FREQFACTOR((currentPitch-refKey)*keyTrack)*cutoff;
 currentCutoff = PITCHOFFSET2FREQFACTOR((currentVelo-64)*cutoffByVel/63)*currentCutoff;

 if(slideSampCount<slideSamples)
 {
  targetCutoff = PITCHOFFSET2FREQFACTOR((targetPitch-refKey)*keyTrack)*cutoff;
  targetCutoff = PITCHOFFSET2FREQFACTOR((currentVelo-64)*cutoffByVel/63)*targetCutoff;
 }
 else
  targetCutoff = currentCutoff;
}
void AggressorFilterSection::setDrive(double Drive)
{
	drive = Drive;
 moogFilter.setDrive(drive);
}
void AggressorFilterSection::setDcOffset(double DcOffset)
{
 dcOffset = DcOffset;
}
void AggressorFilterSection::setRefKey(long RefKey)
{
 if( RefKey>=0 && RefKey<=127)
  refKey = RefKey; 

	//see setCutoff for comments on the following lines:
 currentCutoff = PITCHOFFSET2FREQFACTOR((currentPitch-refKey)*keyTrack)*cutoff;
 currentCutoff = PITCHOFFSET2FREQFACTOR((currentVelo-64)*cutoffByVel/63)*currentCutoff;

 if(slideSampCount<slideSamples)
 {
  targetCutoff = PITCHOFFSET2FREQFACTOR((targetPitch-refKey)*keyTrack)*cutoff;
  targetCutoff = PITCHOFFSET2FREQFACTOR((currentVelo-64)*cutoffByVel/63)*targetCutoff;
 }
 else
  targetCutoff = currentCutoff;
}
void AggressorFilterSection::setOverSamp(long OverSamp)
{
 moogFilter.setOverSamp(OverSamp);
}
//envelope settings:
void AggressorFilterSection::setFiltEnvStart(double Start)
{
 filtEnv.setStart(Start);
}
void AggressorFilterSection::setFiltEnvAttack(double Attack)
{filtEnv.setAttack(Attack);}
void AggressorFilterSection::setFiltEnvPeak(double Peak)
{
 filtEnvPeak = Peak;   //buffer the peak in this class to modify by velocity
 filtEnv.setPeak(Peak);
}
void AggressorFilterSection::setFiltEnvPeakByVel(double PeakByVel)
{filtEnvPeakByVel = PeakByVel;}
void AggressorFilterSection::setFiltEnvHold(double Hold)
{
 filtEnv.setHold(Hold);
}
void AggressorFilterSection::setFiltEnvDecay(double Decay)
{filtEnv.setDecay(Decay);}
void AggressorFilterSection::setFiltEnvRelease(double Release)
{filtEnv.setRelease(Release);}
void AggressorFilterSection::setFiltEnvEnd(double End)
{filtEnv.setEnd(End);}
void AggressorFilterSection::setFiltEnvDurByVel(double EnvDurByVel)
{filtEnvDurByVel = EnvDurByVel;}

//------------------------------------------------------------------------------------------------------------
//others:
void AggressorFilterSection::triggerNote(long NoteNumber, long Velocity, long Detune)
{
 currentNote    = NoteNumber;
 currentVelo    = Velocity;
 currentPitch   = NoteNumber + 0.01 * (double) Detune;
 targetPitch    = currentPitch;
 
 //calculate the current cutoff frequency at it is determined by the value of cutoff,
 //keyTrack,NoteNumber,refKey,cutoffByVel and velocity:
 //scale cutoff by NoteNumber (according to keyTrack and refKey):
 currentCutoff = PITCHOFFSET2FREQFACTOR((currentPitch-refKey)*keyTrack)*cutoff;

 //scale cutoff frequency by velocity (unscaled @ velocity = 64)
 //maximum: +- cutoffByVelo semitones (@ velo = 127 and velo = 1)
 currentCutoff *= PITCHOFFSET2FREQFACTOR((currentVelo-64)*cutoffByVel/63);

 //a new note is triggered, so no slide has to take place:
 targetCutoff   = currentCutoff;

 //calculate and set the resonance as it is determined by reso, resoByVel and Velocity:
 //biQuadCasc.setQ       ( DB2AMP((currentVelo-64)*resoByVel/63) * reso );
 //moogFilter.setFeedback( DB2AMP((currentVelo-64)*resoByVel/63) * reso );
 currentReso = DB2AMP((currentVelo-64)*resoByVel/63) * reso;

 //scale the filter envelopes peak by note velocity:
 filtEnv.setPeak( (currentVelo-64)*filtEnvPeakByVel/63 + filtEnvPeak );

 //scale the filter envelopes duration by note velocity:
 filtEnv.setTimeScale( pow(2, ((64-currentVelo) * filtEnvDurByVel/63) ) );

 //change the dry/wet ratio according to the velocity
 wetVolume = dryWet + currentVelo*dryWetByVel/127;  //may result in values > 1 
                                                    //(dryWet=1, dryWetByVel=1, velo=127)
 if(wetVolume>1)
  wetVolume = 1;
 dryVolume = 1 - wetVolume;

 //test:
 //wetVolume = 1;
 //dryVolume = 1-wetVolume;
 //wetVolume = dryWet + currentVelo*dryWetByVel/127;
 //dryVolume = 0;

 filtEnv.trigger();
}

//
void AggressorFilterSection::slideToNote(long NoteNumber, long Velocity, long Detune)
{
 currentNote      = NoteNumber;
 currentVelo      = Velocity;
 targetPitch      = NoteNumber + 0.01 * (double) Detune;

 //calculation of target cutoff seems to be buggy:

 //calculate the target cutoff frequency at it is determined by the value of cutoff,
 //keyTrack,NoteNumber,refKey,cutoffByVel and velocity:
 //scale cutoff by NoteNumber (according to keyTrack and refKey):
 targetCutoff = PITCHOFFSET2FREQFACTOR((targetPitch-refKey)*keyTrack)*cutoff;

 //scale cutoff frequency by velocity (unscaled @ velocity = 64)
 //maximum: +- cutoffByVelo semitones (@ velo = 127 and velo = 1)
 targetCutoff = PITCHOFFSET2FREQFACTOR((currentVelo-64)*cutoffByVel/63)*targetCutoff;

 //calculate an set the resonance as it is determined by reso, resoByVel and Velocity:
 //biQuadCasc.setQ       ( DB2AMP((currentVelo-64)*resoByVel/63) * reso );
 //moogFilter.setFeedback( DB2AMP((currentVelo-64)*resoByVel/63) * reso );
 currentReso = DB2AMP((currentVelo-64)*resoByVel/63) * reso;
 
	//calculate the difference between the targe-cutoff and the current cutoff in terms of a
	//pitch difference between the two frequnecies:
 double pitchDiff = FREQFACTOR2PITCHOFFSET(targetCutoff) - FREQFACTOR2PITCHOFFSET(currentCutoff);

 if(slideSamples<1)  //immediate slide
 {
  pitchIncPerSamp  = 0;
  freqFacPerSamp   = 1;
  currentCutoff    = targetCutoff;
  currentPitch     = targetPitch;
  biQuadCasc.setFreq(targetCutoff); //-> unnecessary?
 }
 else                //scale frequency sample by sample (is done in getSample)
 {
  pitchIncPerSamp  = pitchDiff/slideSamples;
  freqFacPerSamp   = PITCHOFFSET2FREQFACTOR(pitchIncPerSamp);
 }
 
 //reset the sampleCounter:
 slideSampCount = 0;
}
void AggressorFilterSection::noteOff(long NoteNumber)
{
 filtEnv.noteOff();
 currentNote = -1;
}
void AggressorFilterSection::resetBuffers()
{
	moogFilter.resetBuffers();
	fatFilter.resetBuffers();
	biQuadCasc.resetBuffers();

	prevOutput = 0.0;
}




