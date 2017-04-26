#include "AggressorOscSection.h"

//construction/destruction
AggressorOscSection::AggressorOscSection()
{
 //init member variables:
 slideTime      = 0.5;   //init slieTime to 0.5 seconds
 sampleRate     = 44100;
 slideSamples   = slideTime*sampleRate;
 slideSampCount = 0;
	syncMode       = 0;
 modulationMode = 0;
 currentNote    = -1;
 currentDetune  = 0;
 currentPitch   = 64.0;
 targetPitch    = 64.0;
	freqBend       = 1.0;
 currentFreq    = PITCH2FREQ(currentPitch);
 targetFreq     = PITCH2FREQ(targetPitch);

 osc1IsPlaying  = true;
 osc2IsPlaying  = true;
 osc3IsPlaying  = true;
 vibIsActive    = true;
 pwmIsActive    = true;
 modulationsOff = false;
 syncIsActive   = false;

 tuneCoarse2    = 0.0;
 tuneFine2      = 0.0;
 tuneFactor2    = 1.0;
 tuneCoarse3    = 0.0;
 tuneFine3      = 0.0;
 tuneFactor3    = 1.0;

 pulseWidth1    = 0.5;
 pulseWidth2    = 0.5;
 pulseWidth3    = 0.5;

 pitchIncPerSamp = 0.0;
 freqFacPerSamp  = 1.0;
 accent          = 1.0;
 vibDepth        = 0.0;
 pwmDepth        = 0.0;

 phase1             = 0.0;
 vol1Start          = 0.0;
 vol1StartByVel     = 0.0;
 vol1StartByVelAc   = 0.0;
 pitch1Start        = 0.0;
 pitch1StartByVel   = 0.0;
 pitch1StartByVelAc = 0.0;

 amBy2Start         = 0.0;
 amBy2StartByVel    = 0.0;
 amBy2StartByVelAc  = 0.0;

 amBy3Start         = 0.0;
 amBy3StartByVel    = 0.0;
 amBy3StartByVelAc  = 0.0;

 fmBy1Start         = 0.0;
 fmBy1StartByVel    = 0.0;
 fmBy1StartByVelAc  = 0.0;

 fmBy2Start         = 0.0;
 fmBy2StartByVel    = 0.0;
 fmBy2StartByVelAc  = 0.0;

 fmBy3Start         = 0.0;
 fmBy3StartByVel    = 0.0;
 fmBy3StartByVelAc  = 0.0;

 phase2             = 0.0;
 vol2Start          = 0.0;
 vol2StartByVel     = 0.0;
 vol2StartByVelAc   = 0.0;
 pitch2Start        = 0.0;
 pitch2StartByVel   = 0.0;
 pitch2StartByVelAc = 0.0;

 phase3             = 0.0;
 vol3Start          = 0.0;
 vol3StartByVel     = 0.0;
 vol3StartByVelAc   = 0.0;
 pitch3Start        = 0.0;
 pitch3StartByVel   = 0.0;
 pitch3StartByVelAc = 0.0;

 rm12Start          = 0.0;
 rm12StartByVel     = 0.0;
 rm12StartByVelAc   = 0.0;

 rm13Start          = 0.0;
 rm13StartByVel     = 0.0;
 rm13StartByVelAc   = 0.0;

 rm23Start          = 0.0;
 rm23StartByVel     = 0.0;
 rm23StartByVelAc   = 0.0;

 flt1.setLpfCutoff(20000.0);
 flt1.setHpfCutoff(20.0);
 flt1.setApfCutoff(20000.0);

 flt2.setLpfCutoff(20000.0);
 flt2.setHpfCutoff(20.0);
 flt2.setApfCutoff(20000.0);

 flt3.setLpfCutoff(20000.0);
 flt3.setHpfCutoff(20.0);
 flt3.setApfCutoff(20000.0);
}
AggressorOscSection::~AggressorOscSection()
{}
//------------------------------------------------------------------------------------------------------------
//parameter settings:
//general:
void AggressorOscSection::setSampleRate(double SampleRate)
{
 if(SampleRate>0)
  sampleRate = SampleRate;

	//tell the oscillators the new sample-rate:
 osc1.setSampleRate(sampleRate);
	osc2.setSampleRate(sampleRate);
	osc3.setSampleRate(sampleRate);

 //tell the lfo's the new sample-rate and let them 
 //re-calculate their increments:
 vibLfo.setSampleRate(sampleRate);
 pwmLfo.setSampleRate(sampleRate);
 vibLfo.calcIncrements();
 pwmLfo.calcIncrements();

	//tell the filters the new sample-rate:
 flt1.setSampleRate(sampleRate);
 flt2.setSampleRate(sampleRate);
 flt3.setSampleRate(sampleRate);

	//tell the ramp-envelopes the new sample-rate:
	ampRamp1.setSampleRate(sampleRate);
	ampRamp2.setSampleRate(sampleRate);
	ampRamp3.setSampleRate(sampleRate);
	freqRamp1.setSampleRate(sampleRate);
	freqRamp2.setSampleRate(sampleRate);
	freqRamp3.setSampleRate(sampleRate);
 rm12Ramp.setSampleRate(sampleRate);
 rm13Ramp.setSampleRate(sampleRate);
 rm23Ramp.setSampleRate(sampleRate);
	detuneRamp.setSampleRate(sampleRate);
	amBy2Ramp.setSampleRate(sampleRate);
	amBy3Ramp.setSampleRate(sampleRate);
	fmBy1Ramp.setSampleRate(sampleRate);
	fmBy2Ramp.setSampleRate(sampleRate);
	fmBy3Ramp.setSampleRate(sampleRate);
}

void AggressorOscSection::setMonoPoly(long MonoPoly)
{
}

void AggressorOscSection::setSlideTime(double SlideTime)
{
 if(SlideTime>=0)
  slideTime = SlideTime;
 slideSamples = slideTime*sampleRate;
}

void AggressorOscSection::setAccent(double Accent)
{
	accent             = Accent;
	vol1StartByVelAc   = accent*vol1StartByVel;
	vol2StartByVelAc   = accent*vol2StartByVel;
	vol3StartByVelAc   = accent*vol3StartByVel;
	pitch1StartByVelAc = accent*pitch1StartByVel;
	pitch2StartByVelAc = accent*pitch2StartByVel;
	pitch3StartByVelAc = accent*pitch3StartByVel;
	amBy2StartByVelAc  = accent*amBy2StartByVel;
	amBy3StartByVelAc  = accent*amBy3StartByVel;	
	fmBy1StartByVelAc  = accent*fmBy1StartByVel;
	fmBy2StartByVelAc  = accent*fmBy2StartByVel;
	fmBy3StartByVelAc  = accent*fmBy3StartByVel;
}

void AggressorOscSection::setEnvSlope(double EnvSlope)
{
 ampRamp1.setTauScale(EnvSlope);
 ampRamp2.setTauScale(EnvSlope);
 ampRamp3.setTauScale(EnvSlope);
 freqRamp1.setTauScale(EnvSlope);
 freqRamp2.setTauScale(EnvSlope);
 freqRamp3.setTauScale(EnvSlope);
 rm12Ramp.setTauScale(EnvSlope);
 rm13Ramp.setTauScale(EnvSlope);
 rm23Ramp.setTauScale(EnvSlope);
	detuneRamp.setTauScale(EnvSlope);
	amBy2Ramp.setTauScale(EnvSlope);
	amBy3Ramp.setTauScale(EnvSlope);
	fmBy1Ramp.setTauScale(EnvSlope);
	fmBy2Ramp.setTauScale(EnvSlope);
	fmBy3Ramp.setTauScale(EnvSlope);
}

void AggressorOscSection::setVibDepth(double VibDepth)
{vibDepth = 0.01*VibDepth;}
void AggressorOscSection::setVibRate(double VibRate)
{
 vibLfo.setFreq(VibRate);
 vibLfo.calcIncrements();
}
void AggressorOscSection::setVibWave(long VibWave)
{
 if(VibWave>0)
  vibIsActive = true;
 else
  vibIsActive = false;
 vibLfo.setWaveForm(VibWave);
}
void AggressorOscSection::setPwmDepth(double PwmDepth)
{pwmDepth = 0.01*PwmDepth;}
void AggressorOscSection::setPwmRate(double PwmRate)
{
 pwmLfo.setFreq(PwmRate);
 pwmLfo.calcIncrements();
}
void AggressorOscSection::setPwmWave(long PwmWave)
{
 if(PwmWave>0)
  pwmIsActive = true;
 else
  pwmIsActive = false;
 pwmLfo.setWaveForm(PwmWave);
}

//oscillator 1:
void AggressorOscSection::setWaveForm1(long WaveForm1)
{
 if(WaveForm1>0)
  osc1IsPlaying = true;
 else
  osc1IsPlaying = false;
 osc1.setWaveForm(WaveForm1);
}

void AggressorOscSection::setStartPhase1(double StartPhase1)
{
 phase1 = StartPhase1;
 osc1.setStartPhase(StartPhase1);
}

void AggressorOscSection::setVol1Start(double Vol1Start)
{
	vol1Start = Vol1Start;
	ampRamp1.setStart(DB2AMP(vol1Start));
}
void AggressorOscSection::setVol1StartByVel(double Vol1StartByVel)
{
	vol1StartByVel   = Vol1StartByVel;
	vol1StartByVelAc = accent*vol1StartByVel;
}
void AggressorOscSection::setVol1Time(double Vol1Time)
{ampRamp1.setTime(0.001*Vol1Time);}

void AggressorOscSection::setVol1End(double Vol1End)
{ampRamp1.setEnd(DB2AMP(Vol1End));}

void AggressorOscSection::setPitch1Start(double Pitch1Start)
{
	pitch1Start = Pitch1Start;
	freqRamp1.setStart(PITCHOFFSET2FREQFACTOR(pitch1Start));
}
void AggressorOscSection::setPitch1StartByVel(double Pitch1StartByVel)
{
	pitch1StartByVel   = Pitch1StartByVel;
	pitch1StartByVelAc = accent*pitch1StartByVel;
}
void AggressorOscSection::setPitch1Time(double Pitch1Time)
{freqRamp1.setTime(0.001*Pitch1Time);}

void AggressorOscSection::setPitch1End(double Pitch1End)
{freqRamp1.setEnd(PITCHOFFSET2FREQFACTOR(Pitch1End));}

void AggressorOscSection::setOsc1Lpf(double Osc1Lpf)
{
 flt1.setLpfCutoff(Osc1Lpf);
}

void AggressorOscSection::setOsc1Hpf(double Osc1Hpf)
{
 flt1.setHpfCutoff(Osc1Hpf);
}

void AggressorOscSection::setOsc1Apf(double Osc1Apf)
{
 flt1.setApfCutoff(Osc1Apf);
}

void AggressorOscSection::setPulseWidth1(double PulseWidth1)
{
 if( (PulseWidth1>0.0) && (PulseWidth1<100.0) )
	 pulseWidth1 = 0.01*PulseWidth1;
 osc1.setPulseWidth(pulseWidth1);
}

void AggressorOscSection::setDensity(double Density)
{
	osc1.setNumVoices((long) Density);
}

void AggressorOscSection::setDetuneStart(double DetuneStart)
{
	detuneRamp.setStart(DetuneStart);
}

void AggressorOscSection::setDetuneTime(double DetuneTime)
{
 detuneRamp.setTime(0.001*DetuneTime);	
}

void AggressorOscSection::setDetuneEnd(double DetuneEnd)
{
	detuneRamp.setEnd(DetuneEnd);
}
void AggressorOscSection::setDetuneRatio(double DetuneRatio)
{
 osc1.setDetuneRatio(DetuneRatio);
}
void AggressorOscSection::setDetunePhaseSpread(double DetunePhaseSpread)
{
 osc1.setDetunePhaseSpread(DetunePhaseSpread);
}



void AggressorOscSection::setDetuneSpacing(int DetuneSpacing)
{
	osc1.setFreqSpacing(DetuneSpacing);
}

//modulation:
void AggressorOscSection::setAmBy2Start(double AmBy2Start)
{
	amBy2Start = AmBy2Start;
 amBy2Ramp.setStart(amBy2Start);
}
void AggressorOscSection::setAmBy2StartByVel(double AmBy2StartByVel)
{
	amBy2StartByVel    = AmBy2StartByVel;
	amBy2StartByVelAc  = accent*amBy2StartByVel;
}
void AggressorOscSection::setAmBy2Time(double AmBy2Time)
{
 amBy2Ramp.setTime(0.001*AmBy2Time);
}
void AggressorOscSection::setAmBy2End(double AmBy2End)
{
 amBy2Ramp.setEnd(AmBy2End);
}

void AggressorOscSection::setAmBy3Start(double AmBy3Start)
{
	amBy3Start = AmBy3Start;
 amBy3Ramp.setStart(amBy3Start);
}
void AggressorOscSection::setAmBy3StartByVel(double AmBy3StartByVel)
{
	amBy3StartByVel    = AmBy3StartByVel;
	amBy3StartByVelAc  = accent*amBy3StartByVel;
}
void AggressorOscSection::setAmBy3Time(double AmBy3Time)
{
 amBy3Ramp.setTime(0.001*AmBy3Time);
}
void AggressorOscSection::setAmBy3End(double AmBy3End)
{
 amBy3Ramp.setEnd(AmBy3End);
}




void AggressorOscSection::setFmBy1Start(double FmBy1Start)
{
 if(FmBy1Start>=0 && FmBy1Start<=100)
	 fmBy1Start = 0.5*FmBy1Start;
 fmBy1Ramp.setStart(fmBy1Start);
}
void AggressorOscSection::setFmBy1StartByVel(double FmBy1StartByVel)
{
 if(FmBy1StartByVel>=0 && FmBy1StartByVel<=100)
	 fmBy1StartByVel    = 0.01*FmBy1StartByVel;
	fmBy1StartByVelAc  = accent*fmBy1StartByVel;
}
void AggressorOscSection::setFmBy1Time(double FmBy1Time)
{
 fmBy1Ramp.setTime(0.001*FmBy1Time);
}
void AggressorOscSection::setFmBy1End(double FmBy1End)
{
 if(FmBy1End>=0 && FmBy1End<=100)
  fmBy1Ramp.setEnd(0.5*FmBy1End);
}

void AggressorOscSection::setFmBy2Start(double FmBy2Start)
{
	fmBy2Start = FmBy2Start;
 fmBy2Ramp.setStart(fmBy2Start);
}
void AggressorOscSection::setFmBy2StartByVel(double FmBy2StartByVel)
{
	fmBy2StartByVel    = FmBy2StartByVel;
	fmBy2StartByVelAc  = accent*fmBy2StartByVel;
}
void AggressorOscSection::setFmBy2Time(double FmBy2Time)
{
 fmBy2Ramp.setTime(0.001*FmBy2Time);
}
void AggressorOscSection::setFmBy2End(double FmBy2End)
{
 fmBy2Ramp.setEnd(FmBy2End);
	//fm12Index = FmBy2End;
}

void AggressorOscSection::setFmBy3Start(double FmBy3Start)
{
	fmBy3Start = FmBy3Start;
 fmBy3Ramp.setStart(fmBy3Start);
}
void AggressorOscSection::setFmBy3StartByVel(double FmBy3StartByVel)
{
	fmBy3StartByVel    = FmBy3StartByVel;
	fmBy3StartByVelAc  = accent*fmBy3StartByVel;
}
void AggressorOscSection::setFmBy3Time(double FmBy3Time)
{
 fmBy3Ramp.setTime(0.001*FmBy3Time);
}
void AggressorOscSection::setFmBy3End(double FmBy3End)
{
 fmBy3Ramp.setEnd(FmBy3End);
}

void AggressorOscSection::setModulationMode(long ModulationMode)
{
 if(ModulationMode==0)
  modulationsOff = true;
 else
  modulationsOff = false;
	modulationMode = ModulationMode;
}

void AggressorOscSection::setSyncMode(long SyncMode)
{
 if(SyncMode==0)
  syncIsActive = false;
 else
  syncIsActive = true;
	syncMode = SyncMode;
}

//oscillator 2:
void AggressorOscSection::setWaveForm2(long WaveForm2)
{
 if(WaveForm2>0)
  osc2IsPlaying = true;
 else
  osc2IsPlaying = false;
 osc2.setWaveForm(WaveForm2);
}

void AggressorOscSection::setStartPhase2(double StartPhase2)
{
 phase2 = StartPhase2;
 osc2.setStartPhase(StartPhase2);
}

void AggressorOscSection::setVol2Start(double Vol2Start)
{
	vol2Start = Vol2Start;
	ampRamp2.setStart(DB2AMP(vol2Start));
}
void AggressorOscSection::setVol2StartByVel(double Vol2StartByVel)
{
	vol2StartByVel   = Vol2StartByVel;
	vol2StartByVelAc = accent*vol2StartByVel;
}
void AggressorOscSection::setVol2Time(double Vol2Time)
{ampRamp2.setTime(0.001*Vol2Time);}

void AggressorOscSection::setVol2End(double Vol2End)
{ampRamp2.setEnd(DB2AMP(Vol2End));}

void AggressorOscSection::setPitch2Start(double Pitch2Start)
{
	pitch2Start = Pitch2Start;
	freqRamp2.setStart(PITCHOFFSET2FREQFACTOR(pitch2Start));
}
void AggressorOscSection::setPitch2StartByVel(double Pitch2StartByVel)
{
	pitch2StartByVel   = Pitch2StartByVel;
	pitch2StartByVelAc = accent*pitch2StartByVel;
}
void AggressorOscSection::setPitch2Time(double Pitch2Time)
{freqRamp2.setTime(0.001*Pitch2Time);}

void AggressorOscSection::setPitch2End(double Pitch2End)
{freqRamp2.setEnd(PITCHOFFSET2FREQFACTOR(Pitch2End));}

void AggressorOscSection::setOsc2Lpf(double Osc2Lpf)
{
 flt2.setLpfCutoff(Osc2Lpf);
}

void AggressorOscSection::setOsc2Hpf(double Osc2Hpf)
{
 flt2.setHpfCutoff(Osc2Hpf);
}
void AggressorOscSection::setOsc2Apf(double Osc2Apf)
{
 flt2.setApfCutoff(Osc2Apf);
}
void AggressorOscSection::setPulseWidth2(double PulseWidth2)
{
 if( (PulseWidth2>0.0) && (PulseWidth2<100.0) )
	 pulseWidth2 = 0.01*PulseWidth2;
 osc2.setPulseWidth(pulseWidth2);
}

void AggressorOscSection::setTuneCoarse2(double TuneCoarse2)
{
 tuneCoarse2 = TuneCoarse2;
 tuneFactor2 = PITCHOFFSET2FREQFACTOR(tuneCoarse2 + 0.01*tuneFine2);
}

void AggressorOscSection::setTuneFine2(double TuneFine2)
{
 tuneFine2   = TuneFine2;
 tuneFactor2 = PITCHOFFSET2FREQFACTOR(tuneCoarse2 + 0.01*tuneFine2);
}

//oscillator 3:
void AggressorOscSection::setWaveForm3(long WaveForm3)
{
 if(WaveForm3>0)
  osc3IsPlaying = true;
 else
  osc3IsPlaying = false;
 osc3.setWaveForm(WaveForm3);
}

void AggressorOscSection::setStartPhase3(double StartPhase3)
{
 phase3 = StartPhase3;
 osc3.setStartPhase(StartPhase3);
}

void AggressorOscSection::setVol3Start(double Vol3Start)
{
	vol3Start = Vol3Start;
	ampRamp3.setStart(DB2AMP(vol3Start));
}
void AggressorOscSection::setVol3StartByVel(double Vol3StartByVel)
{
	vol3StartByVel   = Vol3StartByVel;
	vol3StartByVelAc = accent*vol3StartByVel;
}
void AggressorOscSection::setVol3Time(double Vol3Time)
{ampRamp3.setTime(0.001*Vol3Time);}

void AggressorOscSection::setVol3End(double Vol3End)
{ampRamp3.setEnd(DB2AMP(Vol3End));}

void AggressorOscSection::setPitch3Start(double Pitch3Start)
{
	pitch3Start = Pitch3Start;
	freqRamp3.setStart(PITCHOFFSET2FREQFACTOR(pitch3Start));
}
void AggressorOscSection::setPitch3StartByVel(double Pitch3StartByVel)
{
	pitch3StartByVel   = Pitch3StartByVel;
	pitch3StartByVelAc = accent*pitch3StartByVel;
}
void AggressorOscSection::setPitch3Time(double Pitch3Time)
{freqRamp3.setTime(0.001*Pitch3Time);}

void AggressorOscSection::setPitch3End(double Pitch3End)
{freqRamp3.setEnd(PITCHOFFSET2FREQFACTOR(Pitch3End));}

void AggressorOscSection::setOsc3Lpf(double Osc3Lpf)
{
 flt3.setLpfCutoff(Osc3Lpf);
}

void AggressorOscSection::setOsc3Hpf(double Osc3Hpf)
{
 flt3.setHpfCutoff(Osc3Hpf);
}
void AggressorOscSection::setOsc3Apf(double Osc3Apf)
{
 flt3.setApfCutoff(Osc3Apf);
}
void AggressorOscSection::setPulseWidth3(double PulseWidth3)
{
 if( (PulseWidth3>0.0) && (PulseWidth3<100.0) )
	 pulseWidth3 = 0.01*PulseWidth3;
 osc3.setPulseWidth(pulseWidth3);
}

void AggressorOscSection::setTuneCoarse3(double TuneCoarse3)
{
 tuneCoarse3 = TuneCoarse3;
 tuneFactor3 = PITCHOFFSET2FREQFACTOR(tuneCoarse3 + 0.01*tuneFine3);
}

void AggressorOscSection::setTuneFine3(double TuneFine3)
{
 tuneFine3   = TuneFine3;
 tuneFactor3 = PITCHOFFSET2FREQFACTOR(tuneCoarse3 + 0.01*tuneFine3);
}


//ringmodulation:
void AggressorOscSection::setRm12Start(double Rm12Start)
{
 rm12Start = Rm12Start;
 rm12Ramp.setStart(DB2AMP(rm12Start));
}
void AggressorOscSection::setRm12StartByVel(double Rm12StartByVel)
{
 rm12StartByVel   = Rm12StartByVel;
 rm12StartByVelAc = accent*rm12StartByVel;
}
void AggressorOscSection::setRm12Time(double Rm12Time)
{rm12Ramp.setTime(0.001*Rm12Time);}

void AggressorOscSection::setRm12End(double Rm12End)
{rm12Ramp.setEnd(DB2AMP(Rm12End));}

void AggressorOscSection::setRm13Start(double Rm13Start)
{
 rm13Start = Rm13Start;
 rm13Ramp.setStart(DB2AMP(rm13Start));
}
void AggressorOscSection::setRm13StartByVel(double Rm13StartByVel)
{
 rm13StartByVel   = Rm13StartByVel;
 rm13StartByVelAc = accent*rm13StartByVel;
}
void AggressorOscSection::setRm13Time(double Rm13Time)
{rm13Ramp.setTime(0.001*Rm13Time);}

void AggressorOscSection::setRm13End(double Rm13End)
{rm13Ramp.setEnd(DB2AMP(Rm13End));}

void AggressorOscSection::setRm23Start(double Rm23Start)
{
 rm23Start = Rm23Start;
 rm23Ramp.setStart(DB2AMP(rm23Start));
}
void AggressorOscSection::setRm23StartByVel(double Rm23StartByVel)
{
 rm23StartByVel   = Rm23StartByVel;
 rm23StartByVelAc = accent*rm23StartByVel;
}
void AggressorOscSection::setRm23Time(double Rm23Time)
{rm23Ramp.setTime(0.001*Rm23Time);}

void AggressorOscSection::setRm23End(double Rm23End)
{rm23Ramp.setEnd(DB2AMP(Rm23End));}


//-------------------------------------------------------------------------------------------------------
//event processing:

//(re)triggers pitchEnvelope and set oscillator to the
//new frequency:
void AggressorOscSection::triggerNote(long NoteNumber, long Velocity, long Detune)
{
	static double temp;

 currentNote  = NoteNumber;
 currentPitch = NoteNumber + 0.01 * (double) Detune;
 targetPitch  = currentPitch;
 currentFreq  = PITCH2FREQ(currentPitch);
 targetFreq   = currentFreq;

 //calculate the start values of the ramp-envelopes from the unscaled values, and the 
	//velocity dependence (which is in itselt dependent from accent):
	//calculate start value of osc1's amplitude ramp and set it:
	temp = DB2AMP(vol1Start + (Velocity-64)*vol1StartByVelAc/63);
	ampRamp1.setStart(temp);
	//calculate start value of osc1's amplitude ramp and set it:
	temp = PITCHOFFSET2FREQFACTOR(pitch1Start + (Velocity-64)*pitch1StartByVelAc/63);
	freqRamp1.setStart(temp);

	//the same for osc2 and osc3:
	temp = DB2AMP(vol2Start + (Velocity-64)*vol2StartByVelAc/63);
	ampRamp2.setStart(temp);
	temp = PITCHOFFSET2FREQFACTOR(pitch2Start + (Velocity-64)*pitch2StartByVelAc/63);
	freqRamp2.setStart(temp);

	temp = DB2AMP(vol3Start + (Velocity-64)*vol3StartByVelAc/63);
	ampRamp3.setStart(temp);
	temp = PITCHOFFSET2FREQFACTOR(pitch3Start + (Velocity-64)*pitch3StartByVelAc/63);
	freqRamp3.setStart(temp);

	//the same for the AM- and FM-ramp-envelopes:
 temp = amBy2Start + (Velocity-64)*amBy2StartByVelAc/63;
 amBy2Ramp.setStart(temp);

 temp = amBy3Start + (Velocity-64)*amBy3StartByVelAc/63;
 amBy3Ramp.setStart(temp);

 /* velocity->feedback-fm has to be implemented here
 temp = fmBy1Start + (Velocity-64)*fmBy1StartByVelAc/63;
 fmBy1Ramp.setStart(temp);
 */

 temp = fmBy2Start + (Velocity-64)*fmBy2StartByVelAc/63;
 fmBy2Ramp.setStart(temp);

 temp = fmBy3Start + (Velocity-64)*fmBy3StartByVelAc/63;
 fmBy3Ramp.setStart(temp);


	//the same for the ringmod-ramp-envelopes:
 temp = DB2AMP(rm12Start + (Velocity-64)*rm12StartByVelAc/63);
 rm12Ramp.setStart(temp);

 temp = DB2AMP(rm13Start + (Velocity-64)*rm13StartByVelAc/63);
 rm13Ramp.setStart(temp);

 temp = DB2AMP(rm23Start + (Velocity-64)*rm23StartByVelAc/63);
 rm23Ramp.setStart(temp);

	//trigger the ramp-envelops:
	ampRamp1.trigger();
	ampRamp2.trigger();
	ampRamp3.trigger();
	freqRamp1.trigger();
	freqRamp2.trigger();
	freqRamp3.trigger();
 rm12Ramp.trigger();
 rm13Ramp.trigger();
	rm23Ramp.trigger();
	detuneRamp.trigger();
	amBy2Ramp.trigger();
	amBy3Ramp.trigger();
	fmBy1Ramp.trigger();
	fmBy2Ramp.trigger();
	fmBy3Ramp.trigger();
}

//slide to the new pitch witchout retriggering envelope and oscillator:
void AggressorOscSection::slideToNote(long NoteNumber, long Velocity, long Detune)
{
 currentNote      = NoteNumber;
 targetPitch      = NoteNumber + 0.01 * (double) Detune;
 targetFreq       = PITCH2FREQ(targetPitch);

 double pitchDiff = targetPitch - currentPitch;

 if(slideSamples<1)  //immediate slide
 {
  pitchIncPerSamp  = 0;
  freqFacPerSamp   = 1;
  currentFreq      = targetFreq;
  currentPitch     = targetPitch;
  osc1.setFreq(targetFreq);
 }
 else                //scale frequency sample by sample (is done in getSample)
 {
  pitchIncPerSamp  = pitchDiff/slideSamples;
  freqFacPerSamp   = PITCHOFFSET2FREQFACTOR(pitchIncPerSamp);
 }
 
 //reset the sampleCounter:
 slideSampCount = 0;
}

void AggressorOscSection::noteOff(long NoteNumber)
{
 currentNote = -1;
}

void AggressorOscSection::setPitchBend(double PitchBend)
{freqBend = PITCHOFFSET2FREQFACTOR(PitchBend);}

//------------------------------------------------------------------------------------------------------------
//audio processing:
/*
__forceinline sample AggressorOscSection::getSample()
{
 static sample instFreq1, instFreq2, instFreq3; //frequency at this instant of time (for the 3 oscs)

	static sample ampFactor1, ampFactor2, ampFactor3; 
	//factor for the amplitudes (comes from envelope and amplitude modulation)

	static sample vibFactor;

	static sample instPulseWidth1, instPulseWidth2, instPulseWidth3, pulseWidthOffset;

	static sample freqDeviation; //for FM in osc1

	 //the audio signals (at different satges of the signal chain):
 static sample osc1Out,              osc2Out,              osc3Out,
               osc1WithFilter,       osc2WithFilter,       osc3WithFilter,
               osc1WithFilterAndEnv, osc2WithFilterAndEnv, osc3WithFilterAndEnv;

 static sample rm12Out, rm13Out, rm23Out;
 rm12Out = 0;
 rm13Out = 0;
 rm23Out = 0;

	static sample outputSample;

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
 //if(vibIsActive)
	 vibFactor = PITCHOFFSET2FREQFACTOR(vibLfo.getSample()*vibDepth);  //linear vibrato would be more efficient
 //else
  vibFactor = 1.0;

 //calculate an offset for the pulseWidths from the PWM-LFO:
 if(pwmIsActive)
	 pulseWidthOffset = pwmDepth*pwmLfo.getSample();  
 else
  pulseWidthOffset = 0.0;


 if(osc2IsPlaying)
 {
	 instFreq2  = currentFreq*freqRamp2.getSample();   //applies by the pitch ramp:
  instFreq2 *= tuneFactor2;                         //applies tuning-factor:
  instFreq2 *= vibFactor;                           //applies the vibrato
	 instFreq2 *= freqBend;                            //applies pitchbend

	 instPulseWidth2 = pulseWidth2 + pulseWidthOffset; //applies pwm to the static pulseWidth

	 ampFactor2 = ampRamp2.getSample();                //calculates the output amplitude

  //set up the oscillator:
  osc2.freq = instFreq2;
	 osc2.setPulseWidth(instPulseWidth2); //the function-call makes sure that pw is between 0.01 and 0.99
  osc2.calcIncrements();

  //calculate the oscillators output:
  //osc2Out = osc2.getSample();
  osc2Out = osc2.getSampleAntiAliased();

  //apply the static filter:
	 osc2WithFilter = bpf2.getSample(osc2Out);

  //apply the amp-ramp:
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
	 instFreq3  = currentFreq*freqRamp3.getSample();   //applies by the pitch ramp:
  instFreq3 *= tuneFactor3;                         //applies tuning-factor:
  instFreq3 *= vibFactor;                           //applies the vibrato
	 instFreq3 *= freqBend;                            //applies pitchbend

	 instPulseWidth3 = pulseWidth3 + pulseWidthOffset; //applies pwm to the static pulseWidth

	 ampFactor3 = ampRamp3.getSample();                //calculates the output amplitude

  //set up the oscillator:
  osc3.freq = instFreq3;
	 osc3.setPulseWidth(instPulseWidth3); //the function-call makes sure that pw is between 0.01 and 0.99
  osc3.calcIncrements();

  //calculate the oscillators output:
  //osc3Out = osc3.getSample();
  osc3Out = osc3.getSampleAntiAliased();

  //apply the static filter:
	 osc3WithFilter = bpf3.getSample(osc3Out);

  //apply the amp-ramp:
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
  instFreq1 *= vibFactor;                           //applies the vibrato
	 instFreq1 *= freqBend;                            //applies pitchbend

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
  osc1Out = osc1.getSampleAntiAliased();

  //apply the static filter:
	 osc1WithFilter = bpf1.getSample(osc1Out);

  //apply the amp-ramp:
  osc1WithFilterAndEnv = osc1WithFilter*ampFactor1;
 }
 else
 {
  osc1Out              = 0.0;
  osc1WithFilter       = 0.0;
  osc1WithFilterAndEnv = 0.0;
 }

	//sync if sync is active:
 if(syncIsActive)
 {
	 if( osc1.wraparoundOccurred && (syncMode==1 || syncMode==3) )
	 	osc2.resetPhase();
	 if( osc1.wraparoundOccurred && (syncMode==2 || syncMode==3) )
	 	osc3.resetPhase();
 }

 //add the signals:
	outputSample =   osc1WithFilterAndEnv + osc2WithFilterAndEnv + osc3WithFilterAndEnv 
		              + rm12Out + rm13Out + rm23Out;

 return outputSample;
}
*/




//------------------------------------------------------------------------------------------------------------
//others:
void AggressorOscSection::resetOscillators()
{
 if( phase1 <= 359.99 ) 
	 osc1.resetPhase();
 if( phase2 <= 359.99 ) 
	 osc2.resetPhase();
 if( phase3 <= 359.99 ) 
	 osc3.resetPhase();

 // lfo's are always re-trigeered:
	vibLfo.resetPhase();
	pwmLfo.resetPhase();
}


