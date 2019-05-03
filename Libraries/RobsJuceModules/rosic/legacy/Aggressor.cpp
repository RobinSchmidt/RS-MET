#include "Aggressor.h"

//construction/destruction
Aggressor::Aggressor()
{
	currentNote = -1;
	currentVelo = 0;
 filterIsOn  = true;
 distIsOn    = true;
 //antiAliasFilter.setCutoff(16000.0);

 list<MidiNoteEvent> noteList; // for test only
}
Aggressor::~Aggressor()
{
}
//------------------------------------------------------------------------------------------------------------
//parameter settings:
//general:
void Aggressor::setPitchBend(double PitchBend)
{
	oscSection.setPitchBend(PitchBend);
	filterSection.setPitchBend(PitchBend); //filter adjudsts its frequency according to keytrack 
	                                    //not only via Note-Numbers but also via pitch-wheel
}
void Aggressor::setSampleRate(double SampleRate)
{
 oscSection.setSampleRate(2*SampleRate);
 antiAliasFilter.setSampleRate(2*SampleRate);
 filterSection.setSampleRate(SampleRate);
 ampEnv1.setSampleRate(SampleRate);
	ampEnv2.setSampleRate(SampleRate);
	distorter.setSampleRate(SampleRate);
}

void Aggressor::setVolume(double Volume)
{
	volume   = Volume;
 finalVol = DB2AMP(volume + (currentVelo-64)*volByVelAc/63);
}
void Aggressor::setVolByVel(double VolByVel)
{
 volByVel   = VolByVel;
 volByVelAc = accent*volByVel;
 finalVol   = DB2AMP(volume + (currentVelo-64)*volByVelAc/63);
}
void Aggressor::setMonoPoly(long MonoPoly)
{
}
void Aggressor::setSlideTime(double SlideTime)
{
 if(SlideTime>=0)
 {
  oscSection.setSlideTime(SlideTime);
  filterSection.setSlideTime(SlideTime);
 }
}
void Aggressor::setAccent(double Accent)
{
 accent = Accent;

	oscSection.setAccent(accent);

 //oscSection.setPitchEnvPeakByVel(accent*pitchEnvPeakByVel);
 filterSection.setCutoffByVel(accent*cutoffByVel);
 filterSection.setResoByVel(accent*resoByVel);
 filterSection.setFiltEnvPeakByVel(accent*filtEnvPeakByVel);
 filterSection.setFiltEnvDurByVel(accent*filtEnvDurByVel);
 filterSection.setDryWetByVel(accent*filtDryWetByVel);

 distDryWetByVelAc   = accent*distDryWetByVel;
 ampEnv1PeakByVelAc  = accent*ampEnv1PeakByVel;
 ampEnv2PeakByVelAc  = accent*ampEnv2PeakByVel;
 volByVelAc          = accent*volByVel;
 //ampEnv....
}
void Aggressor::setEnvSlope(double EnvSlope)
{
 oscSection.setEnvSlope(EnvSlope);
 filterSection.setEnvSlope(EnvSlope);
 ampEnv1.setTauScale(EnvSlope);
 ampEnv2.setTauScale(EnvSlope);
}
void Aggressor::setVibDepth(double VibDepth)
{oscSection.setVibDepth(VibDepth);}

void Aggressor::setVibSpeed(double VibSpeed)
{oscSection.setVibRate(VibSpeed);}

void Aggressor::setVibWave(long VibWave)
{oscSection.setVibWave(VibWave);}

void Aggressor::setPwmDepth(double PwmDepth)
{oscSection.setPwmDepth(PwmDepth);
}

void Aggressor::setPwmSpeed(double PwmSpeed)
{oscSection.setPwmRate(PwmSpeed);
}

void Aggressor::setPwmWave(long PwmWave)
{oscSection.setPwmWave(PwmWave);
}

//oscillator 1:
void Aggressor::setWaveForm1(long WaveForm1)
{oscSection.setWaveForm1(WaveForm1);
}
void Aggressor::setStartPhase1(double StartPhase1)
{oscSection.setStartPhase1(StartPhase1);
}
void Aggressor::setVol1Start(double Vol1Start)
{oscSection.setVol1Start(Vol1Start);
}
void Aggressor::setVol1StartByVel(double Vol1StartByVel)
{oscSection.setVol1StartByVel(Vol1StartByVel);
}
void Aggressor::setVol1Time(double Vol1Time)
{oscSection.setVol1Time(Vol1Time);
}
void Aggressor::setVol1End(double Vol1End)
{oscSection.setVol1End(Vol1End);
}
void Aggressor::setPitch1Start(double Pitch1Start)
{oscSection.setPitch1Start(Pitch1Start);
}
void Aggressor::setPitch1StartByVel(double Pitch1StartByVel)
{oscSection.setPitch1StartByVel(Pitch1StartByVel);
}
void Aggressor::setPitch1Time(double Pitch1Time)
{oscSection.setPitch1Time(Pitch1Time);
}
void Aggressor::setPitch1End(double Pitch1End)
{oscSection.setPitch1End(Pitch1End);
}
void Aggressor::setOsc1Lpf(double Osc1Lpf)
{oscSection.setOsc1Lpf(Osc1Lpf);
}
void Aggressor::setOsc1Hpf(double Osc1Hpf)
{oscSection.setOsc1Hpf(Osc1Hpf);
}
void Aggressor::setOsc1Apf(double Osc1Apf)
{oscSection.setOsc1Apf(Osc1Apf);
}
void Aggressor::setPulseWidth1(double PulseWidth1)
{oscSection.setPulseWidth1(PulseWidth1);
}
void Aggressor::setDensity(double Density)
{oscSection.setDensity(Density);
}
void Aggressor::setDetuneStart(double DetuneStart)
{oscSection.setDetuneStart(DetuneStart);
}
void Aggressor::setDetuneTime(double DetuneTime)
{oscSection.setDetuneTime(DetuneTime);
}
void Aggressor::setDetuneEnd(double DetuneEnd)
{oscSection.setDetuneEnd(DetuneEnd);
}
void Aggressor::setDetuneSpacing(int DetuneSpacing)
{oscSection.setDetuneSpacing(DetuneSpacing);
}
void Aggressor::setDetuneRatio(double DetuneRatio)
{oscSection.setDetuneRatio(DetuneRatio);
}
void Aggressor::setDetunePhaseSpread(double DetunePhaseSpread)
{oscSection.setDetunePhaseSpread(DetunePhaseSpread);
}

//modulations:
void Aggressor::setAmBy2Start(double AmBy2Start)
{oscSection.setAmBy2Start(AmBy2Start);
}
void Aggressor::setAmBy2StartByVel(double AmBy2StartByVel)
{oscSection.setAmBy2StartByVel(AmBy2StartByVel);
}
void Aggressor::setAmBy2Time(double AmBy2Time)
{oscSection.setAmBy2Time(AmBy2Time);
}
void Aggressor::setAmBy2End(double AmBy2End)
{oscSection.setAmBy2End(AmBy2End);
}
void Aggressor::setAmBy3Start(double AmBy3Start)
{oscSection.setAmBy3Start(AmBy3Start);
}
void Aggressor::setAmBy3StartByVel(double AmBy3StartByVel)
{oscSection.setAmBy3StartByVel(AmBy3StartByVel);
}
void Aggressor::setAmBy3Time(double AmBy3Time)
{oscSection.setAmBy3Time(AmBy3Time);
}
void Aggressor::setAmBy3End(double AmBy3End)
{oscSection.setAmBy3End(AmBy3End);
}
void Aggressor::setFmBy1Start(double FmBy1Start)
{oscSection.setFmBy1Start(FmBy1Start);
}
void Aggressor::setFmBy1StartByVel(double FmBy1StartByVel)
{oscSection.setFmBy1StartByVel(FmBy1StartByVel);
}
void Aggressor::setFmBy1Time(double FmBy1Time)
{oscSection.setFmBy1Time(FmBy1Time);
}
void Aggressor::setFmBy1End(double FmBy1End)
{oscSection.setFmBy1End(FmBy1End);
}
void Aggressor::setFmBy2Start(double FmBy2Start)
{oscSection.setFmBy2Start(FmBy2Start);
}
void Aggressor::setFmBy2StartByVel(double FmBy2StartByVel)
{oscSection.setFmBy2StartByVel(FmBy2StartByVel);
}
void Aggressor::setFmBy2Time(double FmBy2Time)
{oscSection.setFmBy2Time(FmBy2Time);
}
void Aggressor::setFmBy2End(double FmBy2End)
{oscSection.setFmBy2End(FmBy2End);
}
void Aggressor::setFmBy3Start(double FmBy3Start)
{oscSection.setFmBy3Start(FmBy3Start);
}
void Aggressor::setFmBy3StartByVel(double FmBy3StartByVel)
{oscSection.setFmBy3StartByVel(FmBy3StartByVel);
}
void Aggressor::setFmBy3Time(double FmBy3Time)
{oscSection.setFmBy3Time(FmBy3Time);
}
void Aggressor::setFmBy3End(double FmBy3End)
{oscSection.setFmBy3End(FmBy3End);
}
void Aggressor::setModulationMode(long ModulationMode)
{oscSection.setModulationMode(ModulationMode);
}
void Aggressor::setSyncMode(long SyncMode)
{oscSection.setSyncMode(SyncMode);
}

//oscillator 2:
void Aggressor::setWaveForm2(long WaveForm2)
{oscSection.setWaveForm2(WaveForm2);
}
void Aggressor::setStartPhase2(double StartPhase2)
{oscSection.setStartPhase2(StartPhase2);
}
void Aggressor::setVol2Start(double Vol2Start)
{oscSection.setVol2Start(Vol2Start);
}
void Aggressor::setVol2StartByVel(double Vol2StartByVel)
{oscSection.setVol2StartByVel(Vol2StartByVel);
}
void Aggressor::setVol2Time(double Vol2Time)
{oscSection.setVol2Time(Vol2Time);
}
void Aggressor::setVol2End(double Vol2End)
{oscSection.setVol2End(Vol2End);
}
void Aggressor::setPitch2Start(double Pitch2Start)
{oscSection.setPitch2Start(Pitch2Start);
}
void Aggressor::setPitch2StartByVel(double Pitch2StartByVel)
{oscSection.setPitch2StartByVel(Pitch2StartByVel);
}
void Aggressor::setPitch2Time(double Pitch2Time)
{oscSection.setPitch2Time(Pitch2Time);
}
void Aggressor::setPitch2End(double Pitch2End)
{oscSection.setPitch2End(Pitch2End);
}
void Aggressor::setOsc2Lpf(double Osc2Lpf)
{oscSection.setOsc2Lpf(Osc2Lpf);
}
void Aggressor::setOsc2Hpf(double Osc2Hpf)
{oscSection.setOsc2Hpf(Osc2Hpf);
}
void Aggressor::setOsc2Apf(double Osc2Apf)
{oscSection.setOsc2Apf(Osc2Apf);
}
void Aggressor::setPulseWidth2(double PulseWidth2)
{oscSection.setPulseWidth2(PulseWidth2);
}
void Aggressor::setTuneCoarse2(double TuneCoarse2)
{oscSection.setTuneCoarse2(TuneCoarse2);
}
void Aggressor::setTuneFine2(double TuneFine2)
{oscSection.setTuneFine2(TuneFine2);
}

//oscillator 3:
void Aggressor::setWaveForm3(long WaveForm3)
{oscSection.setWaveForm3(WaveForm3);
}
void Aggressor::setStartPhase3(double StartPhase3)
{oscSection.setStartPhase3(StartPhase3);
}
void Aggressor::setVol3Start(double Vol3Start)
{oscSection.setVol3Start(Vol3Start);
}
void Aggressor::setVol3StartByVel(double Vol3StartByVel)
{oscSection.setVol3StartByVel(Vol3StartByVel);
}
void Aggressor::setVol3Time(double Vol3Time)
{oscSection.setVol3Time(Vol3Time);
}
void Aggressor::setVol3End(double Vol3End)
{oscSection.setVol3End(Vol3End);
}
void Aggressor::setPitch3Start(double Pitch3Start)
{oscSection.setPitch3Start(Pitch3Start);
}
void Aggressor::setPitch3StartByVel(double Pitch3StartByVel)
{oscSection.setPitch3StartByVel(Pitch3StartByVel);
}
void Aggressor::setPitch3Time(double Pitch3Time)
{oscSection.setPitch3Time(Pitch3Time);
}
void Aggressor::setPitch3End(double Pitch3End)
{oscSection.setPitch3End(Pitch3End);
}
void Aggressor::setOsc3Lpf(double Osc3Lpf)
{oscSection.setOsc3Lpf(Osc3Lpf);
}
void Aggressor::setOsc3Hpf(double Osc3Hpf)
{oscSection.setOsc3Hpf(Osc3Hpf);
}
void Aggressor::setOsc3Apf(double Osc3Apf)
{oscSection.setOsc3Apf(Osc3Apf);
}
void Aggressor::setPulseWidth3(double PulseWidth3)
{oscSection.setPulseWidth3(PulseWidth3);
}
void Aggressor::setTuneCoarse3(double TuneCoarse3)
{oscSection.setTuneCoarse3(TuneCoarse3);
}
void Aggressor::setTuneFine3(double TuneFine3)
{oscSection.setTuneFine3(TuneFine3);
}

//ringmodulation output volumes:
void Aggressor::setRm12Start(double Rm12Start)
{oscSection.setRm12Start(Rm12Start);
}
void Aggressor::setRm12StartByVel(double Rm12StartByVel)
{oscSection.setRm12StartByVel(Rm12StartByVel);
}
void Aggressor::setRm12Time(double Rm12Time)
{oscSection.setRm12Time(Rm12Time);
}
void Aggressor::setRm12End(double Rm12End)
{oscSection.setRm12End(Rm12End);
}
void Aggressor::setRm13Start(double Rm13Start)
{oscSection.setRm13Start(Rm13Start);
}
void Aggressor::setRm13StartByVel(double Rm13StartByVel)
{oscSection.setRm13StartByVel(Rm13StartByVel);
}
void Aggressor::setRm13Time(double Rm13Time)
{oscSection.setRm13Time(Rm13Time);
}
void Aggressor::setRm13End(double Rm13End)
{oscSection.setRm13End(Rm13End);
}
void Aggressor::setRm23Start(double Rm23Start)
{oscSection.setRm23Start(Rm23Start);
}
void Aggressor::setRm23StartByVel(double Rm23StartByVel)
{oscSection.setRm23StartByVel(Rm23StartByVel);
}
void Aggressor::setRm23Time(double Rm23Time)
{oscSection.setRm23Time(Rm23Time);
}
void Aggressor::setRm23End(double Rm23End)
{oscSection.setRm23End(Rm23End);
}

//filter settings:
void Aggressor::setFiltMode(long Mode)
{
 if(Mode==0)
  filterIsOn = false;
 else
  filterIsOn = true;
 filterSection.setMode(Mode);
}

void Aggressor::setFiltStages(long Stages)
{filterSection.setStages(Stages);}

void Aggressor::setFiltDryWet(double FiltDryWet)
{
 filterSection.setDryWet(FiltDryWet);
}

void Aggressor::setFiltDryWetByVel(double FiltDryWetByVel)
{
 filtDryWetByVel = FiltDryWetByVel;
 filterSection.setDryWetByVel(accent*filtDryWetByVel);
}

void Aggressor::setFiltCutoff(double Cutoff)
{filterSection.setCutoff(Cutoff);}

void Aggressor::setFiltCutoffByVel(double CutoffByVel)
{
 cutoffByVel = CutoffByVel;
 filterSection.setCutoffByVel(accent*cutoffByVel);
}

void Aggressor::setFiltCutoffByInput(double CutoffByInput)
{
 filterSection.setCutoffByInput(CutoffByInput);
}

void Aggressor::setFiltCutoffByOutput(double CutoffByOutput)
{
 filterSection.setCutoffByOutput(CutoffByOutput);
}

void Aggressor::setFiltReso(double Reso)
{filterSection.setReso(Reso);}

void Aggressor::setFiltResoByVel(double ResoByVel)
{
 resoByVel = ResoByVel;
 filterSection.setResoByVel(accent*ResoByVel);
}

void Aggressor::setFiltResoByInput(double ResoByInput)
{filterSection.setResoByInput(ResoByInput);
}

void Aggressor::setFiltResoByOutput(double ResoByOutput)
{filterSection.setResoByOutput(ResoByOutput);
}

void Aggressor::setFiltGain(double Gain)
{filterSection.setGain(Gain);}

void Aggressor::setFiltResoByCutoff(double ResByCut)
{filterSection.setResoByCutoff(ResByCut);}

void Aggressor::setFiltDrive(double Drive)
{filterSection.setDrive(Drive);}

void Aggressor::setFiltDcOffset(double DcOffset)
{filterSection.setDcOffset(DcOffset);}

void Aggressor::setFiltRaiseToPower(double Power)
{//filterSection.setRaiseToPower(Power);
}

void Aggressor::setFiltFatness(double Fatness)
{//filterSection.setFatness(Fatness);
}

void Aggressor::setFiltKeyTrack(double KeyTrack)
{filterSection.setKeyTrack(KeyTrack);}

void Aggressor::setFiltRefKey(long RefKey)
{filterSection.setRefKey(RefKey);}

void Aggressor::setFiltOverSamp(long OverSamp)
{filterSection.setOverSamp(OverSamp);}


//filter envelope settings:
void Aggressor::setFiltEnvStart(double Start)
{filterSection.setFiltEnvStart(Start);}

void Aggressor::setFiltEnvAttack(double Attack)
{filterSection.setFiltEnvAttack(Attack);}

void Aggressor::setFiltEnvPeak(double Peak)
{filterSection.setFiltEnvPeak(Peak);}

void Aggressor::setFiltEnvPeakByVel(double PeakByVel)
{
 filtEnvPeakByVel = PeakByVel;
 filterSection.setFiltEnvPeakByVel(accent*PeakByVel);
}

void Aggressor::setFiltEnvHold(double Hold)
{
 filterSection.setFiltEnvHold(Hold);
}

void Aggressor::setFiltEnvDecay(double Decay)
{filterSection.setFiltEnvDecay(Decay);}

void Aggressor::setFiltEnvRelease(double Release)
{filterSection.setFiltEnvRelease(Release);}

void Aggressor::setFiltEnvEnd(double End)
{filterSection.setFiltEnvEnd(End);}

void Aggressor::setFiltEnvDurByVel(double EnvDurByVel)
{
 filtEnvDurByVel = EnvDurByVel;
 filterSection.setFiltEnvDurByVel(accent*EnvDurByVel);
}


//amplitude envelope 1 settings:
void Aggressor::setAmpEnv1Start(double Start)
{ampEnv1.setStart(Start);}
void Aggressor::setAmpEnv1Attack(double Attack)
{ampEnv1.setAttack(Attack);}
void Aggressor::setAmpEnv1Peak(double Peak)
{
 ampEnv1Peak = Peak;
 ampEnv1.setPeak(Peak);
}
void Aggressor::setAmpEnv1PeakByVel(double PeakByVel)
{
 ampEnv1PeakByVel   = PeakByVel;
 ampEnv1PeakByVelAc = accent*ampEnv1PeakByVel;
}
void Aggressor::setAmpEnv1Hold(double Hold)
{ampEnv1.setHold(Hold);
}

void Aggressor::setAmpEnv1Decay(double Decay)
{ampEnv1.setDecay(Decay);}
void Aggressor::setAmpEnv1Sustain(double Sustain)
{ampEnv1.setSustain(Sustain);}
void Aggressor::setAmpEnv1Release(double Release)
{ampEnv1.setRelease(Release);}

//distortion settings:
void Aggressor::setDistMode(long DistMode)
{
 if(DistMode==0)
  distIsOn = false;
 else
 {
  distIsOn = true;
  distorter.setMode(DistMode);
 }
}
void Aggressor::setDistDryWet(double DryWet)
{
 distDryWet = DryWet;
 distorter.setDryWet(DryWet);
}
void Aggressor::setDistDryWetByVel(double DryWetByVel)
{
 distDryWetByVel   = DryWetByVel;
 distDryWetByVelAc = accent*distDryWetByVel; 
}
void Aggressor::setDistDrive(double DistDrive)
{distorter.setDrive(DistDrive);}
void Aggressor::setDistDcOffset(double DcOffset)
{
 distorter.setDcOffset(DcOffset);
}
void Aggressor::setDistLowThresh(double LowThresh)
{distorter.setLowThresh(LowThresh);}
void Aggressor::setDistHighThresh(double HighThresh)
{distorter.setHighThresh(HighThresh);}
void Aggressor::setDistLowSat(double LowSat)
{distorter.setLowSat(LowSat);}
void Aggressor::setDistHighSat(double HighSat)
{distorter.setHighSat(HighSat);}
void Aggressor::setDistLowClamp(double LowClamp)
{distorter.setLowClamp(LowClamp);}
void Aggressor::setDistHighClamp(double HighClamp)
{distorter.setHighClamp(HighClamp);}
void Aggressor::setDistRectify(double Rectify)
{distorter.setRectify(Rectify);}

//amplitude envelope 2 settings:
void Aggressor::setAmpEnv2Start(double Start)
{ampEnv2.setStart(Start);}
void Aggressor::setAmpEnv2Attack(double Attack)
{ampEnv2.setAttack(Attack);}
void Aggressor::setAmpEnv2Peak(double Peak)
{
 ampEnv2Peak = Peak;
 ampEnv2.setPeak(Peak);
}
void Aggressor::setAmpEnv2PeakByVel(double PeakByVel)
{
 ampEnv2PeakByVel   = PeakByVel;
 ampEnv2PeakByVelAc = accent*ampEnv2PeakByVel;
}
void Aggressor::setAmpEnv2Hold(double Hold)
{ampEnv2.setHold(Hold);
}

void Aggressor::setAmpEnv2Decay(double Decay)
{ampEnv2.setDecay(Decay);}
void Aggressor::setAmpEnv2Sustain(double Sustain)
{ampEnv2.setSustain(Sustain);}
void Aggressor::setAmpEnv2Release(double Release)
{ampEnv2.setRelease(Release);}

//------------------------------------------------------------------------------------------------------------
//others:
void Aggressor::noteOff(long NoteNumber)
{
 // noteOffs are handled in noteOn() now, too
}

void Aggressor::noteOn(long NoteNumber, long Velocity, long Detune)
{
 if( Velocity == 0 )
 {
  MidiNoteEvent releasedNote(NoteNumber, 0);
  noteList.remove(releasedNote);

  // check if the note-list is empty now. if so, trigger a release. if not, 
  // slide to the note at the beginning of the list (this is the most recent
  // one which is still in the list):
  if( noteList.empty() )
  {
   oscSection.noteOff(currentNote);
   filterSection.noteOff(currentNote); 
   ampEnv1.noteOff();
   ampEnv2.noteOff();
   currentNote = -1;
	 	currentVelo = 0;
  }
  else
  {
   oscSection.slideToNote(noteList.front().getKey(), 
                          noteList.front().getVel(), 
                          noteList.front().getDetune());

   filterSection.slideToNote(noteList.front().getKey(), 
                             noteList.front().getVel(), 
                             noteList.front().getDetune());
  }

 } // end of if(Velocity==0)
 else // velocity was not zero
 {
  // check if the note-list is empty. if so, trigger a new note. if not, 
  // slide to the new note:
  if( noteList.empty() )
  {
			// retrigger osc and set filter buffers to zero only if amplitude is near 
   // zero - otherwise clicks will result:
   if(ampEnv1.endIsReached() || ampEnv2.endIsReached())
   {
    oscSection.resetOscillators();
    //filterSection.resetBuffers();
   }
   oscSection.triggerNote(NoteNumber, Velocity, Detune);
   filterSection.triggerNote(NoteNumber, Velocity, Detune);

   // scale the velocity dependent amplitude-envelope parameters:
   ampEnv1.setPeak(ampEnv1Peak + (Velocity-64)*ampEnv1PeakByVelAc/63);
   ampEnv2.setPeak(ampEnv2Peak + (Velocity-64)*ampEnv2PeakByVelAc/63);

   // trigger the amplitude envelopes:
   ampEnv1.trigger();
   ampEnv2.trigger();
  }
  else // noteList was not empty
  {
   oscSection.slideToNote(NoteNumber, Velocity, Detune);
   filterSection.slideToNote(NoteNumber, Velocity, Detune);
  }
  // scale the distortion amount by the note-velocity and adjust the final
  // volume factor:
  distorter.setDryWet(distDryWet + Velocity*distDryWetByVelAc/127);
  finalVol = DB2AMP(volume + (Velocity-64)*volByVelAc/63);

  // and we need to add the new note to our list, of course:
  MidiNoteEvent newNote(NoteNumber, Velocity);
  noteList.push_front(newNote);
 }

 currentNote   = NoteNumber;
	currentVelo   = Velocity;
 currentDetune = Detune;
}


