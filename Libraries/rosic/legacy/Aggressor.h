/*

Aggressor.h: interface for the Aggressor class.

© Braindoc 2002 (www.braindoc.de)

This is a monophonic bass-synth which can be played with a great dynamic because
many parameters can respond to the velocity of the note. This goes far beyond the
standard "accent" parameter usually found in bass-synths which controls the velocity
response of several parameters simultaneously. Here the velocity response of each
parameter, that may be controlled by velocity, can be set seperately. Moreover, the 
velocity not only turns "accent" on and off but determines the amount of accentuation.
Also this synth provides much more parameters than standard bass-synths - for example
quite sophisticated envelopes for cutoff, amplitude and pitch, a multimode filter and
different distortion modes.

Of course it is still possible to simulate the sound of the famous Roland TB-303. 

*/

#if !defined(Aggressor_h_Included)
#define Aggressor_h_Included

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "AudioModule.h"
#include "MidiNoteEvent.h"
#include "AggressorOscSection.h"
#include "AggressorFilterSection.h"
#include "AntiAliasFilter.h"
#include "EllipticHalfbandFilter.h"

#include "AmpEnvRc.h"
#include "Distortion.h"

#include <list>
using namespace std; // for the noteList

class Aggressor
{
public:
 //construction/destruction:
 Aggressor();
 virtual ~Aggressor();
 
 //parameter settings:
 void setSampleRate       (double SampleRate); //may cause error
 void noteOn              (long   NoteNumber, long Velocity, long Detune);
 void noteOff             (long   NoteNumber);
	void setPitchBend        (double PitchBend);  //expects PitchBend in semitones 

 //general settings:
 void setVolume           (double Volume);     //volume in dB
 void setVolByVel         (double VolByVel);
 void setMonoPoly         (long   MonoPoly);
 void setSlideTime        (double SlideTime);
 void setAccent           (double Accent);     //between 0 and 1
 void setEnvSlope         (double EnvSlope);    
 void setVibDepth         (double VibDepth);
 void setVibSpeed         (double VibSpeed);
 void setVibWave          (long   VibWave);
 void setPwmDepth         (double PwmDepth);
 void setPwmSpeed         (double PwmSpeed);
 void setPwmWave          (long   PwmWave);

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

 //filter settings:
 void setFiltMode          (long   Mode);          //filter mode (LPF, HPF, BPF, etc.)
 void setFiltStages        (long   Stages);        //number of filter stages
 void setFiltDryWet        (double FiltDryWet);
 void setFiltDryWetByVel   (double FiltDryWetByVel);
 void setFiltCutoff        (double Cutoff);        //basic cutoff frequency (will be modified
                                                   //by envelope, velocity an keytrack
 void setFiltCutoffByVel   (double CutoffByVel);   //intensity of cutoff modification by velocity )
 void setFiltCutoffByInput (double CutoffByInput); //filter-FM by filter input signal
 void setFiltCutoffByOutput(double CutoffByOutput);//filter-FM by filter output signal (feedback)
 void setFiltReso          (double Reso);          //resonance of the filter
 void setFiltResoByVel     (double ResoByVel);     //intensity of q modification by velocity (in "dB")
 void setFiltResoByInput   (double ResoByInput);   //reso-modulation by filter input signal
 void setFiltResoByOutput  (double ResoByOutput);  //reso-modulation by filter output signal (feedback)
 void setFiltGain          (double Gain);          //gain for peaking an shelving filters
 void setFiltResoByCutoff  (double ResByCut);      //cutoff can affect the resonance
 void setFiltDrive         (double Drive);         //drive the filter into distortion (in dB)
 void setFiltDcOffset      (double DcOffset);      //add dc to the filters input
 void setFiltRaiseToPower  (double Power);         //raise input signal to a power before filtering
	                                                  //and do the inverse operation after filtering
 void setFiltFatness       (double Fatness);       //fatness
 void setFiltKeyTrack      (double KeyTrack);      //keytracking for the cutoff frequency
 void setFiltRefKey        (long   RefKey);        //reference key for keytracking
 void setFiltOverSamp      (long   OverSamp);      //303-filter can operate on a higher double-rate

 //filter envelope settings:
 void setFiltEnvStart    (double Start);
 void setFiltEnvAttack   (double Attack);      //attack time of the envenlope in ms
 void setFiltEnvPeak     (double Peak);        //highest point of the envelope 
 void setFiltEnvPeakByVel(double PeakByVel);  
 void setFiltEnvHold     (double Hold);        //time to stay on the peak 
 void setFiltEnvDecay    (double Decay);       //decay time
 void setFiltEnvRelease  (double Release);     //release time (sustain level is determined
                                               //by cutoff and its modificators)
 void setFiltEnvEnd      (double End);
 void setFiltEnvDurByVel (double EnvDurByVel); //velocity control of the duration of the 
                                               //envelope speed (high vel: short env)

 //amplitude envelope 1 (pre-distortion) settings:
 void setAmpEnv1Start     (double Start);
 void setAmpEnv1Attack    (double Attack);
 void setAmpEnv1Peak      (double Peak);
 void setAmpEnv1PeakByVel (double PeakByVel);
 void setAmpEnv1Hold      (double Hold);
 void setAmpEnv1Decay     (double Decay);
 void setAmpEnv1Sustain   (double Sustain);
 void setAmpEnv1Release   (double Release);

 //distortion settings:
 void setDistMode        (long   DistMode);
 void setDistDryWet      (double DryWet);
 void setDistDryWetByVel (double DryWetByVel);
 void setDistDrive       (double Drive);
 void setDistDcOffset    (double DcOffset);
 void setDistLowThresh   (double LowThresh);
 void setDistHighThresh  (double HighThresh);
 void setDistLowSat      (double LowSat);
 void setDistHighSat     (double HighSat);
 void setDistLowClamp    (double LowClamp);
 void setDistHighClamp   (double HighClamp);
 void setDistRectify     (double Rectify);

	//amplitude envelope 2 (post-distortion) settings:
 void setAmpEnv2Start     (double Start);
 void setAmpEnv2Attack    (double Attack);
 void setAmpEnv2Peak      (double Peak);
 void setAmpEnv2PeakByVel (double PeakByVel);
 void setAmpEnv2Hold      (double Hold);
 void setAmpEnv2Decay     (double Decay);
 void setAmpEnv2Sustain   (double Sustain);
 void setAmpEnv2Release   (double Release);
																																						
 //audio processing:
 __forceinline double getSample(); //has to be define inside the .h file to get __forceinline to work

protected:
 doubleA sampleRate;

 //embedded opcode objects:
 AggressorOscSection    oscSection;
 //AntiAliasFilter antiAliasFilter;
 EllipticHalfbandFilter antiAliasFilter;
 AggressorFilterSection filterSection;
 AmpEnvRc               ampEnv1, ampEnv2;
 Distortion             distorter;

 //state variables:
 bool   filterIsOn;   //indicate if the filter is active
 bool   distIsOn;     //indicate if distortion is active
 long   currentNote;  //note which is currently played (-1 if none)
	long   currentVelo;  //velocity of currently played note
 long   currentDetune;//detuning of the current note
 doubleA volume;       //volume parameter (in dB)
 doubleA accent;       //scales all "byVel" parameters before forwarding them
                      //to the opcodes
 doubleA finalVol;     //final volume factor (not in dB  calculated from volume, 
                      //volByVel, currentVelo and accent)
 
 doubleA pitchEnvPeakByVel, cutoffByVel, resoByVel, filtDryWetByVel, filtEnvPeakByVel, filtEnvDurByVel, 
         ampEnv1PeakByVel, ampEnv2PeakByVel, distDryWetByVel, volByVel; 
        // values which have to be scaled by accent; need to be stored here to remember them when
        // setAccent is called

 doubleA ampEnv1Peak,ampEnv2Peak, distDryWet;

 doubleA ampEnv1PeakByVelAc, ampEnv2PeakByVelAc, volByVelAc,
         distDryWetByVelAc;
        //values scaled by accent
        //pitch and filter values are scaled and transmitted to MonoGenerator and -Filter
        //these objects take care of the velocity of incoming note events themselves

 list<MidiNoteEvent> noteList;

};//end of class


//-----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):
__forceinline double Aggressor::getSample()
{
 static doubleA outSamp;

 if(ampEnv1.outputIsZero || ampEnv2.outputIsZero)
  return 0.0;

 // generate the source-signal (2x oversampled):
 outSamp  = oscSection.getSample();
 outSamp  = antiAliasFilter.getSampleDirect1(outSamp);
 outSamp  = oscSection.getSample();
 outSamp  = antiAliasFilter.getSampleDirect1(outSamp);

 //decimate: use only every other sample generated by monoGenerator - nothing
 //to be done here

 //filter the signal:
 if(filterIsOn)
  outSamp  = filterSection.getSample(outSamp); 

 //apply amplitude envelope:
 outSamp *= ampEnv1.getSample();

 //distort the signal:
 if(distIsOn)
  outSamp  = distorter.getSampleAntiAliased(outSamp);

 //apply 2nd amplitude envelope:
 outSamp *= ampEnv2.getSample();

 return finalVol*outSamp;
}



#endif // !defined(Aggressor_h_Included)
