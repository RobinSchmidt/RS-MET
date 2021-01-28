#ifndef rosic_AciDevil_h
#define rosic_AciDevil_h

//#include "../infrastructure/rosic_MidiNoteEvent.h"
//#include "../generators/rosic_BlendOscillator.h"
//#include "../filters/rosic_TeeBeeFilter.h"
//#include "../modulators/rosic_AnalogEnvelope.h"
//#include "../modulators/rosic_DecayEnvelope.h"
//#include "../filters/rosic_LeakyIntegrator.h"
//#include "../filters/rosic_EllipticQuarterBandFilter.h"
//#include "../sequencing/rosic_AcidSequencer.h"
//
//#include <list>
//using namespace std; // for the noteList

namespace rosic
{

  /**

  This is a monophonic bass-synth that aims to emulate the sound of the famous Roland TB 303 and
  goes a bit beyond.

  © Robin Schmidt (www.rs-met.com)

  */

  class AciDevil
  {

  public:

    //-----------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    AciDevil();

    /** Destructor. */
    ~AciDevil();

    //-----------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate (in Hz). */
    void setSampleRate(double newSampleRate);

    /** Sets the master volume level (in dB). */
    void setMasterLevel(double newLevel);     

    /** Sets the accent (in percent).  */
    void setAccent(double newAccent);

    /** Sets the slide-time (in ms). The TB-303 had a slide time of 60 ms. */
    void setSlideTime(double newSlideTime);

    /** Sets up the waveform continuously between saw and square - the input should be in the range 
    0...1 where 0 means pure saw and 1 means pure square. */
    void setWaveform(double newWaveform) { oscillator.setBlendFactor(newWaveform); }

    /** Sets the pulse-width for the pulse-oscilattor in percent. */
    void setPulseWidth(double newPulseWidth) { oscillator.setPulseWidth(newPulseWidth); }

    /** Sets up the suboscillator's waveform continuously between saw and square - the input should 
    be in the range 0...1 where 0 means pure saw and 1 means pure square. */
    void setSubOscWaveform(double newWaveform) { subOscillator.setBlendFactor(newWaveform); }

    /** Sets the volume level for the suboscillator (in decibels). */
    void setSubOscLevel(double newLevel) { subOscGain = RAPT::rsDbToAmp(newLevel); }

    /** Sets the filter's nominal cutoff frequency (in Hz). */
    void setCutoff(double newCutoff); 

    /** Sets the filter's resonance (in percent). */
    void setResonance(double newResonance) { filter.setResonance(newResonance); }

    /** Sets the filter mode as one of the values defined in TeeBeeFilter::modes. */
    void setFilterMode(int newMode) { filter.setMode(newMode); }

    /** Sets the modulation depth of the filter's cutoff frequency by the filter-envelope generator 
    (in percent of the nominal cutoff frequency). */
    void setEnvMod(double newEnvMod);

    /** Sets the filter envelope's attack time for non-accented notes (in milliseconds). 
    Devil Fish provides range of 0.3...30 ms for this parameter. */
    void setNormalAttack(double newNormalAttack) 
    { 
      normalAttack = newNormalAttack; 
      rc1.setTimeConstant(normalAttack);
    }

    /** Sets the filter envelope's attack time for accented notes (in milliseconds). In the 
    Devil Fish, accented notes have a fixed attack time of 3 ms.  */
    void setAccentAttack(double newAccentAttack) 
    { 
      accentAttack = newAccentAttack; 
      rc2.setTimeConstant(accentAttack);
    }

    /** Sets the filter envelope's decay time for non-accented notes (in milliseconds). 
    Devil Fish provides range of 30...3000 ms for this parameter. On the normal 303, this 
    parameter had a range of 200...2000 ms.  */
    void setNormalDecay(double newNormalDecay) { normalDecay = newNormalDecay; }

    /** Sets the filter envelope's decay time for accented notes (in milliseconds). 
    Devil Fish provides range of 30...3000 ms for this parameter. On the normal 303, this 
    parameter was fixed to 200 ms.  */
    void setAccentDecay(double newAccentDecay) { accentDecay = newAccentDecay; }

    /** Sets the amplitudes envelope's decay time (in milliseconds). Devil Fish provides range of 
    16...3000 ms for this parameter. On the normal 303, this parameter was fixed to 
    approximately 3-4 seconds.  */
    void setAmpDecay(double newAmpDecay) { ampEnv.setDecay(newAmpDecay); }

    /** Sets the amplitudes envelope's sustain level in decibels. Devil Fish uses the second half 
    of the range of the (amplitude) decay pot for this and lets the user adjust it between 0 
    and 100% of the full volume. In the normal 303, this parameter was fixed to zero. */
    void setAmpSustain(double newAmpSustain) { ampEnv.setSustainInDecibels(newAmpSustain); }

    /** Sets the amplitudes envelope's release time (in milliseconds). On the normal 303, this 
    parameter was fixed to .....  */
    void setAmpRelease(double newAmpRelease) 
    { 
      normalAmpRelease = newAmpRelease;
      ampEnv.setRelease(newAmpRelease); 
    }

    /** Sets the amount of drive for the clipper (in dB). */
    void setClipperDrive(double newDrive) { clipperGain = RAPT::rsDbToAmp(newDrive); }

    /** Sets the DC offset for the clipper. */
    void setClipperDC(double newDC) { clipperDC = newDC; }

    // todo: add more parameters: CliperSoftnessLow, ClipperSoftnessHigh - both from 0..1, default
    // is zero and means hardclip

    //-----------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the master volume level (in dB). */
    double getMasterLevel() const { return level; }

    /** Returns the accent (in percent). */
    double getAccent() const { return 100.0 * accent; }

    /** Returns the waveform as a continuous value between 0...1 where 0 means pure saw and 1 means 
    pure square. */
    double getWaveform() const { return oscillator.getBlendFactor(); }

    /** Returns the waveform of the suboscillator as a continuous value between 0...1 where 0 means 
    pure saw and 1 means pure square. */
    double getSubOscWaveform() const { return subOscillator.getBlendFactor(); }

    /** Returns the volume level for the suboscillator (in decibels). */
    double getSubOscLevel() const { return RAPT::rsAmpToDb(subOscGain); }

    /** Returns the filter's nominal cutoff frequency (in Hz). */
    double getCutoff() const { return cutoff; }

    /** Returns the slide-time (in ms). */
    double getSlideTime() const { return slideTime; }

    /** Returns the modulation depth of the filter's cutoff frequency by the filter-envelope 
    generator (in percent of the nominal cutoff frequency). */
    double getEnvMod() const { return envMod; }

    /** Returns the filter envelope's attack time for non-accented notes (in milliseconds). */
    double getNormalAttack() const { return normalAttack; }

    /** Returns the filter envelope's attack time for non-accented notes (in milliseconds). */
    double getAccentAttack() const { return accentAttack; }

    /** Returns the filter envelope's decay time for non-accented notes (in milliseconds). */
    double getNormalDecay() const { return normalDecay; }

    /** Returns the filter envelope's decay time for non-accented notes (in milliseconds). */
    double getAccentDecay() const { return accentDecay; }

    /** Returns the amplitudes envelope's decay time (in milliseconds). */
    double getAmpDecay() const { return ampEnv.getDecay(); }

    /** Returns the amplitudes envelope's sustain level (in dB). */
    double getAmpSustain() const { return RAPT::rsAmpToDb(ampEnv.getSustain()); }

    /** Returns the amplitudes envelope's release time (in milliseconds). */
    double getAmpRelease() const { return normalAmpRelease; }

    /** Returns the amount of drive for the clipper (in dB). */
    double getClipperDrive() const { return RAPT::rsAmpToDb(clipperGain); }

    //-----------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates onse output sample at a time. */
    INLINE double getSample(); 

    //-----------------------------------------------------------------------------------------------
    // event handling:

    /** Accepts note-on events (note offs are also handled here as note ons with velocity zero). */ 
    void noteOn(int noteNumber, int velocity, double detune);

    /** Turns all possibly running notes off. */
    void allNotesOff();

    /** Sets the pitchbend value. Expects a normalized value between -+1. The actual pitch-shift 
    will then be pitchWheelRange times the passed value. */ 
    void setPitchBend(double newPitchBend);  

    //-----------------------------------------------------------------------------------------------
    // embedded objects: 

    MipMappedWaveTable waveTable1, waveTable2, subWaveTable1, subWaveTable2;
    BlendOscillator    oscillator, subOscillator;
    TeeBeeFilter       filter;
    AnalogEnvelope     ampEnv; 
    DecayEnvelope      mainEnv;
    LeakyIntegrator    pitchSlewLimiter, ampDeClicker; // ampEnv3; //ampEnv2, 
    LeakyIntegrator    rc1, rc2;
    rsOnePoleFilterDD  hp1, hp2;                      // highpasses
    AcidSequencer      sequencer;
    rsEllipticQuarterBandFilter antiAliasFilter;


  protected:

    /** Triggers a note (called either directly in noteOn or in getSample when the sequencer is 
    used). */
    void triggerNote(int noteNumber, bool hasAccent);

    /** Slides to a note (called either directly in noteOn or in getSample when the sequencer is 
    used). */
    void slideToNote(int noteNumber, bool hasAccent);

    /** Releases a note (called either directly in noteOn or in getSample when the sequencer is 
    used). */
    void releaseNote(int noteNumber);

    /** Sets the decay-time of the main envelope and updates the normalizers n1, n2 accordingly. */
    void setMainEnvDecay(double newDecay);

    void calculateEnvModScalerAndOffset();

    /** Updates the normalizer n1 according to the time-constant of rc1 and the decay-time of the
    main envelope generator. */
    void updateNormalizer1();

    /** Updates the normalizer n2 according to the time-constant of rc2 and the decay-time of the
    main envelope generator. */
    void updateNormalizer2();

    static const int oversampling = 4;

    double ampScaler;        // final volume as raw factor
    double oscFreq;          // frequecy of the oscillator (without pitchbend)
    double sampleRate;       // the (non-oversampled) sample rate
    double level;            // master volume level (in dB)
    double levelByVel;       // velocity dependence of the level (in dB)
    double accent;           // scales all "byVel" parameters
    double slideTime;        // the time to slide from one note to another (in ms)
    double cutoff;           // nominal cutoff frequency of the filter
    double envMod;           // strength of the envelope modulation in semitones
    double envUpFraction;    // fraction of the envelope that goes upward
    double envOffset;        // offset for the normalized envelope ('bipolarity' parameter)
    double envScaler;        // scale-factor for the normalized envelope (derived from envMod)
    double normalAttack;     // attack time for the filter envelope on non-accented notes
    double accentAttack;     // attack time for the filter envelope on accented notes
    double normalDecay;      // decay time for the filter envelope on non-accented notes
    double accentDecay;      // decay time for the filter envelope on accented notes

    double normalAmpRelease;
    double accentAmpRelease;

    double subOscGain;       // gain factor for the suboscillator
    double accentGain;       // between 0.0...1.0 - to scale the 3rd amp-envelope on accents
    double clipperGain;      // gain factor for the input into the clipper
    double clipperDC;        // DC offset for the clipper
    double pitchWheelRange;  // range of pitch wheen in plus/minus semitones
    double pitchWheelFactor; // scale factor for oscillator frequency from pitch-wheel
    double n1, n2;           // normalizers

    int    currentNote;      // note which is currently played (-1 if none)
    int    currentVel;       // velocity of currently played note
    int    noteOffCountDown; // a countdown variable till next note-off in sequencer mode
    bool   slideToNextNote;  // indicate that we need to slide to the next note in sequencer mode
    bool   idle;             // flag to indicate that we have currently nothing to do in getSample

    std::list<MidiNoteEvent> noteList;

  };

  //-------------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed to 
  // be called at audio-rate (they can't be put into the .cpp file):

  INLINE double AciDevil::getSample()
  {
    //if( sequencer.getSequencerMode() == AcidSequencer::OFF && ampEnv.endIsReached() )
    //  return 0.0;
    if( idle )
      return 0.0;

    // check the sequencer if we have some note to trigger:
    if( sequencer.getSequencerMode() != AcidSequencer::OFF )
    {
      noteOffCountDown--;
      if( noteOffCountDown == 0 || sequencer.isRunning() == false )
        releaseNote(currentNote);

      AcidNote *note = sequencer.getNote();
      if( note != NULL )
      {
        if( note->gate == true && currentNote != -1)
        {
          int key = note->key + 12*note->octave + currentNote;
          key = RAPT::rsClip(key, 0, 127);

          if( !slideToNextNote )
            triggerNote(key, note->accent);
          else
            slideToNote(key, note->accent);

          AcidNote* nextNote = sequencer.getNextScheduledNote();
          if( note->slide && nextNote->gate == true )
          {
            noteOffCountDown = INT_MAX;
            slideToNextNote  = true;
          }
          else
          {
            noteOffCountDown = sequencer.getStepLengthInSamples();
            slideToNextNote  = false;
          }
        }
      }
    }

    // calculate instantaneous oscillator frequency and set up the oscillator:
    double instFreq = pitchSlewLimiter.getSample(oscFreq);
    oscillator.setFrequency(instFreq*pitchWheelFactor);
    oscillator.calculateIncrement();
    subOscillator.setIncrement(0.5*oscillator.getIncrement());

    // calculate instantaneous cutoff frequency from the nominal cutoff and all its modifiers and 
    // set up the filter:
    double mainEnvOut = mainEnv.getSample();
    double tmp1       = n1 * rc1.getSample(mainEnvOut);
    double tmp2       = 0.0;
    if( accentGain > 0.0 )
      tmp2 = mainEnvOut;
    tmp2 = n2 * rc2.getSample(tmp2);
    tmp1              = (envMod/12.0) * ( tmp1 - 0.35 ); // \todo: make offset adjustable
    tmp2              = accentGain*tmp2;
    double instCutoff = cutoff * pow(2.0, tmp1+tmp2);
    filter.setCutoff(instCutoff);
      // \todo: include pitch-tracking later

    double ampEnvOut = ampEnv.getSample();
    //ampEnvOut += 0.45*filterEnvOut + accentGain*6.8*filterEnvOut; 
    if( ampEnv.isNoteOn() )
      ampEnvOut += 0.45*mainEnvOut + accentGain*4.0*mainEnvOut; 
    ampEnvOut = ampDeClicker.getSample(ampEnvOut);
    idle      = sequencer.getSequencerMode() == AcidSequencer::OFF 
                 && ampEnv.endIsReached() && ampEnvOut < 0.000001;

    // oversampled calculations:
    double tmp;
    for(int i=1; i<=oversampling; i++)
    {
      tmp  = oscillator.getSample();                   // the raw oscillator signal 
      tmp += subOscGain*subOscillator.getSample();     // suboscillator signal added
      tmp  = hp1.getSample(tmp);                       // pre-filter highpass
      tmp  = filter.getSample(tmp);                    // now it's filtered
      tmp *= ampEnvOut;                                // amplified
      tmp  = RAPT::rsClip(clipperGain*tmp, -1.0, 1.0); // distorted
      tmp  = hp2.getSample(tmp);                       // may operate without oversampling....
      tmp  = antiAliasFilter.getSample(tmp);           // anti-aliasing filtered
    }

    return tmp * ampScaler;
  }

}

#endif 
