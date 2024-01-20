//#include "rosic_Quadrifex.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

Quadrifex::Quadrifex() : mixMatrix(5, 5)
{
  mutex.lock();

  slotRouting = R_1TO2TO3TO4;
  sampleRate  = 44100.0;
  bpm         = 120.0;
  wetLevel    = 0.0;
  setDryWet(1.0);
  for(int i=0; i<numEffectSlots; i++)
  {
    effectAlgorithmIndices[i] = BYPASS;
    effectModules[i]          = new BypassModule();
      // don't change this - the constructor of QuadrifexAudioModule relies on them to be of type
      // BypassModule in the beginning (using static_cast)
  }
  reset();

  mutex.unlock();
}

Quadrifex::~Quadrifex()
{
  mutex.lock();
  for(int i=0; i<numEffectSlots; i++)
    delete effectModules[i];
  mutex.unlock();
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void Quadrifex::setSampleRate(double newSampleRate)
{
  if( newSampleRate <= 0.0 )
  {
    DEBUG_BREAK;
    return;
  }
  mutex.lock();
  sampleRate = newSampleRate;
  for(int i=0; i<numEffectSlots; i++)
    effectModules[i]->setSampleRate(sampleRate);
  mutex.unlock();
}

void Quadrifex::setTempoInBPM(double newTempoInBPM)
{
  if( newTempoInBPM != bpm )
  {
    bpm = newTempoInBPM;
    mutex.lock();
    for(int i=0; i<numEffectSlots; i++)
      effectModules[i]->setTempoInBPM(bpm);
    mutex.unlock();
  }
}

void Quadrifex::setSlotRouting(int newRoutingIndex)
{
  if( newRoutingIndex < 0 || newRoutingIndex >= NUM_SLOT_ROUTINGS )
    return;
  slotRouting = newRoutingIndex;
}

void Quadrifex::setEffectAlgorithm(int slotIndex, int newAlgorithmIndex)
{
  if( slotIndex < 0 || slotIndex >= numEffectSlots )
    return;

  mutex.lock();

  // delete the old EffectModule:
  delete effectModules[slotIndex];

  // create and set up the new EffectModule for the slot:
  effectAlgorithmIndices[slotIndex] = newAlgorithmIndex;
  Module *newModule = NULL;
  switch( effectAlgorithmIndices[slotIndex] )
  {
  case MUTE:                 newModule = new MuteModule;                      break;
  case BYPASS:               newModule = new BypassModule;                    break;
  case BIT_CRUSHER:          newModule = new BitCrusherModule;                break;
  case CHORUS:               newModule = new ChorusModule;                    break;
  case COMB_BANK:            newModule = new CombBankModule;                  break;
  case COMB_RESONATOR:       newModule = new CombResonatorStereoModule;       break;
  case COMB_STEREOIZER:      newModule = new CombStereoizerModule;            break;
  case COMPRESSOR:           newModule = new SoftKneeCompressorModule;        break;
  case DUAL_TWO_POLE_FILTER: newModule = new DualTwoPoleFilterModule;         break;
  case EQUALIZER:            newModule = new EqualizerModule;                 break;
  case EXPANDER:             newModule = new SoftKneeExpanderModule;          break;
  case FLANGER:              newModule = new FlangerModule;                   break;
  case FORMANT_SHIFTER:      newModule = new FormantShifterModule(8192);      break;
  case FOUR_POLE_FILTER:     newModule = new FourPoleFilterModule;            break;
  case FREQUENCY_SHIFTER:    newModule = new FrequencyShifterStereoModule;    break;
  case HARMONICS:            newModule = new HarmonicsModule;                 break;
  case LADDER_FILTER:        newModule = new LadderFilterModule;              break;
  case LIMITER:              newModule = new LimiterModule;                   break;
  case MODULATED_ALLPASS:    newModule = new ModulatedAllpassModule;          break;
  case NOISE_GATE:           newModule = new NoiseGateModule;                 break;
  case NOISIFIER:            newModule = new NoisifierModule;                 break;
  case PHASER:               newModule = new PhaserModule;                    break;
  case PHASE_STEREOIZER:     newModule = new PhaseStereoizerModule;           break;
  case PINGPONG_ECHO:        newModule = new PingPongEchoModule;              break;
  case PITCH_SHIFTER:        newModule = new PitchShifterModule;              break;
  case REVERB:               newModule = new rsReverbModule;                  break;
  case RINGMODULATOR:        newModule = new RingModulatorModule;             break;
  case SIMPLE_DELAY:         newModule = new SimpleDelayModule;               break;
  case SINE_OSCILLATOR:      newModule = new SineOscillatorModule;            break;
  case SLOPE_FILTER:         newModule = new SlopeFilterModule;               break;
  case SSB_MODULATOR:        newModule = new SingleSidebandModulatorModule;   break;
  case SLEWRATE_LIMITER:     newModule = new SlewRateLimiterStereoModule;     break;
  case STEREO_PAN:           newModule = new StereoPanModule;                 break;
  case STEREO_WIDTH:         newModule = new StereoWidthModule;               break;
  case TREMOLO:              newModule = new TremoloModule;                   break;
  case TWO_POLE_FILTER:      newModule = new TwoPoleFilterModule;             break;
  case VIBRATO:              newModule = new VibratoModule;                   break;
  case WAH_WAH:              newModule = new WahWahModule;                    break;
  case WAVESHAPER:           newModule = new WaveShaperModule;                break;

  default:
    {
      newModule                         = new BypassModule;
      effectAlgorithmIndices[slotIndex] = BYPASS;
    }
  }
  effectModules[slotIndex] = newModule;
  effectModules[slotIndex]->setSampleRate(sampleRate);
  effectModules[slotIndex]->setTempoInBPM(bpm);
  effectModules[slotIndex]->reset();

  mutex.unlock();
}

void Quadrifex::setDryWet(double newDryWet)
{
  dryWet = newDryWet;
  RAPT::rsEqualPowerGainFactors(dryWet, &dryFactor, &wetFactor, 0.0, 1.0);
  wetFactor *= RAPT::rsDbToAmp(wetLevel);
}

void Quadrifex::setWetLevel(double newLevel)
{
  wetLevel = newLevel;
  RAPT::rsEqualPowerGainFactors(dryWet, &dryFactor, &wetFactor, 0.0, 1.0);
  wetFactor *= RAPT::rsDbToAmp(wetLevel);
}

//-------------------------------------------------------------------------------------------------
// inquiry:

int Quadrifex::getEffectAlgorithmIndex(int slotIndex) const
{
  if( slotIndex < 0 || slotIndex >= numEffectSlots )
    return MUTE;
  else
    return effectAlgorithmIndices[slotIndex];
}

//-------------------------------------------------------------------------------------------------
// others:

void Quadrifex::trigger()
{
  mutex.lock();
  for(int i=0; i<numEffectSlots; i++)
    effectModules[i]->trigger();
  mutex.unlock();
}

void Quadrifex::reset()
{
  for(int i=0; i<5; i++)
  {
    inputsL[i]  = 0.0;
    inputsR[i]  = 0.0;
    outputsL[i] = 0.0;
    outputsR[i] = 0.0;
  }
  mutex.lock();
  for(int i=0; i<numEffectSlots; i++)
    effectModules[i]->reset();
  mutex.unlock();
}

//-------------------------------------------------------------------------------------------------
// audio processing:

void Quadrifex::processBlock(double* inOutL, double* inOutR, int numFrames)
{
  mutex.lock();

  double tmp1L = 0.0;
  double tmp1R = 0.0;
  double tmp2L = 0.0;
  double tmp2R = 0.0;
  double tmp3L = 0.0;
  double tmp3R = 0.0;
  double tmp4L = 0.0;
  double tmp4R = 0.0;
  double inL   = inOutL[0];
  double inR   = inOutR[0];

//  // some debug stuff:
//  double mixMatrixDebug[5][5];
//  for(int i=0; i<5; i++)
//  {
//    for(int j=0; j<5; j++)
//      mixMatrixDebug[i][j] = mixMatrix.getMatrixEntry(i, j);
//  }

  int n;
  switch( slotRouting )
  {
  case R_1TO2TO3TO4:
    {
      for(n=0; n<numFrames; n++)
      {
        inL = tmp1L = inOutL[n];
        inR = tmp1R = inOutR[n];
        effectModules[0]->processSampleFrame(&tmp1L, &tmp1R);
        effectModules[1]->processSampleFrame(&tmp1L, &tmp1R);
        effectModules[2]->processSampleFrame(&tmp1L, &tmp1R);
        effectModules[3]->processSampleFrame(&tmp1L, &tmp1R);
        inOutL[n] = (dryFactor*inL + wetFactor*tmp1L);
        inOutR[n] = (dryFactor*inR + wetFactor*tmp1R);
      }
    }
    break;
  case R_1TO2TO3_PLUS4:
    {
      for(n=0; n<numFrames; n++)
      {
        inL = tmp1L = tmp4L = inOutL[n];
        inR = tmp1R = tmp4R = inOutR[n];
        effectModules[0]->processSampleFrame(&tmp1L, &tmp1R);
        effectModules[1]->processSampleFrame(&tmp1L, &tmp1R);
        effectModules[2]->processSampleFrame(&tmp1L, &tmp1R);
        effectModules[3]->processSampleFrame(&tmp4L, &tmp4R);
        tmp1L += tmp4L;
        tmp1R += tmp4R;
        inOutL[n] = (dryFactor*inL + wetFactor*tmp1L);
        inOutR[n] = (dryFactor*inR + wetFactor*tmp1R);
      }
    }
    break;
  case R_1TO2_PLUS_3TO4:
    {
      for(n=0; n<numFrames; n++)
      {
        inL = tmp1L = tmp3L = inOutL[n];
        inR = tmp1R = tmp3R = inOutR[n];
        effectModules[0]->processSampleFrame(&tmp1L, &tmp1R);
        effectModules[1]->processSampleFrame(&tmp1L, &tmp1R);
        effectModules[2]->processSampleFrame(&tmp3L, &tmp3R);
        effectModules[3]->processSampleFrame(&tmp3L, &tmp3R);
        tmp1L += tmp3L;
        tmp1R += tmp3R;
        inOutL[n] = (dryFactor*inL + wetFactor*tmp1L);
        inOutR[n] = (dryFactor*inR + wetFactor*tmp1R);
      }
    }
    break;
  case R_1PLUS2PLUS3PLUS4:
    {
      for(n=0; n<numFrames; n++)
      {
        inL = tmp1L = tmp2L = tmp3L = tmp4L = inOutL[n];
        inR = tmp1R = tmp2R = tmp3R = tmp4R = inOutR[n];
        effectModules[0]->processSampleFrame(&tmp1L, &tmp1R);
        effectModules[1]->processSampleFrame(&tmp2L, &tmp2R);
        effectModules[2]->processSampleFrame(&tmp3L, &tmp3R);
        effectModules[3]->processSampleFrame(&tmp4L, &tmp4R);
        tmp1L += tmp2L+tmp3L+tmp4L;
        tmp1R += tmp2R+tmp3R+tmp4R;
        inOutL[n] = (dryFactor*inL + wetFactor*tmp1L);
        inOutR[n] = (dryFactor*inR + wetFactor*tmp1R);
      }
    }
    break;
  case R_1PLUS2PLUS3_TO_4:
    {
      for(n=0; n<numFrames; n++)
      {
        inL = tmp1L = tmp2L = tmp3L = inOutL[n];
        inR = tmp1R = tmp2R = tmp3R = inOutR[n];
        effectModules[0]->processSampleFrame(&tmp1L, &tmp1R);
        effectModules[1]->processSampleFrame(&tmp2L, &tmp2R);
        effectModules[2]->processSampleFrame(&tmp3L, &tmp3R);
        tmp1L += tmp2L+tmp3L;
        tmp1R += tmp2R+tmp3R;
        effectModules[3]->processSampleFrame(&tmp1L, &tmp1R);
        inOutL[n] = (dryFactor*inL + wetFactor*tmp1L);
        inOutR[n] = (dryFactor*inR + wetFactor*tmp1R);
      }
    }
    break;
  case R_1_TO_2PLUS3_TO_4:
    {
      for(n=0; n<numFrames; n++)
      {
        inL = tmp1L = inOutL[n];
        inR = tmp1R = inOutR[n];
        effectModules[0]->processSampleFrame(&tmp1L, &tmp1R);
        tmp2L = tmp3L = tmp1L;
        tmp2R = tmp3R = tmp1R;
        effectModules[1]->processSampleFrame(&tmp2L, &tmp2R);
        effectModules[2]->processSampleFrame(&tmp3L, &tmp3R);
        tmp1L += tmp2L+tmp3L;
        tmp1R += tmp2R+tmp3R;
        effectModules[3]->processSampleFrame(&tmp1L, &tmp1R);
        inOutL[n] = (dryFactor*inL + wetFactor*tmp1L);
        inOutR[n] = (dryFactor*inR + wetFactor*tmp1R);
      }
    }
    break;
  case R_1PLUS2_TO_3TO4:
    {
      for(n=0; n<numFrames; n++)
      {
        inL = tmp1L = tmp2L = inOutL[n];
        inR = tmp1R = tmp2R = inOutR[n];
        effectModules[0]->processSampleFrame(&tmp1L, &tmp1R);
        effectModules[1]->processSampleFrame(&tmp2L, &tmp2R);
        tmp1L += tmp2L;
        tmp1R += tmp2R;
        effectModules[2]->processSampleFrame(&tmp1L, &tmp1R);
        effectModules[3]->processSampleFrame(&tmp1L, &tmp1R);
        inOutL[n] = (dryFactor*inL + wetFactor*tmp1L);
        inOutR[n] = (dryFactor*inR + wetFactor*tmp1R);
      }

    }
    break;
  case R_1TO2_TO_3PLUS4:
    {
      for(n=0; n<numFrames; n++)
      {
        inL = tmp1L = inOutL[n];
        inR = tmp1R = inOutR[n];
        effectModules[0]->processSampleFrame(&tmp1L, &tmp1R);
        effectModules[1]->processSampleFrame(&tmp1L, &tmp1R);
        tmp3L = tmp1L;
        tmp3R = tmp1R;
        effectModules[2]->processSampleFrame(&tmp1L, &tmp1R);
        effectModules[3]->processSampleFrame(&tmp3L, &tmp3R);
        tmp1L += tmp3L;
        tmp1R += tmp3R;
        inOutL[n] = (dryFactor*inL + wetFactor*tmp1L);
        inOutR[n] = (dryFactor*inR + wetFactor*tmp1R);
      }
    }
    break;
  case R_1PLUS2_TO_3PLUS4:
    {
      for(n=0; n<numFrames; n++)
      {
        inL = tmp1L = tmp2L = inOutL[n];
        inR = tmp1R = tmp2R = inOutR[n];
        effectModules[0]->processSampleFrame(&tmp1L, &tmp1R);
        effectModules[1]->processSampleFrame(&tmp2L, &tmp2R);
        tmp1L += tmp2L;
        tmp1R += tmp2R;
        tmp2L  = tmp1L;
        tmp2R  = tmp1R;
        effectModules[2]->processSampleFrame(&tmp1L, &tmp1R);
        effectModules[3]->processSampleFrame(&tmp2L, &tmp2R);
        tmp1L += tmp2L;
        tmp1R += tmp2R;
        inOutL[n] = (dryFactor*inL + wetFactor*tmp1L);
        inOutR[n] = (dryFactor*inR + wetFactor*tmp1R);
      }
    }
    break;
  case R_1_TO_2PLUS3PLUS4:
    {
      for(n=0; n<numFrames; n++)
      {
        inL = tmp1L = inOutL[n];
        inR = tmp1R = inOutR[n];
        effectModules[0]->processSampleFrame(&tmp1L, &tmp1R);
        tmp2L = tmp3L = tmp4L = tmp1L;
        tmp2R = tmp3R = tmp4R = tmp1L;
        effectModules[1]->processSampleFrame(&tmp2L, &tmp2R);
        effectModules[2]->processSampleFrame(&tmp3L, &tmp3R);
        effectModules[3]->processSampleFrame(&tmp4L, &tmp4R);
        tmp1L = tmp2L+tmp3L+tmp4L;
        tmp1R = tmp2R+tmp3R+tmp4L;
        inOutL[n] = (dryFactor*inL + wetFactor*tmp1L);
        inOutR[n] = (dryFactor*inR + wetFactor*tmp1R);
      }
    }
    break;

  case MATRIX:
    {
      for(n=0; n<numFrames; n++)
      {
        inL   = inOutL[n];
        inR   = inOutR[n];

        // establish input for 1st effect unit:
        tmp1L =   mixMatrix.getMatrixEntryFast(0, 0) * inL
                + mixMatrix.getMatrixEntryFast(1, 0) * outputsL[0]
                + mixMatrix.getMatrixEntryFast(2, 0) * outputsL[1]
                + mixMatrix.getMatrixEntryFast(3, 0) * outputsL[2]
                + mixMatrix.getMatrixEntryFast(4, 0) * outputsL[3];
        tmp1R =   mixMatrix.getMatrixEntryFast(0, 0) * inR
                + mixMatrix.getMatrixEntryFast(1, 0) * outputsR[0]
                + mixMatrix.getMatrixEntryFast(2, 0) * outputsR[1]
                + mixMatrix.getMatrixEntryFast(3, 0) * outputsR[2]
                + mixMatrix.getMatrixEntryFast(4, 0) * outputsR[3];

        // apply the 1st effect unit:
        effectModules[0]->processSampleFrame(&tmp1L, &tmp1R);

        // remember the ouutput of the 1st effect unit for the next call:
        outputsL[0] = tmp1L;
        outputsR[0] = tmp1R;

        // ...and so on:
        tmp1L =   mixMatrix.getMatrixEntryFast(0, 1) * inL
                + mixMatrix.getMatrixEntryFast(1, 1) * outputsL[0]
                + mixMatrix.getMatrixEntryFast(2, 1) * outputsL[1]
                + mixMatrix.getMatrixEntryFast(3, 1) * outputsL[2]
                + mixMatrix.getMatrixEntryFast(4, 1) * outputsL[3];
        tmp1R =   mixMatrix.getMatrixEntryFast(0, 1) * inR
                + mixMatrix.getMatrixEntryFast(1, 1) * outputsR[0]
                + mixMatrix.getMatrixEntryFast(2, 1) * outputsR[1]
                + mixMatrix.getMatrixEntryFast(3, 1) * outputsR[2]
                + mixMatrix.getMatrixEntryFast(4, 1) * outputsR[3];
        effectModules[1]->processSampleFrame(&tmp1L, &tmp1R);
        outputsL[1] = tmp1L;
        outputsR[1] = tmp1R;

        tmp1L =   mixMatrix.getMatrixEntryFast(0, 2) * inL
                + mixMatrix.getMatrixEntryFast(1, 2) * outputsL[0]
                + mixMatrix.getMatrixEntryFast(2, 2) * outputsL[1]
                + mixMatrix.getMatrixEntryFast(3, 2) * outputsL[2]
                + mixMatrix.getMatrixEntryFast(4, 2) * outputsL[3];
        tmp1R =   mixMatrix.getMatrixEntryFast(0, 2) * inR
                + mixMatrix.getMatrixEntryFast(1, 2) * outputsR[0]
                + mixMatrix.getMatrixEntryFast(2, 2) * outputsR[1]
                + mixMatrix.getMatrixEntryFast(3, 2) * outputsR[2]
                + mixMatrix.getMatrixEntryFast(4, 2) * outputsR[3];
        effectModules[2]->processSampleFrame(&tmp1L, &tmp1R);
        outputsL[2] = tmp1L;
        outputsR[2] = tmp1R;

        tmp1L =   mixMatrix.getMatrixEntryFast(0, 3) * inL
                + mixMatrix.getMatrixEntryFast(1, 3) * outputsL[0]
                + mixMatrix.getMatrixEntryFast(2, 3) * outputsL[1]
                + mixMatrix.getMatrixEntryFast(3, 3) * outputsL[2]
                + mixMatrix.getMatrixEntryFast(4, 3) * outputsL[3];
        tmp1R =   mixMatrix.getMatrixEntryFast(0, 3) * inR
                + mixMatrix.getMatrixEntryFast(1, 3) * outputsR[0]
                + mixMatrix.getMatrixEntryFast(2, 3) * outputsR[1]
                + mixMatrix.getMatrixEntryFast(3, 3) * outputsR[2]
                + mixMatrix.getMatrixEntryFast(4, 3) * outputsR[3];
        effectModules[3]->processSampleFrame(&tmp1L, &tmp1R);
        outputsL[3] = tmp1L;
        outputsR[3] = tmp1R;



        tmp1L =   mixMatrix.getMatrixEntryFast(0, 4) * inL
                + mixMatrix.getMatrixEntryFast(1, 4) * outputsL[0]
                + mixMatrix.getMatrixEntryFast(2, 4) * outputsL[1]
                + mixMatrix.getMatrixEntryFast(3, 4) * outputsL[2]
                + mixMatrix.getMatrixEntryFast(4, 4) * outputsL[3];
        tmp1R =   mixMatrix.getMatrixEntryFast(0, 4) * inR
                + mixMatrix.getMatrixEntryFast(1, 4) * outputsR[0]
                + mixMatrix.getMatrixEntryFast(2, 4) * outputsR[1]
                + mixMatrix.getMatrixEntryFast(3, 4) * outputsR[2]
                + mixMatrix.getMatrixEntryFast(4, 4) * outputsR[3];


        inOutL[n] = (dryFactor*inL + wetFactor*tmp1L);
        inOutR[n] = (dryFactor*inR + wetFactor*tmp1R);
      }
    }
    break;

  default:
    {
      for(n=0; n<numFrames; n++)
      {
        inOutL[n] = 0.0;
        inOutR[n] = 0.0;
        //DEBUG_BREAK;
      }
    }
  }

  mutex.unlock();
}

/*
void Quadrifex::getSampleFrameStereo(double *inOutL, double *inOutR)
{
  double tmp1L = 0.0;
  double tmp1R = 0.0;
  double tmp2L = 0.0;
  double tmp2R = 0.0;
  double tmp3L = 0.0;
  double tmp3R = 0.0;
  double tmp4L = 0.0;
  double tmp4R = 0.0;

  double inL   = *inOutL;
  double inR   = *inOutR;

  switch( slotRouting )
  {
  case R_1TO2TO3TO4:
    {
      applySlotEffect(0, &inL,   &inR,   &tmp1L, &tmp1R, false);
      applySlotEffect(1, &tmp1L, &tmp1R, &tmp1L, &tmp1R, false);
      applySlotEffect(2, &tmp1L, &tmp1R, &tmp1L, &tmp1R, false);
      applySlotEffect(3, &tmp1L, &tmp1R, &tmp1L, &tmp1R, false);
    }
    break;
  case R_1TO2TO3_PLUS4:
    {
      applySlotEffect(0, &inL,   &inR,   &tmp1L, &tmp1R, false);
      applySlotEffect(1, &tmp1L, &tmp1R, &tmp1L, &tmp1R, false);
      applySlotEffect(2, &tmp1L, &tmp1R, &tmp1L, &tmp1R, false);
      applySlotEffect(3, &inL,   &inR,   &tmp2L, &tmp2R, false);
      tmp1L += tmp2L;
      tmp1R += tmp2R;
    }
    break;
  case R_1TO2_PLUS_3TO4:
    {
      applySlotEffect(0, &inL,   &inR,   &tmp1L, &tmp1R,  false);
      applySlotEffect(1, &tmp1L, &tmp1R, &tmp1L, &tmp1R,  false);
      applySlotEffect(2, &inL,   &inR,   &tmp2L, &tmp2R,  false);
      applySlotEffect(3, &tmp2L, &tmp2R, &tmp2L, &tmp2R,  false);
      tmp1L += tmp2L;
      tmp1R += tmp2R;
    }
    break;
  case R_1PLUS2PLUS3PLUS4:
    {
      applySlotEffect(0, &inL, &inR, &tmp1L, &tmp1R, false);
      applySlotEffect(1, &inL, &inR, &tmp2L, &tmp2R, false);
      applySlotEffect(2, &inL, &inR, &tmp3L, &tmp3R, false);
      applySlotEffect(3, &inL, &inR, &tmp4L, &tmp4R, false);
      tmp1L += tmp2L+tmp3L+tmp4L;
      tmp1R += tmp2R+tmp3R+tmp4R;
    }
    break;
  case R_1PLUS2PLUS3_TO_4:
    {
      applySlotEffect(0, &inL,   &inR,   &tmp1L, &tmp1R, false);
      applySlotEffect(1, &inL,   &inR,   &tmp2L, &tmp2R, false);
      applySlotEffect(2, &inL,   &inR,   &tmp3L, &tmp3R, false);
      tmp1L += tmp2L+tmp3L;
      tmp1R += tmp2R+tmp3R;
      applySlotEffect(3, &tmp1L, &tmp1R, &tmp1L, &tmp1R, false);
    }
    break;
  case R_1_TO_2PLUS3_TO_4:
    {
      applySlotEffect(0, &inL,   &inR,   &tmp1L, &tmp1R, false);
      applySlotEffect(1, &tmp1L, &tmp1R, &tmp2L, &tmp2R, false);
      applySlotEffect(2, &tmp1L, &tmp1R, &tmp3L, &tmp3R, false);
      tmp1L = tmp2L+tmp3L;
      tmp1R = tmp2R+tmp3R;
      applySlotEffect(3, &tmp1L, &tmp1R, &tmp1L, &tmp1R, false);
    }
    break;
  case R_1PLUS2_TO_3TO4:
    {
      applySlotEffect(0, &inL,   &inR,   &tmp1L, &tmp1R, false);
      applySlotEffect(1, &inL,   &inR,   &tmp2L, &tmp2R, false);
      tmp1L += tmp2L;
      tmp1R += tmp2R;
      applySlotEffect(2, &tmp1L, &tmp1R, &tmp1L, &tmp1R, false);
      applySlotEffect(3, &tmp1L, &tmp1R, &tmp1L, &tmp1R, false);
    }
    break;
  case R_1TO2_TO_3PLUS4:
    {
      applySlotEffect(0, &inL,   &inR,   &tmp1L, &tmp1R, false);
      applySlotEffect(1, &tmp1L, &tmp1R, &tmp1L, &tmp1R, false);
      applySlotEffect(2, &tmp1L, &tmp1R, &tmp2L, &tmp2R, false);
      applySlotEffect(3, &tmp1L, &tmp1R, &tmp3L, &tmp3R, false);
      tmp1L = tmp2L + tmp3L;
      tmp1R = tmp2R + tmp3R;
    }
    break;
  case R_1PLUS2_TO_3PLUS4:
    {
      applySlotEffect(0, &inL,   &inR,   &tmp1L, &tmp1R, false);
      applySlotEffect(1, &inL,   &inR,   &tmp2L, &tmp2R, false);
      tmp1L += tmp2L;
      tmp1R += tmp2R;
      applySlotEffect(2, &tmp1L, &tmp1R, &tmp2L, &tmp2R, false);
      applySlotEffect(3, &tmp1L, &tmp1R, &tmp3L, &tmp3R, false);
      tmp1L = tmp2L+tmp3L;
      tmp1R = tmp2R+tmp3R;
    }
    break;
  case R_1_TO_2PLUS3PLUS4:
    {
      applySlotEffect(0, &inL,   &inR,   &tmp1L, &tmp1R, false);
      applySlotEffect(1, &tmp1L, &tmp1R, &tmp2L, &tmp2R, false);
      applySlotEffect(2, &tmp1L, &tmp1R, &tmp3L, &tmp3R, false);
      applySlotEffect(3, &tmp1L, &tmp1R, &tmp4L, &tmp4R, false);
      tmp1L = tmp2L+tmp3L+tmp4L;
      tmp1R = tmp2R+tmp3R+tmp4L;
    }
    break;
  }

  *inOutL = inL * dryFactor  + tmp1L * wetFactor;
  *inOutR = inR * dryFactor  + tmp1R * wetFactor;
}
*/