#include "Equalizer5.h"

//----------------------------------------------------------------------------
// construction/destruction:

Equalizer5::Equalizer5()
{
 // init member variables:
 sampleRate = 44100.0;
 stereoMode = STEREO_LINKED;

 volL       = 1.0;
 volR       = 1.0;

 hpfL.setNumStages(4);
 hpfR.setNumStages(4);
 eqChainL.setNumStages(5);
 eqChainR.setNumStages(5);
 lpfL.setNumStages(4);
 lpfR.setNumStages(4);

 resetLpfCoeffsL();
 resetLpfCoeffsR();
 resetEqCoeffsL();
 resetEqCoeffsR();
 resetHpfCoeffsL();
 resetHpfCoeffsR();

 int64A s, c; // indices for stage and channel

 // init biquad-coefficients:
 for(s=0; s<numStages; s++)
 {
  b0eqL[s] = 1.0;
  b1eqL[s] = 0.0;
  b2eqL[s] = 0.0;
  a1eqL[s] = 0.0;
  a2eqL[s] = 0.0;
  b0eqR[s] = 1.0;
  b1eqR[s] = 0.0;
  b2eqR[s] = 0.0;
  a1eqR[s] = 0.0;
  a2eqR[s] = 0.0;
 }
 // pass the coefficients to the "BiquadCasccade"-objects:
 eqChainL.setCoeffs(b0eqL, b1eqL, b2eqL, a1eqL, a2eqL);
 eqChainR.setCoeffs(b0eqR, b1eqR, b2eqR, a1eqR, a2eqR);

 for(c=0; c<2; c++)
 {
  hpfDesigners[c].setMode(IirDesigner::HIGHPASS);
  hpfDesigners[c].setSampleRate(sampleRate);
  hpfDesigners[c].setSlope(2);
  hpfDesigners[c].setFreq1(100.0);

  lpfDesigners[c].setMode(IirDesigner::LOWPASS);
  lpfDesigners[c].setSampleRate(sampleRate);
  lpfDesigners[c].setSlope(2);
  lpfDesigners[c].setFreq1(1000.0);
 }

 for(s=1; s<numStages-1; s++)
 {
  for(c=0; c<2; c++)
  {
   eqMode[s][c] = PEAK;
   eqFreq[s][c] = 1000.0;
   eqGain[s][c] = 0.0;
   eqQ[s][c]    = SQRT2_INV;
   eqDesigners[s][c].setMode(BiquadDesigner::RBJ_PEAK);
   eqDesigners[s][c].setSampleRate(sampleRate);
   eqDesigners[s][c].setFreq(1000.0);
   eqDesigners[s][c].setGain(0.0);
   eqDesigners[s][c].setQ(SQRT2_INV);
  }
 }

}

Equalizer5::~Equalizer5()
{

}

//----------------------------------------------------------------------------
// parameter settings:

void Equalizer5::setSampleRate(flt64 newSampleRate)
{
 AudioModule::setSampleRate(newSampleRate);

 int64A c, s; // for indexing the channel and the filter-stage

 // update the sample-rate in the embdedded "IirDesigner" objects:
 for(c=0; c<2; c++)
 {
  hpfDesigners[c].setSampleRate(sampleRate);
  lpfDesigners[c].setMode(IirDesigner::LOWPASS);
 }

 // recalculate the biquad-coefficients:
 lpfDesigners[0].getBiquadCascadeCoeffs(b0lpL, b1lpL, b2lpL, a1lpL, a2lpL);
 lpfDesigners[1].getBiquadCascadeCoeffs(b0lpR, b1lpR, b2lpR, a1lpR, a2lpR);
 hpfDesigners[0].getBiquadCascadeCoeffs(b0hpL, b1hpL, b2hpL, a1hpL, a2hpL);
 hpfDesigners[1].getBiquadCascadeCoeffs(b0hpR, b1hpR, b2hpR, a1hpR, a2hpR);

 // pass the new biquad-coefficients to the "BiquadCascade"-objects:
 lpfL.setCoeffs(b0lpL, b1lpL, b2lpL, a1lpL, a2lpL);
 lpfR.setCoeffs(b0lpR, b1lpR, b2lpR, a1lpR, a2lpR);
 hpfL.setCoeffs(b0hpL, b1hpL, b2hpL, a1hpL, a2hpL);
 hpfR.setCoeffs(b0hpR, b1hpR, b2hpR, a1hpR, a2hpR);

 // update the sample-rate in the embdedded "BiquadDesigner" objects:
 for(s=1; s<numStages-1; s++)
 {
  for(c=0; c<2; c++)
  {
   eqDesigners[s][c].setSampleRate(sampleRate);
  }
  // let the biquad-designer with indices calculate biquad-coefficients
  // and store them at the s-th position in the respective 
  // coefficient-arrays:
  eqDesigners[s][0].getCoeffs(&(b0eqL[s]), &(b1eqL[s]), &(b2eqL[s]), 
                                           &(a1eqL[s]), &(a2eqL[s]));
  eqDesigners[s][1].getCoeffs(&(b0eqR[s]), &(b1eqR[s]), &(b2eqR[s]), 
                                           &(a1eqR[s]), &(a2eqR[s]));
 }

 // pass the new biquad-coefficients to the "BiquadCasccade"-objects:
 eqChainL.setCoeffs(b0eqL, b1eqL, b2eqL, a1eqL, a2eqL);
 eqChainR.setCoeffs(b0eqR, b1eqR, b2eqR, a1eqR, a2eqR);

 // make sure, that the filter-chains use the same set of coefficients in 
 // stereo-linked mode:
 setStereoMode(stereoMode);

 return;
}

void Equalizer5::setStereoMode(int64 newStereoMode)
{
 stereoMode = newStereoMode;

 // copy the coeffs for the left filters into the coeffs for right filters, 
 // too when we are in STERE_LINKED mode:
 if( stereoMode == STEREO_LINKED )
 {
  hpfR.setCoeffs(b0hpL, b1hpL, b2hpL, a1hpL, a2hpL);
  eqChainR.setCoeffs(b0eqL, b1eqL, b2eqL, a1eqL, a2eqL);
  lpfR.setCoeffs(b0lpL, b1lpL, b2lpL, a1lpL, a2lpL);
 }
 // in the other cases, the right-channel filters use their own sets of 
 // coefficients:
 else
 {
  hpfR.setCoeffs(b0hpR, b1hpR, b2hpR, a1hpR, a2hpR);
  eqChainR.setCoeffs(b0eqR, b1eqR, b2eqR, a1eqR, a2eqR);
  lpfR.setCoeffs(b0lpR, b1lpR, b2lpR, a1lpR, a2lpR);
 }
}

void Equalizer5::setLpfOrder(int64 newLpfOrder, int64 channel)
{
 // make sure that we don't try to access an invalid adress:
 if( channel>1 || channel<0 )
  return;

 // make sure that the order is in the valid range:
 if( newLpfOrder >= 0 && newLpfOrder <= 8 )
 {
  if( channel == 0 )
  {
   resetLpfCoeffsL();
   lpfDesigners[0].setSlope(newLpfOrder);
   lpfDesigners[0].getBiquadCascadeCoeffs(b0lpL, b1lpL, b2lpL, a1lpL, a2lpL);
   lpfL.setCoeffs(b0lpL, b1lpL, b2lpL, a1lpL, a2lpL);
   // copy the same coeffs into the right channel-LPF too, when we are in 
   // STERO_LINKED mode:
   if( stereoMode == STEREO_LINKED )
    lpfR.setCoeffs(b0lpL, b1lpL, b2lpL, a1lpL, a2lpL);
  }
  else if( channel == 1 )
  {
   resetLpfCoeffsR();
   lpfDesigners[1].setSlope(newLpfOrder);
   lpfDesigners[1].getBiquadCascadeCoeffs(b0lpR, b1lpR, b2lpR, a1lpR, a2lpR);
   lpfR.setCoeffs(b0lpR, b1lpR, b2lpR, a1lpR, a2lpR);
  }
 }
}

void Equalizer5::setLpfCutoff(flt64 newLpfCutoff, int64 channel)
{
 // make sure that we don't try to access an invalid adress:
 if( channel>1 || channel<0 )
  return;

 // make sure that the cutoff-frequency is in the valid range:
 if( newLpfCutoff>=20.0 && newLpfCutoff<=20000.0 )
 {
  if( channel == 0 )
  {
   resetLpfCoeffsL();
   lpfDesigners[0].setFreq1(newLpfCutoff);
   lpfDesigners[0].getBiquadCascadeCoeffs(b0lpL, b1lpL, b2lpL, a1lpL, a2lpL);
   lpfL.setCoeffs(b0lpL, b1lpL, b2lpL, a1lpL, a2lpL);
   // copy the same coeffs into the right channel-LPF too, when we are in 
   // STERO_LINKED mode:
   if( stereoMode == STEREO_LINKED )
    lpfR.setCoeffs(b0lpL, b1lpL, b2lpL, a1lpL, a2lpL);
  }
  else if( channel == 1 )
  {
   resetLpfCoeffsR();
   lpfDesigners[1].setFreq1(newLpfCutoff);
   lpfDesigners[1].getBiquadCascadeCoeffs(b0lpR, b1lpR, b2lpR, a1lpR, a2lpR);
   lpfR.setCoeffs(b0lpR, b1lpR, b2lpR, a1lpR, a2lpR);
  }
 }
}

void Equalizer5::setHpfOrder(int64 newHpfOrder, int64 channel)
{
 // make sure that we don't try to access an invalid adress:
 if( channel>1 || channel<0 )
  return;

 // make sure that the order is in the valid range:
 if( newHpfOrder >= 0 && newHpfOrder <= 8 )
 {
  if( channel == 0 )
  {
   resetHpfCoeffsL();
   hpfDesigners[0].setSlope(newHpfOrder);
   hpfDesigners[0].getBiquadCascadeCoeffs(b0hpL, b1hpL, b2hpL, a1hpL, a2hpL);
   hpfL.setCoeffs(b0hpL, b1hpL, b2hpL, a1hpL, a2hpL);
   // copy the same coeffs into the right channel-HPF too, when we are in 
   // STERO_LINKED mode:
   if( stereoMode == STEREO_LINKED )
    hpfR.setCoeffs(b0hpL, b1hpL, b2hpL, a1hpL, a2hpL);
  }
  else if( channel == 1 )
  {
   resetHpfCoeffsR();
   hpfDesigners[1].setSlope(newHpfOrder);
   hpfDesigners[1].getBiquadCascadeCoeffs(b0hpR, b1hpR, b2hpR, a1hpR, a2hpR);
   hpfR.setCoeffs(b0hpR, b1hpR, b2hpR, a1hpR, a2hpR);
  }
 }
}

void Equalizer5::setHpfCutoff(flt64 newHpfCutoff, int64 channel)
{
 // make sure that we don't try to access an invalid adress:
 if( channel>1 || channel<0 )
  return;

 // make sure that the cutoff-frequency is in the valid range:
 if( newHpfCutoff>=20.0 && newHpfCutoff<=20000.0 )
 {
  if( channel == 0 )
  {
   resetHpfCoeffsL();
   hpfDesigners[0].setFreq1(newHpfCutoff);
   hpfDesigners[0].getBiquadCascadeCoeffs(b0hpL, b1hpL, b2hpL, a1hpL, a2hpL);
   hpfL.setCoeffs(b0hpL, b1hpL, b2hpL, a1hpL, a2hpL);
   // copy the same coeffs into the right channel-HPF too, when we are in 
   // STERO_LINKED mode:
   if( stereoMode == STEREO_LINKED )
    hpfR.setCoeffs(b0hpL, b1hpL, b2hpL, a1hpL, a2hpL);
  }
  else if( channel == 1 )
  {
   resetHpfCoeffsR();
   hpfDesigners[1].setFreq1(newHpfCutoff);
   hpfDesigners[1].getBiquadCascadeCoeffs(b0hpR, b1hpR, b2hpR, a1hpR, a2hpR);
   hpfR.setCoeffs(b0hpR, b1hpR, b2hpR, a1hpR, a2hpR);
  }
 }
}

void Equalizer5::setEqMode(int64 newMode, int64 stage, int64 channel)
{
 static int64A s, c;
 s = stage;
 c = channel;

 // make sure that we don't try to access an invalid adress:
 if( stage > (numStages-1) || channel > 1 )
  return;

 eqMode[s][c] = newMode;

 // "translate" the mode index of the Equalizer5-class into the corresponding
 // mode index of the BiquadDesigner-class:
 static int64A translatedMode;
 switch(newMode)
 {
  case equalizerModes::BYPASS:
   translatedMode = BiquadDesigner::modes::BYPASS;
  break;
  case equalizerModes::PEAK:
   translatedMode = BiquadDesigner::modes::RBJ_PEAK;
  break;
  case equalizerModes::ONE_POLE_LOW_SHELV:
   translatedMode = BiquadDesigner::modes::ONE_POLE_LOW_SHELV;
  break;
  case equalizerModes::TWO_POLE_LOW_SHELV:
   translatedMode = BiquadDesigner::modes::RBJ_LOW_SHELV;
  break;
  case equalizerModes::ONE_POLE_HIGH_SHELV:
   translatedMode = BiquadDesigner::modes::ONE_POLE_HIGH_SHELV;
  break;
  case equalizerModes::TWO_POLE_HIGH_SHELV:
   translatedMode = BiquadDesigner::modes::RBJ_HIGH_SHELV;
  break;
  case equalizerModes::NOTCH:
   translatedMode = BiquadDesigner::modes::RBJ_BRF;
  break;
 } // end of  switch(newMode)

 // set the mode in the appropriate biquad-designer:
 eqDesigners[s][c].setMode(translatedMode);

 // two cases for the two channels:
 if( c == 0 )
 {
  // let the biquad-designer with indices (s,0) calculate biquad-coefficients
  // and store them at the s-th position in the respective 
  // coefficient-arrays for the left channel:
  eqDesigners[s][c].getCoeffs(&(b0eqL[s]), &(b1eqL[s]), &(b2eqL[s]), 
                                           &(a1eqL[s]), &(a2eqL[s]));

  // pass the new coefficient-arrays to the left biquad cascade:
  eqChainL.setCoeffs(b0eqL, b1eqL, b2eqL, a1eqL, a2eqL);

  // copy the same coeffs into the right channel biquad cascade too, when we
  // are in  STERO_LINKED mode:
  if( stereoMode == STEREO_LINKED )
   eqChainR.setCoeffs(b0eqL, b1eqL, b2eqL, a1eqL, a2eqL);
 }
 else if ( c == 1 )
 {
  eqDesigners[s][c].getCoeffs(&(b0eqR[s]), &(b1eqR[s]), &(b2eqR[s]), 
                                           &(a1eqR[s]), &(a2eqR[s]));
  eqChainR.setCoeffs(b0eqR, b1eqR, b2eqR, a1eqR, a2eqR);
 }

}

void Equalizer5::setEqFreq(flt64 newFreq, int64 stage, int64 channel)
{
 static int64A s, c;
 s = stage;
 c = channel;

 // make sure that we don't try to access an invalid adress:
 if( stage > (numStages-1) || channel > 1 )
  return;

 if( newFreq >= 20.0 && newFreq <= 20000.0 )
 {
  eqFreq[s][c] = newFreq;
  eqDesigners[s][c].setFreq(newFreq);
  if( c == 0 )
  {
   // let the biquad-designer with indices (s,0) calculate biquad-coefficients
   // and store them at the s-th position in the respective 
   // coefficient-arrays for the left channel:
   eqDesigners[s][c].getCoeffs(&(b0eqL[s]), &(b1eqL[s]), &(b2eqL[s]), 
                                            &(a1eqL[s]), &(a2eqL[s]));

   // pass the new coefficient-arrays to the left biquad cascade:
   eqChainL.setCoeffs(b0eqL, b1eqL, b2eqL, a1eqL, a2eqL);

  // copy the same coeffs into the right channel biquad cascade too, when we
  // are in  STERO_LINKED mode:
  if( stereoMode == STEREO_LINKED )
   eqChainR.setCoeffs(b0eqL, b1eqL, b2eqL, a1eqL, a2eqL);
  }
  else if ( c == 1 )
  {
   eqDesigners[s][c].getCoeffs(&(b0eqR[s]), &(b1eqR[s]), &(b2eqR[s]), 
                                            &(a1eqR[s]), &(a2eqR[s]));
   eqChainR.setCoeffs(b0eqR, b1eqR, b2eqR, a1eqR, a2eqR);
  }

 }

 //float test = 1 / (2-stage); // should trigger a crash when stage==2

}

void Equalizer5::setEqGain(flt64 newGain, int64 stage, int64 channel)
{
 static int64A s, c;
 s = stage;
 c = channel;

 // make sure that we don't try to access an invalid adress:
 if( stage > (numStages-1) || channel > 1 )
  return;

 eqGain[s][c] = newGain;
 eqDesigners[s][c].setGain(newGain);
 if( c == 0 )
 {
  eqDesigners[s][c].getCoeffs(&(b0eqL[s]), &(b1eqL[s]), &(b2eqL[s]), 
                                           &(a1eqL[s]), &(a2eqL[s]));
  eqChainL.setCoeffs(b0eqL, b1eqL, b2eqL, a1eqL, a2eqL);
  if( stereoMode == STEREO_LINKED )
   eqChainR.setCoeffs(b0eqL, b1eqL, b2eqL, a1eqL, a2eqL);
 }
 else if ( c == 1 )
 {
  eqDesigners[s][c].getCoeffs(&(b0eqR[s]), &(b1eqR[s]), &(b2eqR[s]), 
                                           &(a1eqR[s]), &(a2eqR[s]));
  eqChainR.setCoeffs(b0eqR, b1eqR, b2eqR, a1eqR, a2eqR);
 }

}

void Equalizer5::setEqQ(flt64 newQ, int64 stage, int64 channel)
{
 static int64A s, c;
 s = stage;
 c = channel;

 // make sure that we don't try to access an invalid adress:
 if( stage > (numStages-1) || channel > 1 )
  return;

 if( newQ > 0.00001 )
 {
  eqQ[s][c] = newQ;
  eqDesigners[s][c].setQ(newQ);
  if( c == 0 )
  {
   // let the biquad-designer with indices (s,0) calculate biquad-coefficients
   // and store them at the s-th position in the respective 
   // coefficient-arrays for the left channel:
   eqDesigners[s][c].getCoeffs(&(b0eqL[s]), &(b1eqL[s]), &(b2eqL[s]), 
                                            &(a1eqL[s]), &(a2eqL[s]));

   // pass the new coefficient-arrays to the left biquad cascade:
   eqChainL.setCoeffs(b0eqL, b1eqL, b2eqL, a1eqL, a2eqL);
   if( stereoMode == STEREO_LINKED )
    eqChainR.setCoeffs(b0eqL, b1eqL, b2eqL, a1eqL, a2eqL);
  }
  else if ( c == 1 )
  {
   eqDesigners[s][c].getCoeffs(&(b0eqR[s]), &(b1eqR[s]), &(b2eqR[s]), 
                                            &(a1eqR[s]), &(a2eqR[s]));
   eqChainR.setCoeffs(b0eqR, b1eqR, b2eqR, a1eqR, a2eqR);
  }
 }
}

void Equalizer5::setLevel(flt64 newLevel, int64 channel)
{
 if( channel == 0 )
  volL = dB2amp(newLevel);
 else if ( channel == 1 )
  volR = dB2amp(newLevel);
}


//----------------------------------------------------------------------------
// others:

void Equalizer5::getMagnitudeResponse(float *omegas, float *magnitudes, 
                                      int numBins, int channel)
{
 static int64A k;
 static flt64A dBOffsetL, dBOffsetR;

 // the global offsets due to the global gain-factors:
 dBOffsetL = (float) amp2dB(volL);
 dBOffsetR = (float) amp2dB(volR);

 switch( channel )
 {
  case 0:
  {
   // accumulate the magnitude responses of the filters:
   hpfL.getMagnitudeResponse(omegas, magnitudes, numBins, true, false);
   eqChainL.getMagnitudeResponse(omegas, magnitudes, numBins, true, true);
   lpfL.getMagnitudeResponse(omegas, magnitudes, numBins, true, true);
   // apply the global gain value:
   for(k=0; k<numBins; k++)
    magnitudes[k] += dBOffsetL;

   // accumulate the magnitude response of the second filter-chain also, 
   // when we are in MONO_10-mode (in this mode the equalizers are chained and 
   // work both on the left input signal:
   if( stereoMode == MONO_10 )
   {
    hpfR.getMagnitudeResponse(omegas, magnitudes, numBins, true, true);
    eqChainR.getMagnitudeResponse(omegas, magnitudes, numBins, true, true);
    lpfR.getMagnitudeResponse(omegas, magnitudes, numBins, true, true);
    for(k=0; k<numBins; k++)
     magnitudes[k] += dBOffsetR;
   }
  }
  break;
  case 1:
  {
   hpfR.getMagnitudeResponse(omegas, magnitudes, numBins, true, false);
   eqChainR.getMagnitudeResponse(omegas, magnitudes, numBins, true, true);
   lpfR.getMagnitudeResponse(omegas, magnitudes, numBins, true, true);
   for(k=0; k<numBins; k++)
    magnitudes[k] += dBOffsetR;
  }
  break;
 }

}

void Equalizer5::resetHpfCoeffsL()
{
 static int64A k;
 for(k=0; k<numStages; k++)
 {
  b0hpL[k] = 1.0;
  b1hpL[k] = 0.0;
  b2hpL[k] = 0.0;
  a1hpL[k] = 0.0;
  a2hpL[k] = 0.0;
 }
}
void Equalizer5::resetHpfCoeffsR()
{
 static int64A k;
 for(k=0; k<numStages; k++)
 {
  b0hpR[k] = 1.0;
  b1hpR[k] = 0.0;
  b2hpR[k] = 0.0;
  a1hpR[k] = 0.0;
  a2hpR[k] = 0.0;
 }
}
void Equalizer5::resetEqCoeffsL()
{
 static int64A k;
 for(k=0; k<numStages; k++)
 {
  b0eqL[k] = 1.0;
  b1eqL[k] = 0.0;
  b2eqL[k] = 0.0;
  a1eqL[k] = 0.0;
  a2eqL[k] = 0.0;
 }
}
void Equalizer5::resetEqCoeffsR()
{
 static int64A k;
 for(k=0; k<numStages; k++)
 {
  b0eqR[k] = 1.0;
  b1eqR[k] = 0.0;
  b2eqR[k] = 0.0;
  a1eqR[k] = 0.0;
  a2eqR[k] = 0.0;
 }
}
void Equalizer5::resetLpfCoeffsL()
{
 static int64A k;
 for(k=0; k<numStages; k++)
 {
  b0lpL[k] = 1.0;
  b1lpL[k] = 0.0;
  b2lpL[k] = 0.0;
  a1lpL[k] = 0.0;
  a2lpL[k] = 0.0;
 }
}
void Equalizer5::resetLpfCoeffsR()
{
 static int64A k;
 for(k=0; k<numStages; k++)
 {
  b0lpR[k] = 1.0;
  b1lpR[k] = 0.0;
  b2lpR[k] = 0.0;
  a1lpR[k] = 0.0;
  a2lpR[k] = 0.0;
 }
}