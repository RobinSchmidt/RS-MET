#include "HighOrderEqualizer.h"

//----------------------------------------------------------------------------
// construction/destruction:

HighOrderEqualizer::HighOrderEqualizer()
{
 // init member variables:
 sampleRate = 44100.0;
 stereoMode = STEREO_LINKED;
 bypassAll  = false;

 int s, c, k; // indices for eq-stage, channel, and frequency-bin

 double freq, gain, q;
 int    mode, order;

 // set up the lowpass- and highpass-filters:
 for(c=0; c<numChannels; c++)
 {
  setLpfOrder(0, c);
  setLpfCutoff(20000.0, c);
  setHpfOrder(0, c);
  setHpfCutoff(20.0, c);
 }

 // set up the equalizer-chains:
 for(s=0; s<numEqStages; s++)
 {
  gain  = 1.0;
  q     = 1.0;
  order = 2;
  switch(s)
  {
  case 0:
   {
    freq = 125.0;
    mode = HighOrderEqualizer::PEAK;
    //gain  = 12.0;
   } break;
  case 1:
   {
    freq = 250.0;
    mode = HighOrderEqualizer::PEAK;
   } break;
  case 2:
   {
    freq = 1000.0;
    mode = HighOrderEqualizer::PEAK;
   } break;
  case 3:
   {
    freq = 4000.0;
    mode = HighOrderEqualizer::PEAK;
   } break;
  case 4:
   {
    freq = 8000.0;
    mode = HighOrderEqualizer::HIGH_SHELV;
   } break;
  } // end of switch(s)

  for(c=0; c<numChannels; c++)
  {
   setEqMode(mode, s, c);
   setEqOrder(order, s, c);
   setEqFreq(freq, s, c);
   setEqGain(gain, s, c);
   setEqQ(q, s, c);
  }
 }

 // set up the global volume factors:
 vol[0]     = 1.0;
 vol[1]     = 1.0;

 // init the arrays for the magnitude response plot:
 double fMin =    15.625;
 double fMax = 32000.0;
 double f;          // current frequency;
 for(k=0; k<numBins; k++)
 {
  // calculate current frequency: 
  f =  mapLinearToExponential(k, 0, numBins-1, fMin, fMax);

  // calculate some dependent quantities and store them into their arrays:
  freqs[k]             = f;
  floatFrequencies[k] = (float) f;
  omegas[k]           = 2*PI*f/sampleRate;
  cosOmegas[k]        = cos(omegas[k]);
  cos2Omegas[k]       = cos(2*omegas[k]);

  // init the individual magnitude-arrays with 0.0 dB:
  for(c=0; c<numChannels; c++)
  {
   hpfResponses[c][k]    = 0.0;
   lpfResponses[c][k]    = 0.0;
   floatMagnitudes[c][k] = 0.0f;
   for(s=0; s<numEqStages; s++)
    eqResponses[s][c][k] = 0.0;
  }
 }
}

HighOrderEqualizer::~HighOrderEqualizer()
{

}

//----------------------------------------------------------------------------
// parameter settings:

void HighOrderEqualizer::setSampleRate(double newSampleRate)
{
 intA c, s, k; // for indexing the eq-stage, channel and bin

 AudioModule::setSampleRate(newSampleRate); // our member "sampleRate" is up 
                                            // to date now

 // tell the embedded IirDesigner-object the new sample-rate:
 designer.setSampleRate(sampleRate);

 // re-calculate coefficients for all the embedded BiquadCascade-objects:
 for(c=0; c<numChannels; c++)
 {
  updateHpfCoeffs(c);
  for(s=0; s<numEqStages; s++)
   updateEqCoeffs(s,c);
  updateLpfCoeffs(c);
 }

 // re-calculate the omega-arrays for the magnitude response plot:
 double f;          // for the current frequency
 for(k=0; k<numBins; k++)
 {
  // get the current frequency from the freq-array: 
  f =  freqs[k];

  // calculate some dependent quantities and store them into their arrays:
  omegas[k]           = 2*PI*f/sampleRate;
  cosOmegas[k]        = cos(omegas[k]);
  cos2Omegas[k]       = cos(2*omegas[k]);
 }

 // update the individual magnitude-arrays:
 for(c=0; c<numChannels; c++)
 {
  updateHpfFreqResponse(c);
  for(s=0; s<numEqStages; s++)
   updateEqFreqResponse(s,c);
  updateLpfFreqResponse(c);

  updateOverallFreqResponse(c);
 }

 // make sure, that the filter-chains use the same set of coefficients in 
 // stereo-linked mode:
 setStereoMode(stereoMode);
}

void HighOrderEqualizer::setStereoMode(int newStereoMode)
{
 intA slope, s;

 hpf[1].initBiquadCoeffs();
 for(s=0; s<numEqStages; s++)
  eq[s][1].initBiquadCoeffs();
 lpf[1].initBiquadCoeffs();

 if( newStereoMode == STEREO_LINKED )
 {
  // tell the BiquadCascade-objects how many biquad stages are to be used:
  slope = hpfOrder[0];
  if( isEven(slope) )
   hpf[1].setNumStages( slope/2 );
  else
   hpf[1].setNumStages( (slope+1)/2 );
  for(s=0; s<numEqStages; s++)
  {
   slope = eqOrder[s][0];
   if( eqMode[s][0] == LOW_SHELV || eqMode[s][0] == HIGH_SHELV )
    // filter order is the same as the slope-value - half as many 
    // biquad-stages are required (up to a left-over one-pole stage):
   {
    if( isEven(slope) )
     eq[s][1].setNumStages( slope/2 );
    else
     eq[s][1].setNumStages( (slope+1)/2 );
   }
   else if(eqMode[s][0] == PEAK || eqMode[s][0] == NOTCH)
    // filter order is twice the slope-value - just as many biquad-stages 
    // are required:
    eq[s][1].setNumStages( slope );
   else if(eqMode[s][0] == BYPASS)
    eq[s][1].setNumStages( 0 );
  }
  slope = lpfOrder[0];
  if( isEven(slope) )
   lpf[1].setNumStages( slope/2 );
  else
   lpf[1].setNumStages( (slope+1)/2 );

  // pass the coeffs for the left filters to the BiquadCascade-objects for
  // right filters when we are in STEREO_LINKED mode:
  hpf[1].setCoeffs(&(b0hp[0][0]), 
                   &(b1hp[0][0]), 
                   &(b2hp[0][0]),
                   &(a1hp[0][0]), 
                   &(a2hp[0][0])  );

  for(s=0; s<numEqStages; s++)
  {
   eq[s][1].setCoeffs(&(b0eq[s][0][0]), 
                      &(b1eq[s][0][0]), 
                      &(b2eq[s][0][0]),
                      &(a1eq[s][0][0]), 
                      &(a2eq[s][0][0])  );
  }

  lpf[1].setCoeffs(&(b0lp[0][0]), 
                   &(b1lp[0][0]), 
                   &(b2lp[0][0]),
                   &(a1lp[0][0]), 
                   &(a2lp[0][0])  );
 }
 // in the other cases, the right-channel filters use their own sets of 
 // coefficients:
 else
 {
  // tell the BiquadCascade-objects how many biquad stages are to be used:
  slope = hpfOrder[1];
  if( isEven(slope) )
   hpf[1].setNumStages( slope/2 );
  else
   hpf[1].setNumStages( (slope+1)/2 );
  for(s=0; s<numEqStages; s++)
  {
   slope = eqOrder[s][1];
   if( eqMode[s][1] == LOW_SHELV || eqMode[s][1] == HIGH_SHELV )
    // filter order is the same as the slope-value - half as many 
    // biquad-stages are required (up to a left-over one-pole stage):
   {
    if( isEven(slope) )
     eq[s][1].setNumStages( slope/2 );
    else
     eq[s][1].setNumStages( (slope+1)/2 );
   }
   else if(eqMode[s][1] == PEAK || eqMode[s][1] == NOTCH)
    // filter order is twice the slope-value - just as many biquad-stages 
    // are required:
    eq[s][1].setNumStages( slope );
   else if(eqMode[s][1] == BYPASS)
    eq[s][1].setNumStages( 0 );
  }
  slope = lpfOrder[1];
  if( isEven(slope) )
   lpf[1].setNumStages( slope/2 );
  else
   lpf[1].setNumStages( (slope+1)/2 );

  // pass the coeffs for the right filters to the BiquadCascade-objects for
  // right filters when we are not in STEREO_LINKED mode:
  hpf[1].setCoeffs(&(b0hp[1][0]), 
                   &(b1hp[1][0]), 
                   &(b2hp[1][0]),
                   &(a1hp[1][0]), 
                   &(a2hp[1][0])  );

  for(s=0; s<numEqStages; s++)
  {
   eq[s][1].setCoeffs(&(b0eq[s][1][0]), 
                      &(b1eq[s][1][0]), 
                      &(b2eq[s][1][0]),
                      &(a1eq[s][1][0]), 
                      &(a2eq[s][1][0])  );
  }

  lpf[1].setCoeffs(&(b0lp[1][0]), 
                   &(b1lp[1][0]), 
                   &(b2lp[1][0]),
                   &(a1lp[1][0]), 
                   &(a2lp[1][0])  );
 }


 // update the frequency-response curves, if the stereo-mode has been changed:
 if( newStereoMode != stereoMode )
 {
  stereoMode = newStereoMode;
  updateOverallFreqResponse(0);
 }

}

void HighOrderEqualizer::setLpfOrder(int newLpfOrder, int channel)
{
 static intA c;
 c = channel;

 // make sure that we don't try to access an invalid adress:
 if( c<0 || c>(numChannels-1) )
  return;
 else if( newLpfOrder >= 0 && newLpfOrder <= 2*maxBiquadsInFilter )
  lpfOrder[c] = newLpfOrder;

 updateLpfCoeffs(c);
 setStereoMode(stereoMode);  // makes sure that the new coeffs are passed to the 
                             // right-channel filters in stereo-linked mode
 updateLpfFreqResponse(c);
 updateOverallFreqResponse(c);
}

void HighOrderEqualizer::setLpfCutoff(double newLpfCutoff, int channel)
{
 static intA c;
 c = channel;

 // make sure that we don't try to access an invalid adress:
 if( c<0 || c>(numChannels-1) )
  return;
 else if( newLpfCutoff >= 20.0 && newLpfCutoff <= 20000.0 )
  lpfCutoff[c] = newLpfCutoff;

 updateLpfCoeffs(c);
 setStereoMode(stereoMode);
 updateLpfFreqResponse(c);
 updateOverallFreqResponse(c);
}

void HighOrderEqualizer::setHpfOrder(int newHpfOrder, int channel)
{
 static intA c;
 c = channel;

 // make sure that we don't try to access an invalid adress:
 if( c<0 || c>(numChannels-1) )
  return;
 else if( newHpfOrder >= 0 && newHpfOrder <= 2*maxBiquadsInFilter )
  hpfOrder[c] = newHpfOrder;

 updateHpfCoeffs(c);
 setStereoMode(stereoMode);
 updateHpfFreqResponse(c);
 updateOverallFreqResponse(c);
}

void HighOrderEqualizer::setHpfCutoff(double newHpfCutoff, int channel)
{
 static intA c;
 c = channel;

 // make sure that we don't try to access an invalid adress:
 if( c<0 || c>(numChannels-1) )
  return;
 else if( newHpfCutoff >= 20.0 && newHpfCutoff <= 20000.0 )
  hpfCutoff[c] = newHpfCutoff;

 updateHpfCoeffs(c);
 setStereoMode(stereoMode);
 updateHpfFreqResponse(c);
 updateOverallFreqResponse(c);
}

void HighOrderEqualizer::setEqMode(int newMode, int stage, int channel)
{
 static intA s, c;
 s = stage;
 c = channel;

 // make sure that we don't try to access an invalid adress:
 if(  s<0 || s>(numEqStages-1) || c<0 || c>(numChannels-1) )
  return;
 else if( newMode >= HighOrderEqualizer::BYPASS && 
          newMode <= HighOrderEqualizer::HIGH_SHELV )
 {
  eqMode[s][c] = newMode;
 }

 updateEqCoeffs(s,c);
 setStereoMode(stereoMode);
 updateEqFreqResponse(s,c);
 updateOverallFreqResponse(c);
}

void HighOrderEqualizer::setEqOrder(int newOrder, int stage, int channel)
{
 static intA s, c;
 s = stage;
 c = channel;

 // make sure that we don't try to access an invalid adress:
 if(  s<0 || s>(numEqStages-1) || c<0 || c>(numChannels-1) )
  return;
 else if( newOrder >= 0 && newOrder <= maxBiquadsInEq )
  eqOrder[s][c] = newOrder;

 updateEqCoeffs(s,c);
 setStereoMode(stereoMode);
 updateEqFreqResponse(s,c);
 updateOverallFreqResponse(c);
}

void HighOrderEqualizer::setEqFreq(double newFreq, int stage, int channel)
{
 static intA s, c;
 s = stage;
 c = channel;

 // make sure that we don't try to access an invalid adress:
 if(  s<0 || s>(numEqStages-1) || c<0 || c>(numChannels-1) )
  return;
 else if( newFreq >= 20.0 && newFreq <= 20000.0 )
  eqFreq[s][c] = newFreq;

 updateEqCoeffs(s,c);
 setStereoMode(stereoMode);
 updateEqFreqResponse(s,c);
 updateOverallFreqResponse(c);
}

void HighOrderEqualizer::setEqGain(double newGain, int stage, int channel)
{
 static intA s, c;
 s = stage;
 c = channel;

 // make sure that we don't try to access an invalid adress:
 if(  s<0 || s>(numEqStages-1) || c<0 || c>(numChannels-1) )
  return;
 else
  eqGain[s][c] = dB2amp(newGain);

 updateEqCoeffs(s,c);
 setStereoMode(stereoMode);
 updateEqFreqResponse(s,c);
 updateOverallFreqResponse(c);
}

void HighOrderEqualizer::setEqQ(double newQ, int stage, int channel)
{
 static intA s, c;
 s = stage;
 c = channel;

 // make sure that we don't try to access an invalid adress:
 if(  s<0 || s>(numEqStages-1) || c<0 || c>(numChannels-1) )
  return;
 else if( newQ > 0.00001 )
  eqQ[s][c] = newQ;

 updateEqCoeffs(s,c);
 setStereoMode(stereoMode);
 updateEqFreqResponse(s,c);
 updateOverallFreqResponse(c);
}

void HighOrderEqualizer::setLevel(double newLevel, int channel)
{
 if( channel == 0 )
  vol[0] = dB2amp(newLevel);
 else if ( channel == 1 )
  vol[1] = dB2amp(newLevel);

 updateOverallFreqResponse(channel);
}

void HighOrderEqualizer::setBypassAll(bool newBypassAll)
{
 bypassAll = newBypassAll;
}


//----------------------------------------------------------------------------
// others:

void HighOrderEqualizer::resetHpfCoeffs(int channel)
{
 static intA c, b; // indices for thechannel and biquad-stage 
 if( channel >= 0 && channel < numChannels)
  c = channel;
 else
  return;
 for(b=0; b<maxBiquadsInFilter; b++)
 {
  b0hp[c][b] = 1.0;
  b1hp[c][b] = 0.0;
  b2hp[c][b] = 0.0;
  a1hp[c][b] = 0.0;
  a2hp[c][b] = 0.0;
 }
}

void HighOrderEqualizer::resetEqCoeffs(int stage, int channel)
{
 static intA s, c, b;
 s = stage;
 c = channel;

 // make sure that we don't try to access an invalid adress:
 if(  s<0 || s>(numEqStages-1) || c<0 || c>(numChannels-1) )
  return;
 else
 {
  for(b=0; b<maxBiquadsInEq; b++)
  {
   b0eq[s][c][b] = 1.0;
   b1eq[s][c][b] = 0.0;
   b2eq[s][c][b] = 0.0;
   a1eq[s][c][b] = 0.0;
   a2eq[s][c][b] = 0.0;
  }
 }

 // let the embedded BiquadCascade-object reset it's internal copies of the
 // coeffcients, too:
 eq[s][c].initBiquadCoeffs();
}

void HighOrderEqualizer::resetLpfCoeffs(int channel)
{
 static intA c, b; // indices for thechannel and biquad-stage 
 if( channel >= 0 && channel < numChannels)
  c = channel;
 else
  return;
 for(b=0; b<maxBiquadsInFilter; b++)
 {
  b0lp[c][b] = 1.0;
  b1lp[c][b] = 0.0;
  b2lp[c][b] = 0.0;
  a1lp[c][b] = 0.0;
  a2lp[c][b] = 0.0;
 }
}

void HighOrderEqualizer::updateHpfCoeffs(int channel)
{
 static doubleA freq;
 static intA    slope;
 static intA    c;
 c = channel;

 // reset the highpass-coefficients to neutral values:
 resetHpfCoeffs(c);

 // select the appropriate parameter-set:
 freq  = hpfCutoff[c];
 slope = hpfOrder[c];

 // tell the BiquadCascade-object which realizes the highpass-filter how many
 // biquad stages are to be used:
 if( isEven(slope) )
  hpf[c].setNumStages( slope/2 );
 else
  hpf[c].setNumStages( (slope+1)/2 );

 // set up the filter-designer:
 designer.setMode(IirDesigner::HIGHPASS);
 designer.setFreq1(freq);
 designer.setSlope(slope);

 // let the designer calculate the equalizer-coefficients:
 designer.getBiquadCascadeCoeffs(&(b0hp[c][0]), 
                                 &(b1hp[c][0]), 
                                 &(b2hp[c][0]),
                                 &(a1hp[c][0]), 
                                 &(a2hp[c][0])  );

 // pass the coefficients to the appropriate BiquadCascade-object:
 hpf[c].setCoeffs(&(b0hp[c][0]), 
                  &(b1hp[c][0]), 
                  &(b2hp[c][0]),
                  &(a1hp[c][0]), 
                  &(a2hp[c][0])  );
}

void HighOrderEqualizer::updateEqCoeffs(int stage, int channel)
{
 static doubleA freq, lowFreq, highFreq, gain, q;
 static intA    mode, slope;
 static intA    s, c;
 s = stage;
 c = channel;

 // reset the equalizer-coefficients to neutral values:
 resetEqCoeffs(s,c);

 // select the appropriate parameter-set:
 freq  = eqFreq[s][c];
 gain  = eqGain[s][c];
 q     = eqQ[s][c];
 mode  = eqMode[s][c];
 slope = eqOrder[s][c];

 // calculate corner frequency for shelving modes or lower and upper
 // corner-frequencies for peaking- and notch-modes:

 if( mode == HighOrderEqualizer::BYPASS )
 {
  // tell the embedded IirDesigner-object the corner frequency:
  designer.setFreq1(freq);
 }
 else if( mode == HighOrderEqualizer::LOW_SHELV ||
          mode == HighOrderEqualizer::HIGH_SHELV)
 {
  // tell the embedded IirDesigner-object the corner frequency:
  designer.setFreq1(freq);
 }
 else if( mode == HighOrderEqualizer::PEAK ||
          mode == HighOrderEqualizer::NOTCH)
 {
  // calculate lower and upper corner-frequencies for peaking and 
  // notch-filters:
  lowFreq    = freq * ( sqrt(1.0+1.0/(4.0*q*q)) - 1.0/(2.0*q) );
  highFreq   = lowFreq + freq/q;

  // restrict the range of the corner-frequencies:
  if( lowFreq > 0.94*0.5*sampleRate )
   lowFreq = 0.94*0.5*sampleRate; // should actually never happen for 
                                  // sampleRates >= 44.1 kHz
  if( highFreq > 0.95*0.5*sampleRate )
   highFreq = 0.95*0.5*sampleRate; // this, however, could easily happen

  // tell the embedded IirDesigner-object the corner-frequencies:
  designer.setFreq1(lowFreq);
  designer.setFreq2(highFreq);
 }

 // "translate" the mode index of the HighOrderEqualizer-class into the 
 // corresponding mode index of the IirDesigner-class and set up the mode in
 // the embedded IirDesigner-object:
 switch( mode )
 {
 case HighOrderEqualizer::BYPASS:
  {  
   designer.setMode(IirDesigner::BYPASS);
   eq[s][c].setNumStages(0);
  }
  break;
 case HighOrderEqualizer::LOW_SHELV:
  {
   designer.setMode(IirDesigner::LOW_SHELV);
   if( isEven(slope) )
    eq[s][c].setNumStages( slope/2 );
   else
    eq[s][c].setNumStages( (slope+1)/2 );
  }
  break;
 case HighOrderEqualizer::PEAK:
  {
   designer.setMode(IirDesigner::PEAK);
   eq[s][c].setNumStages(slope);
  }
  break;
 case HighOrderEqualizer::NOTCH:
  {
   designer.setMode(IirDesigner::BANDREJECT);
   eq[s][c].setNumStages(slope);
  }
  break;
 case HighOrderEqualizer::HIGH_SHELV:
  {
   designer.setMode(IirDesigner::HIGH_SHELV);
   if( isEven(eqOrder[s][c]) )
    eq[s][c].setNumStages( slope/2 );
   else
    eq[s][c].setNumStages( (slope+1)/2 );
  }
  break;
 default:
  {
   designer.setMode(IirDesigner::BYPASS);
   eq[s][c].setNumStages(0);
  }
 } // end of switch(eqMode[stage][channel])

 // set up the order- and gain-parameters in the embedded IirDesigner-object:
 designer.setSlope(slope);
 designer.setGain(gain);

 // let the designer calculate the equalizer-coefficients:
 designer.getBiquadCascadeCoeffs(&(b0eq[s][c][0]), 
                                 &(b1eq[s][c][0]), 
                                 &(b2eq[s][c][0]),
                                 &(a1eq[s][c][0]), 
                                 &(a2eq[s][c][0])  );

 // pass the coefficients to the appropriate BiquadCascade-object:
 eq[s][c].setCoeffs(&(b0eq[s][c][0]), 
                    &(b1eq[s][c][0]), 
                    &(b2eq[s][c][0]),
                    &(a1eq[s][c][0]), 
                    &(a2eq[s][c][0])  );
}

void HighOrderEqualizer::updateLpfCoeffs(int channel)
{
 static doubleA freq;
 static intA    slope;
 static intA    c;
 c = channel;

 // reset the highpass-coefficients to neutral values:
 resetLpfCoeffs(c);

 // select the appropriate parameter-set:
 freq  = lpfCutoff[c];
 slope = lpfOrder[c];

 // tell the BiquadCascade-object which realizes the highpass-filter how many
 // biquad stages are to be used:
 if( isEven(slope) )
  lpf[c].setNumStages( slope/2 );
 else
  lpf[c].setNumStages( (slope+1)/2 );

 // set up the filter-designer:
 designer.setMode(IirDesigner::LOWPASS);
 designer.setFreq1(freq);
 designer.setSlope(slope);

 // let the designer calculate the equalizer-coefficients:
 designer.getBiquadCascadeCoeffs(&(b0lp[c][0]), 
                                 &(b1lp[c][0]), 
                                 &(b2lp[c][0]),
                                 &(a1lp[c][0]), 
                                 &(a2lp[c][0])  );

 // pass the coefficients to the appropriate BiquadCascade-object:
 lpf[c].setCoeffs(&(b0lp[c][0]), 
                  &(b1lp[c][0]), 
                  &(b2lp[c][0]),
                  &(a1lp[c][0]), 
                  &(a2lp[c][0])  );
}


void HighOrderEqualizer::updateHpfFreqResponse(int channel)
{
 static intA c, k, b; // indices for channel, bin and biquad-stage
 static doubleA num, den, accu, tmp;

 c = channel;

 // make sure that we don't try to access an invalid adress:
 if( c<0 || c>(numChannels-1) )
  return;

 for(k=0; k<numBins; k++)
 { 
  // accumulate the squared magnitude responses of the individual 
  // biquad-stages:
  accu = 1.0;
  for(b=0; b<maxBiquadsInFilter; b++)
  {
   // calculate the numerator of the squared magnitude of biquad-stage b:
   num =   b0hp[c][b] * b0hp[c][b] 
         + b1hp[c][b] * b1hp[c][b] 
         + b2hp[c][b] * b2hp[c][b]
         + 2*cosOmegas[k] * (  b0hp[c][b]*b1hp[c][b] 
                             + b1hp[c][b]*b2hp[c][b])
         + 2*cos2Omegas[k] * b0hp[c][b]*b2hp[c][b];

   // calculate the denominator of the squared magnitude of biquad-stage b:
   den =   1.0 * 1.0
         + a1hp[c][b] * a1hp[c][b] 
         + a2hp[c][b] * a2hp[c][b]
         + 2*cosOmegas[k] * (  1.0*a1hp[c][b] 
                             + a1hp[c][b]*a2hp[c][b])
         + 2*cos2Omegas[k] * 1.0*a2hp[c][b];

   // multiply the accumulator with the squared magnitude of biquad-stage b:
   accu *= (num/den);
  } // end of "for(b=0; b<maxBiquadsInEq; b++)"

  // take the square root of the accumulated squared magnitude response - this
  // is the desired magnitude of the biquad cascade at bin k:
  tmp = sqrt(accu);

  // convert this value to decibels:
  tmp = 20.0 * log10(tmp);

  // store the calculated dB-value in bin k of the freq-response of the 
  // highpass-filter for channel c:
  hpfResponses[c][k] = tmp;
 }
}

void HighOrderEqualizer::updateEqFreqResponse(int stage, int channel)
{
 static intA s, c, k, b; // indices for eq-stage, channel, bin and 
                         // biquad-stage inside the eq-stage

 static doubleA num, den, accu, tmp;

 s = stage;
 c = channel;

 // make sure that we don't try to access an invalid adress:
 if(  s<0 || s>(numEqStages-1) || c<0 || c>(numChannels-1) )
  return;

 for(k=0; k<numBins; k++)
 { 
  // accumulate the squared magnitude responses of the individual 
  // biquad-stages:
  accu = 1.0;
  for(b=0; b<maxBiquadsInEq; b++)
  {
   // calculate the numerator of the squared magnitude of biquad-stage b:
   num =   b0eq[s][c][b] * b0eq[s][c][b] 
         + b1eq[s][c][b] * b1eq[s][c][b] 
         + b2eq[s][c][b] * b2eq[s][c][b]
         + 2*cosOmegas[k] * (  b0eq[s][c][b]*b1eq[s][c][b] 
                             + b1eq[s][c][b]*b2eq[s][c][b])
         + 2*cos2Omegas[k] * b0eq[s][c][b]*b2eq[s][c][b];

   // calculate the denominator of the squared magnitude of biquad-stage b:
   den =   1.0 * 1.0
         + a1eq[s][c][b] * a1eq[s][c][b] 
         + a2eq[s][c][b] * a2eq[s][c][b]
         + 2*cosOmegas[k] * (  1.0*a1eq[s][c][b] 
                             + a1eq[s][c][b]*a2eq[s][c][b])
         + 2*cos2Omegas[k] * 1.0*a2eq[s][c][b];

   // multiply the accumulator with the squared magnitude of biquad-stage b:
   accu *= (num/den);
  } // end of "for(b=0; b<maxBiquadsInEq; b++)"

  // take the square root of the accumulated squared magnitude response - this
  // is the desired magnitude of the biquad cascade at bin k:
  tmp = sqrt(accu);

  // convert this value to decibels:
  tmp = 20.0 * log10(tmp);

  // store the calculated dB-value in bin k of the freq-response of eq-stage 
  // s for channel c:
  eqResponses[s][c][k] = tmp;
 }
}

void HighOrderEqualizer::updateLpfFreqResponse(int channel)
{
 static intA c, k, b; // indices for channel, bin and biquad-stage
 static doubleA num, den, accu, tmp;

 c = channel;

 // make sure that we don't try to access an invalid adress:
 if( c<0 || c>(numChannels-1) )
  return;

 for(k=0; k<numBins; k++)
 { 
  // accumulate the squared magnitude responses of the individual 
  // biquad-stages:
  accu = 1.0;
  for(b=0; b<maxBiquadsInFilter; b++)
  {
   // calculate the numerator of the squared magnitude of biquad-stage b:
   num =   b0lp[c][b] * b0lp[c][b] 
         + b1lp[c][b] * b1lp[c][b] 
         + b2lp[c][b] * b2lp[c][b]
         + 2*cosOmegas[k] * (  b0lp[c][b]*b1lp[c][b] 
                             + b1lp[c][b]*b2lp[c][b])
         + 2*cos2Omegas[k] * b0lp[c][b]*b2lp[c][b];

   // calculate the denominator of the squared magnitude of biquad-stage b:
   den =   1.0 * 1.0
         + a1lp[c][b] * a1lp[c][b] 
         + a2lp[c][b] * a2lp[c][b]
         + 2*cosOmegas[k] * (  1.0*a1lp[c][b] 
                             + a1lp[c][b]*a2lp[c][b])
         + 2*cos2Omegas[k] * 1.0*a2lp[c][b];

   // multiply the accumulator with the squared magnitude of biquad-stage b:
   accu *= (num/den);
  } // end of "for(b=0; b<maxBiquadsInEq; b++)"

  // take the square root of the accumulated squared magnitude response - this
  // is the desired magnitude of the biquad cascade at bin k:
  tmp = sqrt(accu);

  // convert this value to decibels:
  tmp = 20.0 * log10(tmp);

  // store the calculated dB-value in bin k of the freq-response of the 
  // highpass-filter for channel c:
  lpfResponses[c][k] = tmp;
 }
}

void HighOrderEqualizer::updateOverallFreqResponse(int channel)
{
 static intA    s, c, k;
 static doubleA accu;

 c = channel;

 // make sure that we don't try to access an invalid adress:
 if( c<0 || c>(numChannels-1) )
  return;

 for(k=0; k<numBins; k++)
 { 
  if( stereoMode != MONO_10 ) 
  {
   accu  = amp2dB(vol[c]);
   accu += hpfResponses[c][k];
   accu += lpfResponses[c][k];
   for(s=0; s<numEqStages; s++)
    accu += eqResponses[s][c][k];
   // typecast and store:
   floatMagnitudes[c][k] = (float) accu;
  }
  else if( stereoMode == MONO_10 ) 
   // we need to add both channels freq responses in this case
  {
   accu  = amp2dB(vol[0]);
   accu += amp2dB(vol[1]);
   accu += hpfResponses[0][k];
   accu += hpfResponses[1][k];
   accu += lpfResponses[0][k];
   accu += lpfResponses[1][k];
   for(s=0; s<numEqStages; s++)
   {
    accu += eqResponses[s][0][k];
    accu += eqResponses[s][1][k];
   }
   // typecast and store:
   floatMagnitudes[0][k] = (float) accu;
  }
 } // end of for(k=0; k<numBins; k++)
}

void HighOrderEqualizer::getMagnitudeResponse(int channel, int* arrayLengths, 
                                              float** freqArray, float** magArray)
{
 // assign the number of bins to the output-slot "arrayLengths":
 *arrayLengths = numBins;

 // assign the output-slot *freqArray to the single-precision frequency-array:
 *freqArray = floatFrequencies;

 // assign the output-slot *magArray to the single-precision magnitude-array 
 // for the requested channel:
 *magArray = &(floatMagnitudes[channel][0]);
}