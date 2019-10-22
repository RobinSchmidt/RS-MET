#include "MagicCarpetEqualizer.h"

//----------------------------------------------------------------------------
// construction/destruction:

MagicCarpetEqualizer::MagicCarpetEqualizer()
{
 // init member variables:
 sampleRate = 44100.0;
 isOff      = false;

 eqChainL.setNumStages(3);
 eqChainR.setNumStages(3);

 resetEqCoeffs();

 intA s; // indices for stage and channel

 // init biquad-coefficients:
 for(s = 0; s < numStages; s++)
 {
  b0[s] = 1.0;
  b1[s] = 0.0;
  b2[s] = 0.0;
  a1[s] = 0.0;
  a2[s] = 0.0;
 }
 // pass the coefficients to the "BiquadCasccade"-objects:
 eqChainL.setCoeffs(b0, b1, b2, a1, a2);
 eqChainR.setCoeffs(b0, b1, b2, a1, a2);

 for(s = 0; s < numStages; s++)
 {
  eqFreq[s] = 1000.0;
  eqGain[s] = 0.0;
  eqQ[s]    = SQRT2_INV;
  eqDesigners[s].setMode(BiquadDesigner::RBJ_PEAK);
  eqDesigners[s].setSampleRate(sampleRate);
  eqDesigners[s].setFreq(1000.0);
  eqDesigners[s].setGain(0.0);
  eqDesigners[s].setQ(SQRT2_INV);
 }

}

MagicCarpetEqualizer::~MagicCarpetEqualizer()
{

}

//----------------------------------------------------------------------------
// parameter settings:

void MagicCarpetEqualizer::setSampleRate(double newSampleRate)
{
 if( newSampleRate > 0.01 )
  sampleRate = newSampleRate;

 intA s; // for indexing the channel and the filter-stage

 // update the sample-rate in the embdedded "BiquadDesigner" objects:
 for(s = 0; s < numStages; s++)
 {
  eqDesigners[s].setSampleRate(sampleRate);

  // let the biquad-designer with indices calculate biquad-coefficients
  // and store them at the s-th position in the respective 
  // coefficient-arrays:
  eqDesigners[s].getCoeffs(&(b0[s]), &(b1[s]), &(b2[s]), &(a1[s]), &(a2[s]));
 }

 // pass the new biquad-coefficients to the "BiquadCasccade"-objects:
 eqChainL.setCoeffs(b0, b1, b2, a1, a2);
 eqChainR.setCoeffs(b0, b1, b2, a1, a2);
}


void MagicCarpetEqualizer::setEqFreq(double newFreq, int stage)
{
 static intA s;
 s = stage;

 // make sure that we don't try to access an invalid adress:
 if( stage > (numStages-1) )
  return;

 if( newFreq >= 20.0 && newFreq <= 20000.0 )
 {
  eqFreq[s] = newFreq;
  eqDesigners[s].setFreq(newFreq);

  // let the biquad-designer with indices (s,0) calculate biquad-coefficients
  // and store them at the s-th position in the respective 
  // coefficient-arrays for the left channel:
  eqDesigners[s].getCoeffs(&(b0[s]), &(b1[s]), &(b2[s]), &(a1[s]), &(a2[s]));

  // pass the new coefficient-arrays to the left biquad cascade:
  eqChainL.setCoeffs(b0, b1, b2, a1, a2);
  eqChainR.setCoeffs(b0, b1, b2, a1, a2);
 }
}

void MagicCarpetEqualizer::setEqGain(double newGain, int stage)
{
 static intA s;
 s = stage;

 // make sure that we don't try to access an invalid adress:
 if( stage > (numStages-1) )
  return;

 eqGain[s] = newGain;
 eqDesigners[s].setGain(newGain);

 eqDesigners[s].getCoeffs(&(b0[s]), &(b1[s]), &(b2[s]), &(a1[s]), &(a2[s]));
 eqChainL.setCoeffs(b0, b1, b2, a1, a2);
 eqChainR.setCoeffs(b0, b1, b2, a1, a2);

 // set the "isOff"-flag to true when all gains are 0 dB (with 0.05 dB margin)
 /*
 isOff = true;
 if( abs(eqGain[0]) < 0.05 ) 
  isOff = false;
 if( abs(eqGain[1]) < 0.05 ) 
  isOff = false;
 if( abs(eqGain[2]) < 0.05 ) 
  isOff = false;
  */
}

void MagicCarpetEqualizer::setEqQ(double newQ, int stage)
{
 static intA s;
 s = stage;

 // make sure that we don't try to access an invalid adress:
 if( stage > (numStages-1) )
  return;

 if( newQ > 0.00001 )
 {
  eqQ[s] = newQ;
  eqDesigners[s].setQ(newQ);

  // let the biquad-designer with indices (s,0) calculate biquad-coefficients
  // and store them at the s-th position in the respective 
  // coefficient-arrays for the left channel:
  eqDesigners[s].getCoeffs(&(b0[s]), &(b1[s]), &(b2[s]), &(a1[s]), &(a2[s]));

  // pass the new coefficient-arrays to the left biquad cascade:
  eqChainL.setCoeffs(b0, b1, b2, a1, a2);
  eqChainR.setCoeffs(b0, b1, b2, a1, a2);
 }
}

//----------------------------------------------------------------------------
// others:

void MagicCarpetEqualizer::resetEqCoeffs()
{
 static intA k;
 for(k=0; k<numStages; k++)
 {
  b0[k] = 1.0;
  b1[k] = 0.0;
  b2[k] = 0.0;
  a1[k] = 0.0;
  a2[k] = 0.0;
 }
}
