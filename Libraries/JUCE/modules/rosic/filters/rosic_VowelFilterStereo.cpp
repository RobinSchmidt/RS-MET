#include "rosic_VowelFilterStereo.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

VowelFilterStereo::VowelFilterStereo()
{
  vowel         = 0.4;   // vowel 'a'
  amount        = 1.0;
  shiftFactor   = 1.0;
  sampleRate    = 44100.0;
  sampleRateRec = 1.0 / sampleRate;

  // initialize the vowel-data:

  gains[0]           =    1.0;          // u, global gain

  vowelData[0][0][0] =  251.4;          // u, 1st formant, frequency
  vowelData[0][0][1] =    1.65;         // u, 1st formant, bandwidth
  vowelData[0][0][2] =   16.2;          // u, 1st formant, gain

  vowelData[0][1][0] =  603.6;          // u, 2nd formant, frequency
  vowelData[0][1][1] =    0.18;         // u, 2nd formant, bandwidth
  vowelData[0][1][2] =   11.6;          // u, 2nd formant, gain

  gains[1]           =    1.0;          // o, global gain

  vowelData[1][0][0] =  329.6;          // o, 1st formant, frequency
  vowelData[1][0][1] =    0.76;         // o, 1st formant, bandwidth
  vowelData[1][0][2] =   19.2;          // o, 1st formant, gain

  vowelData[1][1][0] =  588.9;          // o, 2nd formant, frequency
  vowelData[1][1][1] =    0.16;         // o, 2nd formant, bandwidth
  vowelData[1][1][2] =   12.0;          // o, 2nd formant, gain

  gains[2]           =    1.0;          // a, global gain

  vowelData[2][0][0] =  730.1;          // a, 1st formant, frequency
  vowelData[2][0][1] =    0.59;         // a, 1st formant, bandwidth
  vowelData[2][0][2] =   17.9;          // a, 1st formant, gain

  vowelData[2][1][0] = 1000.1;          // a, 2nd formant, frequency
  vowelData[2][1][1] =    0.32;         // a, 2nd formant, bandwidth
  vowelData[2][1][2] =   15.0;          // a, 2nd formant, gain

  gains[3]           =    1.0;          // e, global gain

  vowelData[3][0][0] =  302.5;          // e, 1st formant, frequency
  vowelData[3][0][1] =    1.30;         // e, 1st formant, bandwidth
  vowelData[3][0][2] =   24.1;          // e, 1st formant, gain

  vowelData[3][1][0] = 2226.6;          // e, 2nd formant, frequency
  vowelData[3][1][1] =    0.19;         // e, 2nd formant, bandwidth
  vowelData[3][1][2] =   25.9;          // e, 2nd formant, gain

  gains[4]           =    1.0;          // i, global gain

  vowelData[4][0][0] =  167.6;          // i, 1st formant, frequency
  vowelData[4][0][1] =    2.24;         // i, 1st formant, bandwidth
  vowelData[4][0][2] =   20.0;          // i, 1st formant, gain (measured: 27.0)

  vowelData[4][1][0] = 2210.0;          // i, 2nd formant, frequency
  vowelData[4][1][1] =    0.12;         // i, 2nd formant, bandwidth
  vowelData[4][1][2] =   17.5;          // i, 2nd formant, gain (measured: 17.5)

  gains[5]           =    1.0;          // u, global gain

  vowelData[5][0][0] =  251.4;          // u, 1st formant, frequency
  vowelData[5][0][1] =    1.65;         // u, 1st formant, bandwidth
  vowelData[5][0][2] =   16.2;          // u, 1st formant, gain

  vowelData[5][1][0] =  603.6;          // u, 2nd formant, frequency
  vowelData[5][1][1] =    0.18;         // u, 2nd formant, bandwidth
  vowelData[5][1][2] =   11.6;          // u, 2nd formant, gain




  resetBuffers();                // reset the filters memory buffers (to 0)
  updateFilterCoefficients();    // calculate coefficients from default specification
}

VowelFilterStereo::~VowelFilterStereo()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void VowelFilterStereo::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;

  sampleRateRec = 1.0 / sampleRate;
  updateFilterCoefficients();
}

void VowelFilterStereo::setVowel(double newVowel)
{
  if( newVowel >= 0.0 && newVowel <= 1.0 )
    vowel = newVowel;
}

void VowelFilterStereo::setAmount(double newAmount)
{
  if( newAmount >= 0.0 )
    amount = newAmount;
}

void VowelFilterStereo::setShift(double newShift)
{
  shiftFactor = pitchOffsetToFreqFactor(newShift);
}

//-------------------------------------------------------------------------------------------------
// inquiry:

double VowelFilterStereo::getVowel()
{
  return vowel;
}

double VowelFilterStereo::getAmount()
{
  return amount;
}

//-------------------------------------------------------------------------------------------------
// others:

void VowelFilterStereo::resetBuffers()
{
  xL_s1_d1 = 0.0;
  xR_s1_d1 = 0.0;
  xL_s1_d2 = 0.0;  
  xR_s1_d2 = 0.0;   
  yL_s1_d1 = 0.0;  
  yR_s1_d1 = 0.0;   
  yL_s1_d2 = 0.0;    
  yR_s1_d2 = 0.0;

  xL_s2_d1 = 0.0;
  xR_s2_d1 = 0.0;
  xL_s2_d2 = 0.0;  
  xR_s2_d2 = 0.0;   
  yL_s2_d1 = 0.0;  
  yR_s2_d1 = 0.0;   
  yL_s2_d2 = 0.0;    
  yR_s2_d2 = 0.0;
}



