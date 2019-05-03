#include "BiquadCascadeStereoSSE2.h"

//----------------------------------------------------------------------------
//construction/destruction:

BiquadCascadeStereoSSE2::BiquadCascadeStereoSSE2()
{
 initBiquadCoeffs();        // initializes the filter-coefficents
 reset();                   // resets the filters memory buffers (to 0)
 numStages = maxNumStages;  // 
 gain      = 1.0;
}

BiquadCascadeStereoSSE2::~BiquadCascadeStereoSSE2()
{}

//----------------------------------------------------------------------------
//parameter settings:

void BiquadCascadeStereoSSE2::setSampleRate(double newSampleRate)
{
 if( newSampleRate > 0.0 )
  sampleRate = newSampleRate;

 sampleRateRec = 1.0 / sampleRate;
}

void BiquadCascadeStereoSSE2::setNumStages(int newNumStages)
{
 if( (newNumStages >= 0 ) && (newNumStages <= maxNumStages) )
  numStages = newNumStages;
 reset();
}

//----------------------------------------------------------------------------
//others:

void BiquadCascadeStereoSSE2::initBiquadCoeffs()
{
 int i;
 for (i=0; i<maxNumStages; i++)
 {
  b0[i] = 1.0;
  b1[i] = 0.0;
  b2[i] = 0.0;
  a0[i] = 1.0;
  a1[i] = 0.0;
  a2[i] = 0.0;
 }
}

void BiquadCascadeStereoSSE2::reset()
{
 int i;
 for (i=0; i<maxNumStages; i++)
 {
  x1[i] = _mm_set1_pd(0.0);
  x2[i] = _mm_set1_pd(0.0);
  x1[i] = _mm_set1_pd(0.0);
  y2[i] = _mm_set1_pd(0.0);
  y1[i] = _mm_set1_pd(0.0);
 }
}

void BiquadCascadeStereoSSE2::getMagnitudeResponse(double *frequencies, 
                                         double *magnitudes,
                                         int     numBins, 
                                         bool    inDecibels,
                                         bool    accumulate)
{
 static doubleA num,  // numerator of the squared magnitude of one stage
                den,  // denominator of the squared magnitude of one stage
                accu; // accumulator which multiplicatively accumulates
                      // the squared magnitude of the individual stages

 static doubleA omega;     // current normalized radian frequency
 static doubleA cosOmega;  // cosine of the current omega
 static doubleA cos2Omega; // cosine of twice the current omega

 static intA k;    // index for the current frequency bin    
 static intA s;    // index for the current biquad-stage

 for(k=0; k<numBins; k++)
 {
  // convert frequency in Hz to normalized radian frequency omega:
  omega = (2.0*PI)*(frequencies[k]*sampleRateRec);

  // calculate cos(omega) and cos(2*omega) because these are needed in the 
  // numerator and denominator as well:
  cosOmega  = cos(omega);
  cos2Omega = cos(2*omega);

  // accumulate the squared magnitude responses of the individual 
  // biquad-stages:
  accu = 1.0;
  for(s=0; s<numStages; s++)
  {
   // calculate the numerator of the squared magnitude of stage s:
   num = b0[s]*b0[s] + b1[s]*b1[s] + b2[s]*b2[s]
         + 2*cosOmega*(b0[s]*b1[s] + b1[s]*b2[s])
         + 2*cos2Omega*b0[s]*b2[s];

   // calculate the denominator of the squared magnitude of stage s:
   den = a0[s]*a0[s] + a1[s]*a1[s] + a2[s]*a2[s]
         + 2*cosOmega*(a0[s]*a1[s] + a1[s]*a2[s])
         + 2*cos2Omega*a0[s]*a2[s];

   // multiply the accumulator with the squared magnitude of stage s:
   accu *= (num/den);

  } // end of "for(s=0; s<numStages; s++)"

  // take the square root of the accumulated squared magnitude response - this
  // is the desired magnitude of the biquad cascade at frequencies[k]:
  accu = sqrt(accu);
  
  // store the calculated value in the array "magnitudes", taking into acount,
  // if the value should be in dB or not and if it should accumulate to what's 
  // already there or not:
  if( !inDecibels && !accumulate )
   magnitudes[k]  = (float) accu;
  else if( !inDecibels && accumulate )
   magnitudes[k] *= (float) accu;
  else if( inDecibels && !accumulate )
   magnitudes[k]  = (float) 20*log10(accu);
  else if( inDecibels && accumulate )
   magnitudes[k] += (float) 20*log10(accu);

 } // end of " for(k=0; k<numBins; k++)"

 /*
 // convert to decibels if this is desired:
 if( inDecibels == true)
  for(k=0; k<numBins; k++)
   magnitudes[k] = 20 * log10(magnitudes[k]);
   //magnitudes[k] = 16.0; // test
 */

}

void BiquadCascadeStereoSSE2::getMagnitudeResponse(float *frequencies, 
                                         float *magnitudes,
                                         int    numBins, 
                                         bool   inDecibels,
                                         bool   accumulate)
{
 static doubleA num,  // numerator of the squared magnitude of one stage
               den,  // denominator of the squared magnitude of one stage
               accu; // accumulator which multiplicatively accumulates
                     // the squared magnitude of the individual stages

 static doubleA omega;     // current normalized radian frequency
 static doubleA cosOmega;  // cosine of the current omega
 static doubleA cos2Omega; // cosine of twice the current omega

 static intA k;    // index for the current frequency bin    
 static intA s;    // index for the current biquad-stage

 for(k=0; k<numBins; k++)
 {
  // convert frequency in Hz to normalized radian frequency omega:
  omega = 2*PI*frequencies[k]/sampleRate;

  // calculate cos(omega) and cos(2*omega) because these are needed in the 
  // numerator and denominator as well:
  cosOmega  = cos(omega);
  cos2Omega = cos(2*omega);

  // accumulate the squared magnitude responses of the individual 
  // biquad-stages:
  accu = 1.0;
  for(s=0; s<numStages; s++)
  {
   // calculate the numerator of the squared magnitude of stage s:
   num = b0[s]*b0[s] + b1[s]*b1[s] + b2[s]*b2[s]
         + 2*cosOmega*(b0[s]*b1[s] + b1[s]*b2[s])
         + 2*cos2Omega*b0[s]*b2[s];

   // calculate the denominator of the squared magnitude of stage s:
   den = a0[s]*a0[s] + a1[s]*a1[s] + a2[s]*a2[s]
         + 2*cosOmega*(a0[s]*a1[s] + a1[s]*a2[s])
         + 2*cos2Omega*a0[s]*a2[s];

   // multiply the accumulator with the squared magnitude of stage s:
   accu *= (num/den);

  } // end of "for(s=0; s<numStages; s++)"

  // take the square root of the accumulated squared magnitude response - this
  // is the desired magnitude of the biquad cascade at frequencies[k]:
  accu = sqrt(accu);
  
  // store the calculated value in the array "magnitudes", taking into acount,
  // if the value should be in dB or not and if it should accumulate to what's 
  // already there or not:
  if( !inDecibels && !accumulate )
   magnitudes[k]  = (float) accu;
  else if( !inDecibels && accumulate )
   magnitudes[k] *= (float) accu;
  else if( inDecibels && !accumulate )
   magnitudes[k]  = (float) (20*log10(accu));
  else if( inDecibels && accumulate )
   magnitudes[k] += (float) (20*log10(accu));

 } // end of " for(k=0; k<numBins; k++)"

 /*
 // convert to decibels if this is desired:
 if( inDecibels == true)
  for(k=0; k<numBins; k++)
   magnitudes[k] = 20 * log10(magnitudes[k]);
   //magnitudes[k] = 16.0; // test
 */

}



