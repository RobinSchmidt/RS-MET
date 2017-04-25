#include "IirDesigner.h"

IirDesigner::IirDesigner()
{
 //init member variables:
 sampleRate = 44100.0;
 mode       = LOWPASS;     // LPF by default
 slope      = 4;           // 24 dB/Oct slope filter by default 
 method     = BUTTERWORTH; // Butterworth approximation by default
 freq1      = 4410.0;      // cutoff for LPF/HPF, lower cutoff for BPF/BRF
 freq2      = 8820.0;      // upper cutoff for BPF/BRF
 gain       = 2.0;

 numPrototypePoles        = 4;
 numPrototypeZeros        = 0;
 numPrototypeBiquads      = 2;
 firstOrderStageIsPresent = false;

 passbandRipple   = 1.0;   // 1 dB passband ripple by default (unused in
                           // the default Butterworth-filter)
 stopbandRipple   = 1.0;   // 1 dB stopband ripple by default (unused in
                           // the default Butterworth-filter)
}

IirDesigner::~IirDesigner()
{

}

//----------------------------------------------------------------------------
//parameter settings:

void IirDesigner::setSampleRate(double  newSampleRate)
{
 if( newSampleRate>0 )
  sampleRate = newSampleRate;
}

void IirDesigner::setMode(int newMode)
{
 if( (newMode>=0) && (newMode<=7) )
  mode = newMode;
}

void IirDesigner::setSlope(int newSlope)
{
 //if( (newSlope>=1) && (newSlope<=(maxOrder/2)) )
 slope = newSlope;

 if( slope > 0 && isOdd(slope) )
 {
  numPrototypeBiquads      = (slope-1)/2;
  firstOrderStageIsPresent = true;
 }
 else if( slope > 0 && isEven(slope) )
 {
  numPrototypeBiquads      = slope/2;
  firstOrderStageIsPresent = false;
 }
 else
 {
  numPrototypeBiquads      = 0;
  firstOrderStageIsPresent = false;
 }
}

void IirDesigner::setFreq1(double newFreq1)
{
 if( (newFreq1>=0) && (newFreq1<=(0.5*sampleRate)) )
  freq1 = newFreq1;
}
void IirDesigner::setFreq2(double newFreq2)
{
 if( (newFreq2>=0) && (newFreq2<=(0.5*sampleRate)) )
  freq2 = newFreq2;
}
void IirDesigner::setGain(double newGain)
{
 if( newGain>=0.0 )
  gain = newGain;
}


//----------------------------------------------------------------------------
// coefficient calculation:

void IirDesigner::getDirectFormCoeffs(double *FeedforwardCoeffs, 
                                      double *FeedbackCoeffs)
{
 // do nothing in the case, that the filters slope is zero:
 if( slope <= 0)
  return;

 prewarpFreqs();

 // use a switch statement here later to switch between the
 // approximation methods:
 //calcButterUnitLpfAna();
 //assignPrototypePolesAndZeros();
 calcPrototypePolesAndZeros();

 frequencyTransform();

 bilinearTransform();

 calcGainFactor();

 calcDirectFormCoeffs();

 //copy the calculated coefficients into the dedicated buffers:
 for(long k=0; k<(numPolesDig+1); k++)
 {
  FeedforwardCoeffs[k] = ffCoeffs[k];
  FeedbackCoeffs[k]    = fbCoeffs[k];
 }

}

void IirDesigner::getBiquadCascadeCoeffs(double *b0, double *b1, double *b2, 
                                                     double *a1, double *a2)
{
 // do nothing in the case, that the filters slope is zero:
 if( slope <= 0)
  return;

 // clear the internal arrays:
 for(int i=0; i<maxOrder; i++)
 {
  prototypePoles[i] = Complex(0.0, 0.0);
  prototypeZeros[i] = Complex(0.0, 0.0);

  polesAna[i] = Complex(0.0, 0.0);
  zerosAna[i] = Complex(0.0, 0.0);

  polesDig[i] = Complex(0.0, 0.0);
  zerosDig[i] = Complex(0.0, 0.0);
 }

 // do all the stuff which calculates the postions of the poles and zeros in
 // the z-plane - the member-variables "polesDig" and "zerosDig" will contain
 // the correct values after that:
 prewarpFreqs();
 //assignPrototypePolesAndZeros();
 calcPrototypePolesAndZeros();
 frequencyTransform();
 bilinearTransform();

 // combine every two successive poles and zeros to one biquad-stage (the 
 // poles and zeros are in an order such that two successive array-entries
 // are complex-conjugate to each other):
 intA k;
 for(k=0; k<numPolesDig; k+=2)
 {
  b0[k/2] = 1.0;
  b1[k/2] = -(zerosDig[k] + zerosDig[k+1]).re;
  b2[k/2] =  (zerosDig[k] * zerosDig[k+1]).re;

  a1[k/2] = -(polesDig[k] + polesDig[k+1]).re;
  a2[k/2] =  (polesDig[k] * polesDig[k+1]).re;
 }

 // if there is a single left over pole (and zero), it will be realized by 
 // the last biquad-stage:
 if( isOdd((int)numPolesDig) )
 {

  k = numPolesDig-1;
  b0[k/2] = 1.0;
  b1[k/2] = -zerosDig[k].re;
  b2[k/2] = 0.0;

  a1[k/2] = -polesDig[k].re;
  a2[k/2] = 0.0;
 }

 // normalize the biquad stages such that each has unit magnitude at the 
 // passband-center frequency (use dc for bandreject filters):
 static doubleA g, normalizer, omegaDigC;
 static intA numStages;
 if( isEven((int)numPolesDig) )
  numStages = numPolesDig/2;
 else
  numStages = (numPolesDig+1)/2;

 calcGainFactor(); 
  // the member-variable "gainFactor" now contains a value, to be multiplied
  // with the filters input- or output-signal to achieve unity gain at some 
  // frequency at which the filter is assumed to be neutral

 // we realize this gain factor at the input-stage:
 b0[0] *= gainFactor;
 b1[0] *= gainFactor;
 b2[0] *= gainFactor;

 /*
 switch( mode )
 {
 case LOWPASS:
  {
   for(k=0; k<numStages; k++)
   {
    // get gain at DC:
    g = getBiquadMagnitudeAt(b0[k], b1[k], b2[k], a1[k], a2[k], 0.0);

    // divide the b-coefficients by that factor:
    normalizer = 1.0/g;
    b0[k]      = normalizer * b0[k];
    b1[k]      = normalizer * b1[k];
    b2[k]      = normalizer * b2[k];
   }
  }
  break;
 case HIGHPASS:
  {
   for(k=0; k<numStages; k++)
   {
    // get gain at Nyquist-frequency (PI):
    g = getBiquadMagnitudeAt(b0[k], b1[k], b2[k], a1[k], a2[k], PI);

    // divide the b-coefficients by that factor:
    normalizer = 1.0/g;
    b0[k]      = normalizer * b0[k];
    b1[k]      = normalizer * b1[k];
    b2[k]      = normalizer * b2[k];
   }
  }
  break;
 case BANDPASS:
  {  
   for(k=0; k<numStages; k++)
   {
    // get gain at the normalized radian center-frequency:
    omegaDigC = (2*PI/sampleRate) * sqrt(freq1*freq2);
    g = getBiquadMagnitudeAt(b0[k], b1[k], b2[k], a1[k], a2[k], omegaDigC);

    // divide the b-coefficients by that factor:
    normalizer = 1.0/g;
    b0[k]      = normalizer * b0[k];
    b1[k]      = normalizer * b1[k];
    b2[k]      = normalizer * b2[k];
   }
  }
  break;
 case BANDREJECT:
  {
   for(k=0; k<numStages; k++)
   {
    // get gain at DC:
    g = getBiquadMagnitudeAt(b0[k], b1[k], b2[k], a1[k], a2[k], 0.0);

    // divide the b-coefficients by that factor:
    normalizer = 1.0/g;
    b0[k]      = normalizer * b0[k];
    b1[k]      = normalizer * b1[k];
    b2[k]      = normalizer * b2[k];
   }
  }
  break;
 case LOW_SHELV:
  {
   for(k=0; k<numStages; k++)
   {
    // get gain at Nyquist-frequency (PI):
    g = getBiquadMagnitudeAt(b0[k], b1[k], b2[k], a1[k], a2[k], PI);

    // divide the b-coefficients by that factor:
    normalizer = 1.0/g;
    b0[k]      = normalizer * b0[k];
    b1[k]      = normalizer * b1[k];
    b2[k]      = normalizer * b2[k];
   }
  }
  break;
 case HIGH_SHELV:
  {
   for(k=0; k<numStages; k++)
   {
    // get gain at DC:
    g = getBiquadMagnitudeAt(b0[k], b1[k], b2[k], a1[k], a2[k], 0.0);

    // divide the b-coefficients by that factor:
    normalizer = 1.0/g;
    b0[k]      = normalizer * b0[k];
    b1[k]      = normalizer * b1[k];
    b2[k]      = normalizer * b2[k];
   }
  }
  break;
 case PEAK:
  {
   for(k=0; k<numStages; k++)
   {
    // get gain at DC:
    g = getBiquadMagnitudeAt(b0[k], b1[k], b2[k], a1[k], a2[k], 0.0);

    // the last biquad-stage in an eq based on an odd-order prototype does not
    // need to have unity gain at dc and does not need any normalizing
    if( k==numStages && isOdd(numPrototypePoles) )
     g = 1.0;

    // divide the b-coefficients by that factor:
    normalizer = 1.0/g;
    b0[k]      = normalizer * b0[k];
    b1[k]      = normalizer * b1[k];
    b2[k]      = normalizer * b2[k];
   }
  }
  break;
 }  // end of "switch( mode )
 */
}


/*
void IirDesigner::getBiquadCascadeCoeffs(double *b0, double *b1, double *b2, 
                                         double *a1, double *a2)
{

}
*/

//----------------------------------------------------------------------------
//internal calculations:
void IirDesigner::prewarpFreqs()
{
 // convert the cutoff-frequency(ies) into a normalized cutoff-frequency(ies)
 // in radians/sec (a frequency between 0 and pi):
 omegaDig1   = (2*PI*freq1)/sampleRate;
 omegaDig2   = (2*PI*freq2)/sampleRate;

 // prewarp the desired digital cutoff radian frequency to the desired analog
 // cutoff radian frequency:
 omegaAna1   = 2*sampleRate * tan(0.5*omegaDig1);
 omegaAna2   = 2*sampleRate * tan(0.5*omegaDig2);

 //for bandpass and bandreject:
 omegaAnaC = sqrt(omegaAna1*omegaAna2);
 bwAna     = omegaAna2-omegaAna1;
}

void IirDesigner::calcPrototypePolesAndZeros()
{
 switch(method)
 {
  case BUTTERWORTH:
   calcButterworthPrototypePolesAndZeros();
  break;

  // ... more to come

  default:
   calcButterworthPrototypePolesAndZeros();
 }; // end of switch(method)
}

void IirDesigner::calcButterworthPrototypePolesAndZeros()
{
 static intA  k;
 static Complex p, z;  // for the current pole and zero

 // variables for some intermediate results
 static doubleA g_b;      // gain at the bandwidth ( = sqrt(gain) )
 static doubleA epsilon;
 static doubleA g, beta, phi, s, c; 

 if( mode==LOWPASS || mode==HIGHPASS || mode==BANDPASS || mode==BANDREJECT )
 // generate a lowpass prototype:
 {
  numPrototypePoles = slope;
  numPrototypeZeros = 0;   // no zeros in the numerator of H(s) -> all zeros     
                           // of H(s) at infinity
  k = 0;
  while( k < numPrototypeBiquads )
  {
   p.re = -sin( ((2*(k+1)-1)*PI)/(2.0*numPrototypePoles) );
   p.im =  cos( ((2*(k+1)-1)*PI)/(2.0*numPrototypePoles) );

   prototypePoles[2*k]   = p;
   prototypePoles[2*k+1] = p.getConj();

   k += 1;
  }
  if( firstOrderStageIsPresent )
   prototypePoles[numPrototypePoles-1] = Complex(-1.0, 0.0);
 }
 else if( mode==LOW_SHELV || mode==HIGH_SHELV || mode==PEAK )
 // generate a low-shelv prototype:
 {
  numPrototypePoles = slope;
  numPrototypeZeros = slope;   
  if( fabs(gain-1.0) < 0.0001 )
  {
   // gain is too close to 1.0 (all poles and zeros coincide, when gain == 1.0)
   // make a "bypass-prototype":
   numPrototypePoles = 0;
   numPrototypeZeros = 0;
  }
  else // gain has reasonable value
  {
   g_b     = sqrt(gain);
   epsilon = sqrt( (gain*gain-g_b*g_b)/(g_b*g_b-1.0) ); // Eq. 12 with G_0=1.0
   g       = pow(gain, 1.0/numPrototypePoles );     // Eq. 94
   beta    = pow(epsilon, -1.0/numPrototypePoles);  // Eq. 94 with OmegaB=1.0
   k       = 0;
   while( k < numPrototypeBiquads )
   {
    phi  = ((2*(k+1)-1)*PI) / (2.0*numPrototypePoles);  // Eq.95
    s    = sin(phi);  // Eq. 95
    c    = cos(phi);  // Eq. 95

    p.re = -s*beta;   // Eq. 93
    p.im =  c*beta;   // Eq. 93
    z.re = -s*g*beta; // Eq. 93
    z.im =  c*g*beta; // Eq. 93

    prototypePoles[2*k]   = p;
    prototypePoles[2*k+1] = p.getConj();
    prototypeZeros[2*k]   = z;
    prototypeZeros[2*k+1] = z.getConj();

    k += 1;
   }
   if( firstOrderStageIsPresent )
   {
    z = Complex(-g*beta, 0.0);
    p = Complex(-beta, 0.0);
    prototypePoles[numPrototypePoles-1] = p;
    prototypeZeros[numPrototypePoles-1] = z;
   }
  } // end of else
 } // end of else if( mode==LOW_SHELV || mode==HIGH_SHELV || mode==PEAK )
 else
 {
  // mode is none of the defined ones - make a "bypass-prototype":
  numPrototypePoles = 0;
  numPrototypeZeros = 0;
 }

}


void IirDesigner::frequencyTransform()
{
 intA k;           // for indexing the poles and zeros
 Complex temp1, temp2; // for some intermediate values in the BPF and BRF case

 switch(mode)
 {

  case LOWPASS: //lowpass->lowpass transform
  {
   numPolesAna = numPrototypePoles;
   numZerosAna = numPrototypeZeros;
   for(k=0; k<numPolesAna; k++)
    polesAna[k] = prototypePoles[k] * omegaAna1;
   for(k=0; k<numZerosAna; k++)
    zerosAna[k] = prototypeZeros[k] * omegaAna1;
  }
  break;

  case HIGHPASS: //lowpass->highpass transform
  {
   numPolesAna = numPrototypePoles;
   numZerosAna = numPrototypeZeros;
   for(k=0; k<numPolesAna; k++)
    polesAna[k] = Complex(omegaAna1)/prototypePoles[k];
   for(k=0; k<numZerosAna; k++)
    zerosAna[k] = Complex(omegaAna1)/prototypeZeros[k];
   // if there are no zeros in the numerator (indicated by numZerosAna==0),
   // then all zeros of H(s) are at s=inf - they map to s=0:
   if(numZerosAna==0)
   {
    numZerosAna = numPolesAna;
    for(k=0; k<numZerosAna; k++)
     zerosAna[k] = 0;
   }
  }
  break;

  case BANDPASS: //lowpass->bandpass transform
  {
   numPolesAna = 2*numPrototypePoles; //twice as much poles as in the LPF
   numZerosAna = 2*numPrototypeZeros; 
   for(k=0; k<numPrototypePoles; k++)
   {
    temp1 = prototypePoles[k]*bwAna*0.5;
    temp2 = -prototypePoles[k]*bwAna;
    temp2 = temp2*temp2*0.25;
    temp2 = Complex::sqrt(temp2 - Complex(omegaAnaC*omegaAnaC));
    polesAna[k]               = temp1 + temp2;
    polesAna[numPolesAna-1-k] = temp1 - temp2;
   }
   for(k=0; k<numPrototypeZeros; k++)
   {
    temp1 = prototypeZeros[k]*bwAna*0.5;
    temp2 = -prototypeZeros[k]*bwAna;
    temp2 = temp2*temp2*0.25;
    temp2 = Complex::sqrt(temp2 - Complex(omegaAnaC*omegaAnaC));
    zerosAna[k]               = temp1 + temp2;
    zerosAna[numZerosAna-1-k] = temp1 - temp2;
   }
   // if there are no zeros in the numerator (indicated by numZerosAna==0),
   // then we will get numPolesAna/2 zeros at s=inf and numPolesAna/2 zeros
   // at s=0:
   if(numPrototypeZeros==0)
   {
    numZerosAna = numPolesAna/2;
    for(k=0; k<numZerosAna; k++)
     zerosAna[k] = 0;
   }
  }
  break;

  case BANDREJECT: //lowpass->bandreject transform
  {
   numPolesAna = 2*numPrototypePoles; //twice as much poles as in the LPF
   numZerosAna = 2*numPrototypeZeros; 
   for(k=0; k<numPrototypePoles; k++)
   {
    temp1 = Complex(bwAna*0.5)/prototypePoles[k];
    temp2 = Complex(-bwAna)/prototypePoles[k];
    temp2 = temp2*temp2*0.25;
    temp2 = Complex::sqrt(temp2 - Complex(omegaAnaC*omegaAnaC));
    polesAna[k]               = temp1 + temp2;
    polesAna[numPolesAna-1-k] = temp1 - temp2;
   }
   for(k=0; k<numPrototypeZeros; k++)
   {
    temp1 = Complex(bwAna*0.5)/prototypeZeros[k];
    temp2 = Complex(-bwAna)/prototypeZeros[k];
    temp2 = temp2*temp2*0.25;
    temp2 = Complex::sqrt(temp2 - Complex(omegaAnaC*omegaAnaC));
    zerosAna[k]               = temp1 + temp2;
    zerosAna[numZerosAna-1-k] = temp1 - temp2;
   }
   // if there are no zeros in the numerator (indicated by numZerosAna==0),
   // then we will get numPolesAna/2 zeros at s=j*omegaAnaC and 
   // numPolesAna/2 zeros at s=-j*omegaAnaC:
   if(numPrototypeZeros==0)
   {
    numZerosAna = numPolesAna;
    for(k=0; k<(numZerosAna/2); k++)
    {
     temp1    = Complex(0,0);
     temp1.im = omegaAnaC;
     zerosAna[2*k]   =  temp1;
     zerosAna[2*k+1] = -temp1;
     //zerosAna[k]               =  temp1;
     //zerosAna[numZerosAna-1-k] = -temp1;
    }
   }
  }
  break;

  case LOW_SHELV: // low-shelv -> low-shelv transform
  {
   numPolesAna = numPrototypePoles;
   numZerosAna = numPrototypeZeros;
   for(k=0; k<numPolesAna; k++)
    polesAna[k] = prototypePoles[k] * omegaAna1;
   for(k=0; k<numZerosAna; k++)
    zerosAna[k] = prototypeZeros[k] * omegaAna1;
  }
  break;

  case HIGH_SHELV: // low-shelv -> high-shelv transform
  {
   numPolesAna = numPrototypePoles;
   numZerosAna = numPrototypeZeros;
   for(k=0; k<numPolesAna; k++)
    polesAna[k] = prototypeZeros[k] * omegaAna1;
   for(k=0; k<numZerosAna; k++)
    zerosAna[k] = prototypePoles[k] * omegaAna1;
  }
  break;

  case PEAK: // low-shelv -> peaking transform
  {
   numPolesAna = 2*numPrototypePoles; //twice as much poles as in the LS
   numZerosAna = 2*numPrototypeZeros; 
   for(k=0; k<numPrototypePoles; k++)
   {
    temp1 = prototypePoles[k]*bwAna*0.5;
    temp2 = -prototypePoles[k]*bwAna;
    temp2 = temp2*temp2*0.25;
    temp2 = Complex::sqrt(temp2 - Complex(omegaAnaC*omegaAnaC));
    polesAna[k]               = temp1 + temp2;
    polesAna[numPolesAna-1-k] = temp1 - temp2;
   }
   for(k=0; k<numPrototypeZeros; k++)
   {
    temp1 = prototypeZeros[k]*bwAna*0.5;
    temp2 = -prototypeZeros[k]*bwAna;
    temp2 = temp2*temp2*0.25;
    temp2 = Complex::sqrt(temp2 - Complex(omegaAnaC*omegaAnaC));
    zerosAna[k]               = temp1 + temp2;
    zerosAna[numZerosAna-1-k] = temp1 - temp2;
   }
  }
  break;
 } //end of switch

}

void IirDesigner::bilinearTransform()
{
 long    k;
 Complex x;

 numPolesDig = numPolesAna;
 numZerosDig = numPolesAna; // numZerosAna may contain only half of the zeros 
                            // in the digital filter when there are some at 
                            // s=inf

 //init all zeros with z=-1 (this is where zeros at s=inf map to):
 for(k=0; k<numZerosDig; k++)
  zerosDig[k] = Complex(-1);

 //transform the poles from the s-plane to the z-plane:
 for(k=0; k<numPolesAna; k++)
 {
  x           = polesAna[k] * (0.5/sampleRate);
  polesDig[k] = (Complex(1.0)+x)/(Complex(1.0)-x); //Complex(1.0) creates the complex number 1.0 + 0.0j
  //polesDig[k].print();
 }
 //transform the zeros from the s-plane to the z-plane:
 for(k=0; k<numZerosAna; k++)
 {
  x           = zerosAna[k] * (0.5/sampleRate);
  zerosDig[k] = (Complex(1.0)+x)/(Complex(1.0)-x); //Complex(1.0) creates the complex number 1.0 + 0.0j
  //zerosDig[k].print();
 }

 //in some filters the numerator of H(s) has no zeros at all. in this case all zeros
 //of H(s) are at s=inf. these zeros map to z=-1 in the lowpass case and z=+1 in
 //the highpass case
 if(numZerosAna==0)
 {
  switch(mode)
  {
   case 1: //lowpass mode - all zeros are at -1 (and there are as many zeros as there are poles)
   {
    for(k=0; k<numZerosDig; k++)
    {
     zerosDig[k] = Complex(-1.0); //Complex(-1) creates the complex number -1.0 + 0.0j
     //zerosDig[k].print();
    }
   }
   break;

  } //end of switch(mode)
 } //end of if(numZeros==0)

}

void IirDesigner::calcGainFactor()
{
 Complex num(1,0); // accumulates the numerator of H(z)
 Complex den(1,0); // accumulates the denominator of H(z)
 Complex z;        // value of z at which H(z) is evaluated
 double  temp;

 // decide at which value of z the transfer function H(z) should
 // be evaluated:
 switch(mode)
 {
  case LOWPASS: // evaluate gain at dc (z = 1 + 0j)
   z = Complex(1,0);
  break;
  case HIGHPASS: // evaluate gain at half the sample rate (z = -1 + 0j)
   z = Complex(-1,0);
  break;
  case BANDPASS: // evaluate gain at the center frequency:
  {
   temp = (2*PI/sampleRate) * sqrt(freq1*freq2); 
   z    = Complex(0,0);
   z.im = temp;
   z    = Complex::exp(z);
  }
  break;
  case BANDREJECT: // evaluate gain at dc (z = 1 + 0j)
   z = Complex(1,0);
  break;
  case LOW_SHELV: // evaluate gain at nyquist-frequency (z = -1 + 0j)
   z = Complex(-1,0);
  break;
  case HIGH_SHELV: // evaluate gain at dc (z = 1 + 0j)
   z = Complex(1,0);
  break;
  case PEAK:       // evaluate gain at dc (z = 1 + 0j)
   z = Complex(1,0);
  break;
 }
 //z.print();

 //evaluate H(z):
 static intA k;
 for(k=0; k<numPolesDig; k++)
  den = den * ( z - polesDig[k] );
 for(k=0; k<numZerosDig; k++)
  num = num * ( z - zerosDig[k] );

 Complex H    = num/den;
 double  gain = H.getMagnitude();

 gainFactor = 1/gain;
}

double IirDesigner::getBiquadMagnitudeAt(double b0, double b1, double b2, 
                                        double a1, double a2, double omega)
{
 static doubleA num, den, mag, cosOmega, cos2Omega;

 cosOmega  = cos(omega);
 cos2Omega = cos(2*omega);

 num = b0*b0 + b1*b1 + b2*b2 + 2*cosOmega*(b0*b1 + b1*b2) + 2*cos2Omega*b0*b2;
 den = 1.0   + a1*a1 + a2*a2 + 2*cosOmega*(a1    + a1*a2) + 2*cos2Omega*a2;
 mag = sqrt(num/den);

 return mag;
}

void IirDesigner::calcDirectFormCoeffs()
{
 Complex *polyCoeffs = new Complex[(int)numPolesDig+1]; 
  // this buffer accumulates the final direct-form filter coefficients via a
  // recursively applied convolution of the individual linear factors of H(z).
  // Polynomials are multiplied by convolving their coefficients - here we
  // start with the constant polynomial p(x) = 1, multiply it with (x - x0)
  // via convolution of the coeffs where x0 is one of the poles/zeros, then
  // multiply the result with (x-x1), then the result with (x-x2), etc. ; the
  // index ranges from 0 to "slope" so that we can store slope+1 coeffs

 Complex *linearCoeffs = new Complex[2];     
  // this holds the polynomial coefficients of the current linear factor 
  // (x - x0) where x0 is the zero of the linear factor - the coefficients are
  // 1 for x^1-term and -zerosDig[k] or -polesDig[k] for the x^0-term

 linearCoeffs[0] = Complex(1.0, 0.0);        
  //this one is always the same (see above)

 Complex *buffer = new Complex[(int)numPolesDig+1];     
  // this is needed temporarily to store the result of one convolution (of the
  // current polyCoeffs array with one the linear factor) - then the array is
  // copied into polyCoeffs to apply the next convolution
                                        
 long k = 0;
 long i = 0;

 // initialize the polynomial coefficients:
 polyCoeffs[0] = Complex(1.0, 0.0);
 for(k=1; k<(numPolesDig+1); k++)
  polyCoeffs[k] = Complex(0.0 , 0.0);

 // recursively convolve the polynomial coefficients with all the linear
 // factors:
 for(k=0; k<numPolesDig; k++)
 {
  // choose the k-th linear factor and represent it by its coefficients (the
  // first coefficient linearCoeffs[0] corresponding to x^1 is always 1):
  linearCoeffs[1] = -polesDig[k];

  // convolve the polynomial coefficients currently stored in polyCoeffs with
  // these new linear factor coefficients and store the result in the buffer:
  convolve(polyCoeffs, k+1, linearCoeffs, 2, buffer);

  // copy the new polynomial coefficients which are in the buffer now, into
  // the polyCoeffs array for the next iteration:
  for(i=0; i<=(k+1); i++)   //i goes from 0 to slope-1+1 when k=slope-1
   polyCoeffs[i] = buffer[i];
 }

 // this has lead us to the final polynomial coefficients the imaginary part
 // of which should ideally cancel out to zero. the real part is used as the
 // filter feedback coefficients (as we have worked with the poles here)
 for(k=0; k<(numPolesDig+1); k++)
  fbCoeffs[k] = polyCoeffs[k].re;

 //do the same procedure for the zeros of H(z) - this time without comments:
 polyCoeffs[0] = Complex(1.0, 0.0);
 for(k=1; k<(numZerosDig+1); k++)
  polyCoeffs[k] = Complex(0.0 , 0.0);
 for(k=0; k<numZerosDig; k++)
 {
  linearCoeffs[1] = -zerosDig[k];
  convolve(polyCoeffs, k+1, linearCoeffs, 2, buffer);
  for(i=0; i<=(k+1); i++)   
   polyCoeffs[i] = buffer[i];
 }
 for(k=0; k<(numZerosDig+1); k++)
  ffCoeffs[k] = gainFactor * polyCoeffs[k].re; //the feedforward coeffs include the 
                                               //overall gain factor
}

void IirDesigner::calcBiquadCascadeCoeffs()
{

}

void IirDesigner::convolve(Complex *Seq1, int Length1,
                           Complex *Seq2, int Length2,
                           Complex *Result)
{
 int lengthR = Length1+Length2-1; //length of the resulting sequence

 //init the result-buffer with zeros:
 int n, k;
 for(n=0; n<lengthR; n++) 
  Result[n] = Complex(0.0, 0.0);

 for(n=0; n<lengthR; n++)
  for(k=0; k<Length2; k++)
   if((n-k)>=0)
    Result[n] = Result[n] + (Seq2[k] * Seq1[n-k]);
    //Seq2 corresponds to a impulse response and Seq1 to the input signal
}

void IirDesigner::calculateAndPrint()
{
 prewarpFreqs();

 // use a switch statement here later to switch between the
 // approximation methods:
 //calcButterUnitLpfAna();
 calcPrototypePolesAndZeros();

 frequencyTransform();

 bilinearTransform();

 calcGainFactor();

 calcDirectFormCoeffs();

 // print:
 long i;

 // print direct-form coeffs:
 for(i=0; i<slope; i++)
  Complex(ffCoeffs[i]).print();
 for(i=0; i<slope; i++)
  Complex(fbCoeffs[i]).print();
}



//----------------------------------------------------------------------------
// This function acts merely like a table - it assigns the poles and zeros of
// the analog unit lowpass prototype with hardcoded values:

void IirDesigner::assignPrototypePolesAndZeros()
{
 switch( method )
 {
  case BUTTERWORTH:
  {
   switch( slope )
   {
    case 1: // 1st order Butterworth prototype
    {
     numPrototypePoles = 1;
     numPrototypeZeros = 0;

     prototypePoles[0].re = -1.0; // the single pole
     prototypePoles[0].im =  0.0;
    }
    break;
    case 2: // 2nd order Butterworth prototype
    {
     numPrototypePoles = 2;
     numPrototypeZeros = 0;

     prototypePoles[0].re = -7.071067811865475e-001; // the pole-pair
     prototypePoles[0].im = +7.071067811865476e-001;
     prototypePoles[1]    = prototypePoles[0].getConj();
    }
    break;
    case 3: // 3rd order Butterworth prototype
    {
     numPrototypePoles = 3;
     numPrototypeZeros = 0;

     prototypePoles[0].re = -4.999999999999998e-001; // the pole-pair
     prototypePoles[0].im = +8.660254037844387e-001;
     prototypePoles[1]    = prototypePoles[0].getConj();

     prototypePoles[2].re = -1.0; // the additional single pole
     prototypePoles[2].im =  0.0;
    }
    break;
    case 4: // 4th order Butterworth prototype
    {
     numPrototypePoles = 4;
     numPrototypeZeros = 0;

     prototypePoles[0].re = -3.826834323650897e-001; // 1st pole-pair
     prototypePoles[0].im = +9.238795325112867e-001;
     prototypePoles[1]    = prototypePoles[0].getConj();

     prototypePoles[2].re = -9.238795325112867e-001; // 2nd pole-pair
     prototypePoles[2].im = +3.826834323650899e-001;
     prototypePoles[3]    = prototypePoles[2].getConj();
    }
    break;
    case 5: // 5th order Butterworth prototype
    {
     numPrototypePoles = 5;
     numPrototypeZeros = 0;

     prototypePoles[0].re = -3.090169943749473e-001; // 1st pole-pair
     prototypePoles[0].im = +9.510565162951536e-001;
     prototypePoles[1]    = prototypePoles[0].getConj();

     prototypePoles[2].re = -8.090169943749473e-001; // 2nd pole-pair
     prototypePoles[2].im = +5.877852522924733e-001;
     prototypePoles[3]    = prototypePoles[2].getConj();

     prototypePoles[4].re = -1.0; // additional single pole
     prototypePoles[4].im =  0.0;
    }
    break;
    case 6: // 6th order Butterworth prototype
    {
     numPrototypePoles = 6;
     numPrototypeZeros = 0;

     prototypePoles[0].re = -2.588190451025206e-001; // 1st pole-pair
     prototypePoles[0].im = +9.659258262890683e-001;
     prototypePoles[1]    = prototypePoles[0].getConj();

     prototypePoles[2].re = -7.071067811865475e-001; // 2nd pole-pair
     prototypePoles[2].im = +7.071067811865476e-001;
     prototypePoles[3]    = prototypePoles[2].getConj();

     prototypePoles[4].re = -9.659258262890683e-001; // 3rd pole-pair
     prototypePoles[4].im = +2.588190451025206e-001;
     prototypePoles[5]    = prototypePoles[4].getConj();
    }
    break;

    case 7: // 7th order Butterworth prototype
    {
     numPrototypePoles = 7;
     numPrototypeZeros = 0;

     prototypePoles[0].re = -2.225209339563143e-001; // 1st pole-pair
     prototypePoles[0].im = +9.749279121818236e-001;
     prototypePoles[1]    = prototypePoles[0].getConj();


     prototypePoles[2].re = -6.234898018587335e-001; // 2nd pole-pair
     prototypePoles[2].im = +7.818314824680299e-001;
     prototypePoles[3]    = prototypePoles[2].getConj();

     prototypePoles[4].re = -9.009688679024190e-001; // 3rd pole-pair
     prototypePoles[4].im = +4.338837391175582e-001;
     prototypePoles[5]    = prototypePoles[4].getConj();

     prototypePoles[6].re = -1.0; // additional single pole
     prototypePoles[6].im =  0.0;
    }
    break;
    case 8: // 8th order Butterworth prototype
    {
     numPrototypePoles = 8;
     numPrototypeZeros = 0;

     prototypePoles[0].re = -1.950903220161282e-001; // 1st pole-pair
     prototypePoles[0].im = +9.807852804032304e-001;
     prototypePoles[1]    = prototypePoles[0].getConj();

     prototypePoles[2].re = -5.555702330196020e-001; // 2nd pole-pair
     prototypePoles[2].im = +8.314696123025454e-001;
     prototypePoles[3]    = prototypePoles[2].getConj();


     prototypePoles[4].re = -8.314696123025454e-001; // 3rd pole-pair
     prototypePoles[4].im = +5.555702330196022e-001;
     prototypePoles[5]    = prototypePoles[4].getConj();

     prototypePoles[6].re = -9.807852804032304e-001; // 4th pole-pair
     prototypePoles[6].im = +1.950903220161286e-001;
     prototypePoles[7]    = prototypePoles[6].getConj();
    }
    break;

    //.....

   } // end of "switch( order )"

  } // end of "case BUTTERWORTH"
  break;

 } // end of "switch( method )"
}
