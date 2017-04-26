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

void IirDesigner::setSampleRate(flt64  newSampleRate)
{
 if( newSampleRate>0 )
  sampleRate = newSampleRate;
}

void IirDesigner::setMode(int64 newMode)
{
 if( (newMode>=1) && (newMode<=4) )
  mode = newMode;
}

void IirDesigner::setSlope(int64 newSlope)
{
 //if( (newSlope>=1) && (newSlope<=(maxOrder/2)) )
  slope = newSlope;
}

void IirDesigner::setFreq1(flt64 newFreq1)
{
 if( (newFreq1>=0) && (newFreq1<=(0.5*sampleRate)) )
  freq1 = newFreq1;
}

void IirDesigner::setFreq2(flt64 newFreq2)
{
 if( (newFreq2>=0) && (newFreq2<=(0.5*sampleRate)) )
  freq2 = newFreq2;
}

//----------------------------------------------------------------------------
// coefficient calculation:

void IirDesigner::getDirectFormCoeffs(flt64 *FeedforwardCoeffs, 
                                      flt64 *FeedbackCoeffs)
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

void IirDesigner::getBiquadCascadeCoeffs(flt64 *b0, flt64 *b1, flt64 *b2, 
                                                    flt64 *a1, flt64 *a2)
{
 // do nothing in the case, that the filters slope is zero:
 if( slope <= 0)
  return;

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
 int64A k;
 for(k=0; k<numPolesDig; k+=2)
 {
  b0[k/2] = 1.0;
  b1[k/2] = -2*zerosDig[k].re;
  b2[k/2] = zerosDig[k].getMagnitude() * zerosDig[k].getMagnitude();

  a1[k/2] = -2*polesDig[k].re;
  a2[k/2] = polesDig[k].getMagnitude() * polesDig[k].getMagnitude();
 }

 // if there is a singles left over pole (and zero), it will be realized by 
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
 //....
 static flt64A gain, normalizer, omegaDigC;
 static int64A numStages;
 if( isEven((int)numPolesDig) )
  numStages = numPolesDig/2;
 else
  numStages = (numPolesDig+1)/2;
 switch( mode )
 {
 case LOWPASS:
  {
   for(k=0; k<numStages; k++)
   {
    // get gain at DC:
    gain = getBiquadMagnitudeAt(b0[k], b1[k], b2[k], a1[k], a2[k], 0.0);

    // divide the b-coefficients by that factor:
    normalizer = 1.0/gain;
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
    gain = getBiquadMagnitudeAt(b0[k], b1[k], b2[k], a1[k], a2[k], PI);

    // divide the b-coefficients by that factor:
    normalizer = 1.0/gain;
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
    gain = getBiquadMagnitudeAt(b0[k], b1[k], b2[k], a1[k], a2[k], omegaDigC);

    // divide the b-coefficients by that factor:
    normalizer = 1.0/gain;
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
    gain = getBiquadMagnitudeAt(b0[k], b1[k], b2[k], a1[k], a2[k], 0.0);

    // divide the b-coefficients by that factor:
    normalizer = 1.0/gain;
    b0[k]      = normalizer * b0[k];
    b1[k]      = normalizer * b1[k];
    b2[k]      = normalizer * b2[k];
   }
  }
 }  // end of "switch( mode )"


}


/*
void IirDesigner::getBiquadCascadeCoeffs(flt64 *b0, flt64 *b1, flt64 *b2, 
                                         flt64 *a1, flt64 *a2)
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
 numPolesLpfAna = slope;
 numZerosLpfAna = 0;   // no zeros in the numerator of H(s) -> all zeros
                       // of H(s) at infinity

 sample real = 0;  //for the real part of the pole/zero
 sample imag = 0;  //for the imaginary part of the pole/zero

 for(long k=0; k<numPolesLpfAna; k++)
 {
  real = -sin( ((2*(k+1)-1)*PI)/(2*numPolesLpfAna) ); //can be optimized: preCalc 1/(2*numPolesAna) and 2*PI
  imag =  cos( ((2*(k+1)-1)*PI)/(2*numPolesLpfAna) );

  unitLpfPolesAna[k].re = real;
  unitLpfPolesAna[k].im = imag;
 }
}


void IirDesigner::frequencyTransform()
{
 long k;               // for indexing the poles and zeros
 Complex temp1, temp2; // for some intermediate in the BPF and BRF case

 switch(mode)
 {

  case 1: //lowpass->lowpass transform
  {
   numPolesAna = numPolesLpfAna;
   numZerosAna = numZerosLpfAna;
   for(k=0; k<numPolesAna; k++)
    polesAna[k] = unitLpfPolesAna[k] * omegaAna1;
   for(k=0; k<numZerosAna; k++)
    zerosAna[k] = unitLpfZerosAna[k] * omegaAna1;
  }
  break;

  case 2: //lowpass->highpass transform
  {
   numPolesAna = numPolesLpfAna;
   numZerosAna = numZerosLpfAna;
   for(k=0; k<numPolesAna; k++)
    polesAna[k] = Complex(omegaAna1)/unitLpfPolesAna[k];
   for(k=0; k<numZerosAna; k++)
    zerosAna[k] = Complex(omegaAna1)/unitLpfZerosAna[k];
   //if there are no zeros in the numerator (indicated by numZerosAna==0),
   //then all zeros of H(s) are at s=inf - they map to s=0:
   if(numZerosAna==0)
   {
    numZerosAna = numPolesAna;
    for(k=0; k<numZerosAna; k++)
     zerosAna[k] = 0;
   }
  }
  break;

  case 3: //lowpass->bandpass transform
  {
   numPolesAna = 2*numPolesLpfAna; //twice as much poles as in the LPF
   numZerosAna = 2*numZerosLpfAna; 
   for(k=0; k<numPolesLpfAna; k++)
   {
    temp1 = unitLpfPolesAna[k]*bwAna*0.5;
    temp2 = -unitLpfPolesAna[k]*bwAna;
    temp2 = temp2*temp2*0.25;
    temp2 = Complex::sqrt(temp2 - Complex(omegaAnaC*omegaAnaC));
    polesAna[k]               = temp1 + temp2;
    polesAna[numPolesAna-1-k] = temp1 - temp2;
    //k2 = 2*k;
    //polesAna[k2]   = temp1 + temp2;
    //polesAna[k2+1] = temp1 - temp2;
   }
   for(k=0; k<numZerosLpfAna; k++)
   {
    temp1 = unitLpfZerosAna[k]*bwAna*0.5;
    temp2 = -unitLpfZerosAna[k]*bwAna;
    temp2 = temp2*temp2*0.25;
    temp2 = Complex::sqrt(temp2 - Complex(omegaAnaC*omegaAnaC));
    zerosAna[k]               = temp1 + temp2;
    zerosAna[numZerosAna-1-k] = temp1 - temp2;
   }
   //if there are no zeros in the numerator (indicated by numZerosAna==0),
   //then we will get numPolesAna/2 zeros at s=inf and numPolesAna/2 zeros
   //at s=0:
   if(numZerosLpfAna==0)
   {
    numZerosAna = numPolesAna/2;
    for(k=0; k<numZerosAna; k++)
     zerosAna[k] = 0;
   }
  }
  break;

  case 4: //lowpass->bandreject transform
  {
   numPolesAna = 2*numPolesLpfAna; //twice as much poles as in the LPF
   numZerosAna = 2*numZerosLpfAna; 
   for(k=0; k<numPolesLpfAna; k++)
   {
    temp1 = Complex(bwAna*0.5)/unitLpfPolesAna[k];
    temp2 = Complex(-bwAna)/unitLpfPolesAna[k];
    temp2 = temp2*temp2*0.25;
    temp2 = Complex::sqrt(temp2 - Complex(omegaAnaC*omegaAnaC));
    polesAna[k]               = temp1 + temp2;
    polesAna[numPolesAna-1-k] = temp1 - temp2;
    //k2 = 2*k;
    //polesAna[k2]   = temp1 + temp2;
    //polesAna[k2+1] = temp1 - temp2;
   }
   for(k=0; k<numZerosLpfAna; k++)
   {
    temp1 = Complex(bwAna*0.5)/unitLpfZerosAna[k];
    temp2 = Complex(-bwAna)/unitLpfZerosAna[k];
    temp2 = temp2*temp2*0.25;
    temp2 = Complex::sqrt(temp2 - Complex(omegaAnaC*omegaAnaC));
    zerosAna[k]               = temp1 + temp2;
    zerosAna[numZerosAna-1-k] = temp1 - temp2;
   }
   //if there are no zeros in the numerator (indicated by numZerosAna==0),
   //then we will get numPolesAna/2 zeros at s=j*omegaAnaC and numPolesAna/2 zeros
   //at s=-j*omegaAnaC:
   if(numZerosLpfAna==0)
   {
    numZerosAna = numPolesAna;
    for(k=0; k<(numZerosAna/2); k++)
    {
    temp1    = Complex(0,0);
    temp1.im = omegaAnaC;
    zerosAna[k]               =  temp1;
    zerosAna[numZerosAna-1-k] = -temp1;
    }
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
 numZerosDig = numPolesAna; //numZerosAna may contain only half of the zeros in the digital filter
                            //when there are some at s=inf

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
 Complex num(1,0); //accumulates the numerator of H(z)
 Complex den(1,0); //accumulates the denominator of H(z)
 Complex z;        //value of z at which H(z) is evaluated
 sample  temp;

 //decide at which value of z the transfer function H(z) should
 //be evaluated:
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
 }
 //z.print();

 //evaluate H(z):
 static int64A k;
 for(k=0; k<numPolesDig; k++)
  den = den * ( z - polesDig[k] );
 for(k=0; k<numZerosDig; k++)
  num = num * ( z - zerosDig[k] );

 Complex H    = num/den;
 sample  gain = H.getMagnitude();

 gainFactor = 1/gain;
}

flt64 IirDesigner::getBiquadMagnitudeAt(flt64 b0, flt64 b1, flt64 b2, 
                                        flt64 a1, flt64 a2, flt64 omega)
{
 static flt64A num, den, mag, cosOmega, cos2Omega;

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

void IirDesigner::convolve(Complex *Seq1, int64 Length1,
                           Complex *Seq2, int64 Length2,
                           Complex *Result)
{
 int64 lengthR = Length1+Length2-1; //length of the resulting sequence

 //init the result-buffer with zeros:
 int64 n, k;
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
     numPolesLpfAna = 1;
     numZerosLpfAna = 0;

     unitLpfPolesAna[0].re = -1.0; // the single pole
     unitLpfPolesAna[0].im =  0.0;
    }
    break;
    case 2: // 2nd order Butterworth prototype
    {
     numPolesLpfAna = 2;
     numZerosLpfAna = 0;

     unitLpfPolesAna[0].re = -7.071067811865475e-001; // the pole-pair
     unitLpfPolesAna[0].im = +7.071067811865476e-001;
     unitLpfPolesAna[1]    = unitLpfPolesAna[0].getConj();
    }
    break;
    case 3: // 3rd order Butterworth prototype
    {
     numPolesLpfAna = 3;
     numZerosLpfAna = 0;

     unitLpfPolesAna[0].re = -4.999999999999998e-001; // the pole-pair
     unitLpfPolesAna[0].im = +8.660254037844387e-001;
     unitLpfPolesAna[1]    = unitLpfPolesAna[0].getConj();

     unitLpfPolesAna[2].re = -1.0; // the additional single pole
     unitLpfPolesAna[2].im =  0.0;
    }
    break;
    case 4: // 4th order Butterworth prototype
    {
     numPolesLpfAna = 4;
     numZerosLpfAna = 0;

     unitLpfPolesAna[0].re = -3.826834323650897e-001; // 1st pole-pair
     unitLpfPolesAna[0].im = +9.238795325112867e-001;
     unitLpfPolesAna[1]    = unitLpfPolesAna[0].getConj();

     unitLpfPolesAna[2].re = -9.238795325112867e-001; // 2nd pole-pair
     unitLpfPolesAna[2].im = +3.826834323650899e-001;
     unitLpfPolesAna[3]    = unitLpfPolesAna[2].getConj();
    }
    break;
    case 5: // 5th order Butterworth prototype
    {
     numPolesLpfAna = 5;
     numZerosLpfAna = 0;

     unitLpfPolesAna[0].re = -3.090169943749473e-001; // 1st pole-pair
     unitLpfPolesAna[0].im = +9.510565162951536e-001;
     unitLpfPolesAna[1]    = unitLpfPolesAna[0].getConj();

     unitLpfPolesAna[2].re = -8.090169943749473e-001; // 2nd pole-pair
     unitLpfPolesAna[2].im = +5.877852522924733e-001;
     unitLpfPolesAna[3]    = unitLpfPolesAna[2].getConj();

     unitLpfPolesAna[4].re = -1.0; // additional single pole
     unitLpfPolesAna[4].im =  0.0;
    }
    break;
    case 6: // 6th order Butterworth prototype
    {
     numPolesLpfAna = 6;
     numZerosLpfAna = 0;

     unitLpfPolesAna[0].re = -2.588190451025206e-001; // 1st pole-pair
     unitLpfPolesAna[0].im = +9.659258262890683e-001;
     unitLpfPolesAna[1]    = unitLpfPolesAna[0].getConj();

     unitLpfPolesAna[2].re = -7.071067811865475e-001; // 2nd pole-pair
     unitLpfPolesAna[2].im = +7.071067811865476e-001;
     unitLpfPolesAna[3]    = unitLpfPolesAna[2].getConj();

     unitLpfPolesAna[4].re = -9.659258262890683e-001; // 3rd pole-pair
     unitLpfPolesAna[4].im = +2.588190451025206e-001;
     unitLpfPolesAna[5]    = unitLpfPolesAna[4].getConj();
    }
    break;

    case 7: // 7th order Butterworth prototype
    {
     numPolesLpfAna = 7;
     numZerosLpfAna = 0;

     unitLpfPolesAna[0].re = -2.225209339563143e-001; // 1st pole-pair
     unitLpfPolesAna[0].im = +9.749279121818236e-001;
     unitLpfPolesAna[1]    = unitLpfPolesAna[0].getConj();


     unitLpfPolesAna[2].re = -6.234898018587335e-001; // 2nd pole-pair
     unitLpfPolesAna[2].im = +7.818314824680299e-001;
     unitLpfPolesAna[3]    = unitLpfPolesAna[2].getConj();

     unitLpfPolesAna[4].re = -9.009688679024190e-001; // 3rd pole-pair
     unitLpfPolesAna[4].im = +4.338837391175582e-001;
     unitLpfPolesAna[5]    = unitLpfPolesAna[4].getConj();

     unitLpfPolesAna[6].re = -1.0; // additional single pole
     unitLpfPolesAna[6].im =  0.0;
    }
    break;
    case 8: // 8th order Butterworth prototype
    {
     numPolesLpfAna = 8;
     numZerosLpfAna = 0;

     unitLpfPolesAna[0].re = -1.950903220161282e-001; // 1st pole-pair
     unitLpfPolesAna[0].im = +9.807852804032304e-001;
     unitLpfPolesAna[1]    = unitLpfPolesAna[0].getConj();

     unitLpfPolesAna[2].re = -5.555702330196020e-001; // 2nd pole-pair
     unitLpfPolesAna[2].im = +8.314696123025454e-001;
     unitLpfPolesAna[3]    = unitLpfPolesAna[2].getConj();


     unitLpfPolesAna[4].re = -8.314696123025454e-001; // 3rd pole-pair
     unitLpfPolesAna[4].im = +5.555702330196022e-001;
     unitLpfPolesAna[5]    = unitLpfPolesAna[4].getConj();

     unitLpfPolesAna[6].re = -9.807852804032304e-001; // 4th pole-pair
     unitLpfPolesAna[6].im = +1.950903220161286e-001;
     unitLpfPolesAna[7]    = unitLpfPolesAna[6].getConj();
    }
    break;

    //.....

   } // end of "switch( order )"

  } // end of "case BUTTERWORTH"
  break;

 } // end of "switch( method )"
}
