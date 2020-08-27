//#include "romos_FilterDesign.h"
//#include "../Modules/romos_FilterModules.h"
//using namespace romos;

void biquadBypassCoeffs(double *coeffs)
{
  coeffs[0] = 1.0;    // b0
  coeffs[1] = 0.0;    // b1
  coeffs[2] = 0.0;    // b2
  coeffs[3] = 0.0;    // a1
  coeffs[4] = 0.0;    // a2
}

void biquadLowpassCoeffsBilinear1(double *coeffs, double f)
{
  f = RAPT::rsClip(f, NANO, 0.5*processingStatus.getSystemSampleRate()-NANO);

  // can perhaps be simplified - see DAFX

  double omegaPreWarped = tan(PI * f * processingStatus.getSystemSamplePeriod());
  double pAnalog        = -omegaPreWarped;                                  // pole in the s-plane
  coeffs[3]             = -(1.0+pAnalog)/(1.0-pAnalog);                     // a1
  coeffs[4]             = 0.0;                                              // a2
  double g              = 0.5*sqrt(1.0+coeffs[3]*coeffs[3]+2.0*coeffs[3]);  // gain-factor for normalization at DC 
  coeffs[0]             = g;                                                // b0
  coeffs[1]             = g;                                                // b1
  coeffs[2]             = 0.0;                                              // b2
}

void biquadHighpassCoeffsBilinear1(double *coeffs, double f)
{
  f = RAPT::rsClip(f, NANO, 0.5*processingStatus.getSystemSampleRate()-NANO);

  // can perhaps be simplified - see DAFX

  double omegaPreWarped = tan(PI * f * processingStatus.getSystemSamplePeriod());
  double pAnalog        = -omegaPreWarped;
  coeffs[3]             = -(1.0+pAnalog)/(1.0-pAnalog);  
  coeffs[4]             = 0.0;
  double g              = 0.5*sqrt(1.0+coeffs[3]*coeffs[3]-2.0*coeffs[3]); 
  coeffs[0]             = g;
  coeffs[1]             = -g;  
  coeffs[2]             = 0.0;
}

void biquadLowShelfCoeffsBilinear1(double *coeffs, double f, double g)
{
  f = RAPT::rsClip(f, NANO, 0.5*processingStatus.getSystemSampleRate()-NANO);
  g = RAPT::rsMax(g, NANO);

  double omegaPreWarped = tan(PI * f * processingStatus.getSystemSamplePeriod());
  double A              = sqrt(g) / omegaPreWarped;
  double dr             = 1  / (1+A);  // reciprocal of denominator

  coeffs[0] = dr * (g+A);  // b0
  coeffs[1] = dr * (g-A);  // b1
  coeffs[2] = 0.0;         // b1 
  coeffs[3] = dr * (1-A);  // a1
  coeffs[4] = 0.0;         // a2 
}

void biquadHighShelfCoeffsBilinear1(double *coeffs, double f, double g)
{
  f = RAPT::rsClip(f, NANO, 0.5*processingStatus.getSystemSampleRate()-NANO);
  g = RAPT::rsMax(g, NANO);

  double omegaPreWarped = tan(PI * f * processingStatus.getSystemSamplePeriod());
  double A              = 1.0 / (omegaPreWarped*sqrt(g));
  double B              = g*A;
  double dr             = 1.0 / (1.0+A);  // reciprocal of denominator

  coeffs[0] = dr * (1+B);  // b0
  coeffs[1] = dr * (1-B);  // b1
  coeffs[2] = 0.0;         // b1 
  coeffs[3] = dr * (1-A);  // a1
  coeffs[4] = 0.0;         // a2 
}

void biquadAllpassCoeffsBilinear1(double *coeffs, double f)
{
  double t = tan(PI * f * processingStatus.getSystemSamplePeriod());
  double c = (t-1.0) / (t+1.0);

  coeffs[0] = c;    // b0
  coeffs[1] = 1.0;  // b1
  coeffs[2] = 0.0;  // b2
  coeffs[3] = c;    // a1
  coeffs[4] = 0.0;  // a2
}

void biquadLowpassCoeffsBilinear2(double *coeffs, double f, double q)
{
  f = RAPT::rsClip(f, NANO, 0.5*processingStatus.getSystemSampleRate()-NANO);
  q = RAPT::rsMax(q, NANO);
   
  double omega  = f * processingStatus.getFreqToOmegaFactor();
  double s, c;  
  RAPT::rsSinCos(omega, &s, &c);
  double alpha = s   / (2.0*q);
  double a0r   = 1.0 / (1.0+alpha);

  coeffs[3] = -2.0*c        * a0r;  // a1
  coeffs[4] = -(alpha-1.0)  * a0r;  // a2
  coeffs[1] = (1.0-c)       * a0r;  // b1
  coeffs[0] = 0.5*coeffs[1];        // b0
  coeffs[2] = coeffs[0];            // b2
}

void biquadHighpassCoeffsBilinear2(double *coeffs, double f, double q)
{
  f = RAPT::rsClip(f, NANO, 0.5*processingStatus.getSystemSampleRate()-NANO);
  q = RAPT::rsMax(q, NANO);

  double omega  = f * processingStatus.getFreqToOmegaFactor();
  double s, c;  
  RAPT::rsSinCos(omega, &s, &c);
  double alpha = s   / (2.0*q);
  double a0r   = 1.0 / (1.0+alpha);

  coeffs[3] = -2.0*c         * a0r;  // a1
  coeffs[4] = -(alpha-1.0)   * a0r;  // a2
  coeffs[1] = (1.0+c)        * a0r;  // b1
  coeffs[0] = -0.5*coeffs[1];        // b0
  coeffs[2] = coeffs[0];             // b2
}

void biquadBandpassConstSkirtCoeffs(double *coeffs, double f, double q)
{    
  f = RAPT::rsClip(f, NANO, 0.5*processingStatus.getSystemSampleRate()-NANO);
  q = RAPT::rsMax(q, NANO);

  double omega = f * processingStatus.getFreqToOmegaFactor();
  double s, c;    
  RAPT::rsSinCos(omega, &s, &c);
  double alpha = s   / (2.0*q);
  double a0r   = 1.0 / (1.0+alpha);

  coeffs[3] = - 2.0*c     * a0r;   // a1
  coeffs[4] = (1.0-alpha) * a0r;   // a2 
  coeffs[1] = 0.0;                 // b1
  coeffs[0] = q*alpha     * a0r;   // b0
  coeffs[2] = -coeffs[0];          // b2
}

void biquadBandpassConstPeakCoeffs(double *coeffs, double f, double q)
{    
  f = RAPT::rsClip(f, NANO, 0.5*processingStatus.getSystemSampleRate()-NANO);
  q = RAPT::rsMax(q, NANO);

  double omega = f * processingStatus.getFreqToOmegaFactor();
  double s, c;    
  RAPT::rsSinCos(omega, &s, &c);
  double alpha = s   / (2.0*q);
  double a0r   = 1.0 / (1.0+alpha);

  coeffs[3] = - 2.0*c     * a0r;   // a1
  coeffs[4] = (1.0-alpha) * a0r;   // a2
  coeffs[1] = 0.0;                 // b1
  coeffs[0] = alpha       * a0r;   // b0
  coeffs[2] = -coeffs[0];          // b2
}

void biquadBandrejectCoeffs(double *coeffs, double f, double q)
{
  f = RAPT::rsClip(f, NANO, 0.5*processingStatus.getSystemSampleRate()-NANO);
  q = RAPT::rsMax(q, NANO);

  double omega = f * processingStatus.getFreqToOmegaFactor();
  double s, c;    
  RAPT::rsSinCos(omega, &s, &c);
  double alpha = s   / (2.0*q);
  double a0r   = 1.0 / (1.0+alpha);

  coeffs[3] = -2.0*c       * a0r;  // a1
  coeffs[4] = -(alpha-1.0) * a0r;  // a2
  coeffs[0] = 1.0          * a0r;  // b0
  coeffs[1] = -2.0*c       * a0r;  // b1
  coeffs[2] = 1.0          * a0r;  // b2
}

void biquadPeakCoeffs(double *coeffs, double f, double q, double g)
{
  f = RAPT::rsClip(f, NANO, 0.5*processingStatus.getSystemSampleRate()-NANO);
  q = RAPT::rsMax(q, NANO);
  g = RAPT::rsMax(g, NANO);

  // maybe precompute alpha/A (it's used twice), likewise alpha*g
  double omega = f * processingStatus.getFreqToOmegaFactor();
  double s, c;    
  RAPT::rsSinCos(omega, &s, &c);
  double A     = sqrt(g);
  double alpha = s   / (2.0*q);
  double a0r   = 1.0 / (1.0+alpha/A);

  coeffs[3] = -2.0*c        * a0r;  // a1
  coeffs[4] = (1.0-alpha/A) * a0r;  // a1
  coeffs[0] = (1.0+alpha*A) * a0r;  // b0
  coeffs[1] = -2.0*c        * a0r;  // b1
  coeffs[2] = (1.0-alpha*A) * a0r;  // b2
}

void biquadLowShelfCoeffsBilinear2(double *coeffs, double f, double q, double g)
{
  f = RAPT::rsClip(f, NANO, 0.5*processingStatus.getSystemSampleRate()-NANO);
  q = RAPT::rsMax(q, NANO);
  g = RAPT::rsMax(g, NANO);

  double omega = f * processingStatus.getFreqToOmegaFactor();  
  double s, c;    
  RAPT::rsSinCos(omega, &s, &c);  
  double A     = sqrt(g);
  double beta  = sqrt(A) / q;
  double a0r   = 1.0 / ( (A+1.0) + (A-1.0)*c + beta*s);

  coeffs[3] = - 2.0 *     ( (A-1.0) + (A+1.0)*c          ) * a0r;  // a1
  coeffs[4] =             ( (A+1.0) + (A-1.0)*c - beta*s ) * a0r;  // a2
  coeffs[0] =         A * ( (A+1.0) - (A-1.0)*c + beta*s ) * a0r;  // b0
  coeffs[1] =   2.0 * A * ( (A-1.0) - (A+1.0)*c          ) * a0r;  // b1
  coeffs[2] =         A * ( (A+1.0) - (A-1.0)*c - beta*s ) * a0r;  // b2
}

void biquadHighShelfCoeffsBilinear2(double *coeffs, double f, double q, double g)
{
  f = RAPT::rsClip(f, NANO, 0.5*processingStatus.getSystemSampleRate()-NANO);
  q = RAPT::rsMax(q, NANO);
  g = RAPT::rsMax(g, NANO);

  double omega = f * processingStatus.getFreqToOmegaFactor();  
  double s, c;    
  RAPT::rsSinCos(omega, &s, &c);  
  double A     = sqrt(g);
  double beta  = sqrt(g) / q;
  double a0r   = 1.0 / ( (A+1.0) - (A-1.0)*c + beta*s);

  coeffs[3] =  2.0 *     ( (A-1.0) - (A+1.0)*c          ) * a0r;  // a1
  coeffs[4] =            ( (A+1.0) - (A-1.0)*c - beta*s ) * a0r;  // a2
  coeffs[0] =        A * ( (A+1.0) + (A-1.0)*c + beta*s ) * a0r;  // b0
  coeffs[1] = -2.0 * A * ( (A-1.0) + (A+1.0)*c          ) * a0r;  // b1
  coeffs[2] =        A * ( (A+1.0) + (A-1.0)*c - beta*s ) * a0r;  // b2
}

void biquadAllpassCoeffsBilinear2(double *coeffs, double f, double q)
{
  f = RAPT::rsClip(f, NANO, 0.5*processingStatus.getSystemSampleRate()-NANO);
  q = RAPT::rsMax(q, NANO);

  double omega = f * processingStatus.getFreqToOmegaFactor();  
  double s, c;    
  RAPT::rsSinCos(omega, &s, &c);  
  double alpha = s   / (2.0*q);
  double a0r   = 1.0 / (1.0+alpha);

  coeffs[3] = -2.0*c      * a0r;
  coeffs[4] = (1.0-alpha) * a0r;
  coeffs[0] = (1.0-alpha) * a0r;
  coeffs[1] = (-2.0*c)    * a0r;
  coeffs[2] = (1.0+alpha) * a0r;
}

void ladderCoeffs(double *coeffs, int mode, double f, double r, double ag)
{
  f  = RAPT::rsClip(f,  NANO, 0.5*processingStatus.getSystemSampleRate()-NANO);
  ag = RAPT::rsMax(ag, 0.0);

  // calculate intermediate variables:
  double wc = f * processingStatus.getFreqToOmegaFactor();
  double s, c;
  RAPT::rsSinCos(wc, &s, &c);             // c = cos(wc); s = sin(wc);
  double t = tan(0.25*(wc-PI));

  // map resonance to more intuitive behavior:
  r = (1.0-exp(-3.0*r)) / (1.0-exp(-3.0)); // use rational approximation later - for example: r = (3.0*r) / (1.32*r+0.68*r^2+1.0);

  // maybe get rid of this double-calculation and weighted sum:

  // calculate filter a1-coefficient tuned such the resonance frequency is just right:
  double a1_fullRes = t / (s-c*t);

  // calculate filter a1-coefficient as if there were no resonance:
  double x        = exp(-wc);
  double a1_noRes = -x;

  // use a weighted sum between the resonance-tuned and no-resonance coefficient:
  double a1 = r*a1_fullRes + (1.0-r)*a1_noRes;







  // calculate the b0-coefficient from the condition that each stage should be a leaky
  // integrator:
  double b0 = 1.0+a1;

  // calculate feedback factor by dividing the resonance parameter by the magnitude at the
  // resonant frequency:
  double gsq = b0*b0 / (1.0 + a1*a1 + 2.0*a1*c);
  double k   = r / (gsq*gsq);


  // compute gain compensation:
  double compensationFactor = 1.0;
  if( ag != 0.0 )
  {
    // evaluate the magnitude response at DC:
    double b0_4   = b0*b0*b0*b0;
    double dcGain = b0_4 / ( (((a1+4.0)*a1+6.0)*a1+4.0)*a1 + k*b0_4 + 1.0);

    // derive the required makeup gain from the dc-gain (with optimized calculation for ag == 1.0):
    if( ag == 1.0 )
      compensationFactor = 1.0 / dcGain;
    else
      compensationFactor = pow(dcGain, -ag);
  }

  // compute weighting coefficients for the taps that are required to obtain the different modes:
  double c0, c1, c2, c3, c4;
  switch(mode)
  {
  case romos::LadderFilter::FLAT:      c0 =  1.0; c1 =  0.0; c2 =  0.0; c3 =  0.0; c4 =  0.0;  break;
  case romos::LadderFilter::LP_6:      c0 =  0.0; c1 =  1.0; c2 =  0.0; c3 =  0.0; c4 =  0.0;  break;
  case romos::LadderFilter::LP_12:     c0 =  0.0; c1 =  0.0; c2 =  1.0; c3 =  0.0; c4 =  0.0;  break;
  case romos::LadderFilter::LP_18:     c0 =  0.0; c1 =  0.0; c2 =  0.0; c3 =  1.0; c4 =  0.0;  break;
  case romos::LadderFilter::LP_24:     c0 =  0.0; c1 =  0.0; c2 =  0.0; c3 =  0.0; c4 =  1.0;  break;
  case romos::LadderFilter::HP_6:      c0 =  1.0; c1 = -1.0; c2 =  0.0; c3 =  0.0; c4 =  0.0;  break;
  case romos::LadderFilter::HP_12:     c0 =  1.0; c1 = -2.0; c2 =  1.0; c3 =  0.0; c4 =  0.0;  break;
  case romos::LadderFilter::HP_18:     c0 =  1.0; c1 = -3.0; c2 =  3.0; c3 = -1.0; c4 =  0.0;  break;
  case romos::LadderFilter::HP_24:     c0 =  1.0; c1 = -4.0; c2 =  6.0; c3 = -4.0; c4 =  1.0;  break;  
  case romos::LadderFilter::BP_12_12:  c0 =  0.0; c1 =  0.0; c2 =  1.0; c3 = -2.0; c4 =  1.0;  break;
  case romos::LadderFilter::BP_6_18:   c0 =  0.0; c1 =  0.0; c2 =  0.0; c3 =  1.0; c4 = -1.0;  break;
  case romos::LadderFilter::BP_18_6:   c0 =  0.0; c1 =  1.0; c2 = -3.0; c3 =  3.0; c4 = -1.0;  break;  
  case romos::LadderFilter::BP_6_12:   c0 =  0.0; c1 =  0.0; c2 =  1.0; c3 = -1.0; c4 =  0.0;  break;
  case romos::LadderFilter::BP_12_6:   c0 =  0.0; c1 =  1.0; c2 = -2.0; c3 =  1.0; c4 =  0.0;  break;
  case romos::LadderFilter::BP_6_6:    c0 =  0.0; c1 =  1.0; c2 = -1.0; c3 =  0.0; c4 =  0.0;  break;
  default:                             c0 =  1.0; c1 =  0.0; c2 =  0.0; c3 =  0.0; c4 =  0.0;  // flat
  }

  // assign output variables:
  coeffs[0] = a1;
  coeffs[1] = k;
  coeffs[2] = compensationFactor * c0;
  coeffs[3] = compensationFactor * c1;
  coeffs[4] = compensationFactor * c2;
  coeffs[5] = compensationFactor * c3;
  coeffs[6] = compensationFactor * c4;
}


