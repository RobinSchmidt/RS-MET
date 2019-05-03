//#include "rosic_BiquadDesigner.h"
//using namespace rosic;

double BiquadDesigner::getBiquadMagnitudeAt(const double &b0, const double &b1, const double &b2, 
  const double &a1, const double &a2, const double &frequency, const double &sampleRate)
{
  double omega = 2.0 * PI * frequency / sampleRate; // optimize to one mul
  double c1  = cos(omega);
  double c2  = cos(2.0*omega);

  double a1m  = -a1; // we use the other sign-convention in the getSample-function
  double a2m  = -a2;

  double num = b0*b0 + b1*b1   + b2*b2   + 2.0*(b0*b1 + b1*b2)  *c1 + 2.0*b0*b2*c2;
  double den = 1.0   + a1m*a1m + a2m*a2m + 2.0*(  a1m + a1m*a2m)*c1 + 2.0*  a2m*c2;
  double mag = sqrt(num/den);

  return mag;
}

void BiquadDesigner::calculateFirstOrderLowpassCoeffsPrescribedNyquist(double &b0, double &b1, 
                                                                       double &b2, double &a1, 
                                                                       double &a2, 
                                                                       const double &sampleRate, 
                                                                       const double &frequency)
{
  double wc = 2.0*PI*frequency/sampleRate;        // normalized radian cutoff frequency
  double Wp = tan(wc/2.0);                        // pre-warped analog cutoff frequency                     
  double Wn = PI*sampleRate;                      // radian Nyquist frequency
  double Wc = 2.0*PI*frequency;                   // non-pre-warped analog cutoff frequency
  double k2 = 1.0 / (1.0 + ((Wn*Wn)/(Wc*Wc)) );   // gain of prototype at the Nyquist frequency

  // compute analog filter coefficients:
  double A2  = 1.0 / (Wp*Wp*(1.0-2.0*k2));        // A^2
  double B2  = k2*A2;                             // B^2
  double G02 = 1.0;                               // G0^2                          
  double A   = sqrt(A2);
  double B   = sqrt(B2);
  double G0  = sqrt(G02);                         // == 1.0 -> optimize out
  double rD  = 1.0 / (1.0+A);                     // reciprocal of denominator

  // compute digital filter coefficients:
  b0 =  rD * (G0+B);
  b1 =  rD * (G0-B);
  b2 =  0.0;
  a1 = -rD * (1 -A); 
  a2 =  0.0;
}

void BiquadDesigner::calculateFirstOrderHighpassCoeffsPrescribedNyquist(double &b0, double &b1, 
                                                                        double &b2, double &a1, 
                                                                        double &a2, 
                                                                        const double &sampleRate, 
                                                                        const double &frequency)
{
  double wc = 2.0*PI*frequency/sampleRate;
  double Wp = tan(wc/2.0);
  double Wn = PI*sampleRate;
  double Wc = 2.0*PI*frequency;  

  double r   = (Wn*Wn) / (Wc*Wc);
  double k2  = r / (1+r);
  double A2  = -1 / (Wp*Wp*(1-2*k2));
  double B2  = k2*A2;
  double G02 = 0.0;

  double A  = sqrt(A2);
  double B  = sqrt(B2);
  double G0 = sqrt(G02);
  double rD = 1 / (1+A);  

  b0 =  rD * (G0+B);
  b1 =  rD * (G0-B);
  b2 =  0.0;
  a1 = -rD * (1 -A); 
  a2 =  0.0;
}

void BiquadDesigner::calculateFirstOrderHighShelvCoeffsPrescribedNyQuist(double& b0, double& b1, 
                                                                         double& b2, double& a1, 
                                                                         double& a2,       
                                                                         const double& sampleRate, 
                                                                         const double& frequency, 
                                                                         const double& gainFactor)
{
  // todo: catch special case when gainFactor == 1

  double wc = 2.0*PI*frequency/sampleRate;
  double Wp = tan(wc/2.0);
  double Wn = PI*sampleRate;
  double Wc = 2.0*PI*frequency;  

  // protoype coeffs (unmodified):
  double g0 = 1;
  double a  = 1 / (Wc*sqrt(gainFactor));
  double b  = gainFactor*a;

  // desired Nyquist frequency gain (squared):
  double k2 = (g0*g0 + b*b*Wn*Wn) / (1 + a*a*Wn*Wn);

  // modified analog filter coeffs:
  double G0 = g0;
  double A2 = (gainFactor-G0*G0) / (Wp*Wp*(k2-gainFactor));
  double B2 = k2*A2;
  double A  = sqrt(A2);
  double B  = sqrt(B2);

  // reciprocal of denominator:
  double rD =  1.0 / (1.0+A); 

  // digital filter coeffs:
  b0 =  rD * (G0+B);
  b1 =  rD * (G0-B);
  b2 =  0.0;
  a1 = -rD * (1.0 -A);
  a2 =  0.0;
}

void BiquadDesigner::calculateFirstOrderLowShelvCoeffsPrescribedNyQuist(double& b0, double& b1, 
                                                                        double& b2, double& a1, 
                                                                        double& a2,
                                                                        const double& sampleRate, 
                                                                        const double& frequency,                                                                         
                                                                        const double& gainFactor)
{  
  // todo: catch special case when gainFactor == 1

  double wc = 2.0*PI*frequency/sampleRate;
  double Wp = tan(wc/2.0);
  double Wn = PI*sampleRate;
  double Wc = 2.0*PI*frequency;  

  // protoype coeffs (unmodified):
  double g0 = gainFactor;
  double a  = sqrt(gainFactor) / Wc;
  double b  = a;

  // desired Nyquist frequency gain (squared):
  double k2 = (g0*g0 + b*b*Wn*Wn) / (1 + a*a*Wn*Wn);

  // modified analog filter coeffs:
  double G0 = g0;
  double A2 = (gainFactor-G0*G0) / (Wp*Wp*(k2-gainFactor));
  double B2 = k2*A2;
  double A  = sqrt(A2);
  double B  = sqrt(B2);

  // reciprocal of denominator:
  double rD =  1.0 / (1.0+A); 

  // digital filter coeffs:
  b0 =  rD * (G0+B);
  b1 =  rD * (G0-B);
  b2 =  0.0;
  a1 = -rD * (1.0 -A);
  a2 =  0.0;
}
