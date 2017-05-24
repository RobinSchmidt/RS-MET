//#include <limits>
//#include "rosic_ModalSynthesizer.h"
//using namespace rosic;

//=================================================================================================
// class ModalFilter:

ModalFilter::ModalFilter()
{
  g  = 1.0;
  b1 = a1 = a2 = 0.0;
  reset();
}

void ModalFilter::setModalParameters(double frequency, double amplitude, double decayTime, 
                                     double startPhase, double sampleRate)
{
  // calculate intermediate variables:
  double w      = 2*PI*frequency/sampleRate;
  double alpha  = 1.0 / (decayTime*sampleRate); 
  double A      = amplitude;
  double phi    = wrapToInterval(degreeToRadiant(startPhase), 0, 2*PI);
  double varphi = phi - PI/2;
  double P      = exp(-alpha);
  double r_im   = 0.5 * tan(varphi);
  double R      = sqrt(0.25 + r_im*r_im);
  
  // calculate coefficients:
  a2 = P*P;
  a1 = -2*P*cos(w);
  b1 = -(P/2) * ( 2*(1-cos(2*w))*r_im + sin(2*w) ) / sin(w);
  g  = A/(2*R);
  if( phi > PI )
    g = -g;
}

void ModalFilter::setAmplitude(double newAmplitude)
{
  if( g < 0.0 )
    g = -newAmplitude;
  else
    g = newAmplitude;
}

void ModalFilter::copyCoefficientsFrom(const ModalFilter &other)
{
  g  = other.g;
  b1 = other.b1;
  a1 = other.a1;
  a2 = other.a2;
}

void ModalFilter::reset()
{
  x1 = y1 = y2 = 0.0;
}








//=================================================================================================
// class ModalFilter:

ModalFilter2::ModalFilter2()
{  
  a  = 0.0; 
  cr = 1.0;
  ci = 0.0;

  phaseMod   = 0.0;
  startPhase = 0.0;
  amplitude  = 1.0;
  reset();
}

void ModalFilter2::setModalParameters(double frequency, double amplitude, double decayTime,                                  
                                      double startPhase, double sampleRate)
{
  // compute recursion coeff:
  double w     = 2*PI*frequency/sampleRate;
  double alpha = 1.0 / (decayTime*sampleRate); 
  double r     = exp(-alpha);
  a.setRadiusAndAngle(r, w);


  this->amplitude  = amplitude;
  this->startPhase = startPhase;


  /*
  // compute output weights:
  sinCos((PI/180.0)*startPhase, &cr, &ci);
  //cr *= amplitude;
  //ci *= amplitude;
  */
}

void ModalFilter2::setAmplitude(double newAmplitude)
{
  amplitude = newAmplitude;
}

void ModalFilter2::setPhaseModulation(double newPhaseModulation)
{
  phaseMod = newPhaseModulation;
}

void ModalFilter2::copyCoefficientsFrom(const ModalFilter2 &other)
{
  a  = other.a;
  cr = other.cr;
  ci = other.ci;
}

void ModalFilter2::reset()
{
  z = 0.0;
}






//=================================================================================================
// class ModalFilterWithAttack:

void ModalFilterWithAttack::setModalParameters(double frequency, double amplitude, 
                                               double attackTime, double decayTime, 
                                               double startPhase, double sampleRate, 
                                               double detuneFactor)
{
  rassert(attackTime < decayTime);  // attackTime >= decayTime will not work (because of math)

  double tau1, tau2, scaler;
  tau1 = decayTime;
  getScalerAndSecondTimeConstant(tau1, attackTime, tau2, scaler);
  amplitude *= scaler;

  modalFilter1.setModalParameters(frequency,              amplitude, tau1, startPhase, sampleRate);
  modalFilter2.setModalParameters(frequency*detuneFactor, amplitude, tau2, startPhase, sampleRate);

  // \todo "detuneFactor" does not really work well because when attack and decay are very similar, 
  // the amplitude explodes (due to a high value of "scaler") - either remove this parameter or 
  // find a way to alleviate this (maybe the amplitude access can be computed - math has to be 
  // worked out) the detune is supposed to introduce some roughness into the transient by detuning 
  // the quickly decaying sinusoid
}

void ModalFilterWithAttack::reset()
{
  modalFilter1.reset();
  modalFilter2.reset();
}

double ModalFilterWithAttack::findDecayScaler(double c)
{
  if( c <= 0.0 || c >= 1.0 )
  {
    rassert(false); // assumes 0 < c < 1
    return 1.0;
  }

  // precomputations:
  double kp  = 1/c;                                    // location of the the peak of g(k)
  double k   = 1 + 2*(kp-1);                           // initial guess for the zero of g(k)
  double eps = std::numeric_limits<double>::epsilon(); // relative tolerance
  int    i   = 0;                                      // iteration counter

  // Newton iteration:
  double g, gp;      // g(k), g'(k)
  double kOld = 2*k; // ensure to enter the loop
  while( fabs(k-kOld) > k*eps && i < 1000 )
  {
    kOld = k;
    g    = log(k) + c*(1-k); // g(k)
    gp   = 1/k - c;          // g'(k)
    k    = k - g/gp;         // Newton step
    i++;                     // count iteration
  }
 
  return k;

  // \todo: check this function in the range 0...1, if all works well and the iteration count is 
  // always low, get rid of the iteration counter - it serves a purpose only during development
}

void ModalFilterWithAttack::getScalerAndSecondTimeConstant(double tau1, double peakTime, 
                                                           double &tau2, double &scaler)
{
  if( peakTime >= tau1 )
  {
    rassert(false); // assumes peakTime < tau1
    tau2   = tau1;
    scaler = 1.0;
    return;
  }
  double alpha1 = 1/tau1;
  double c      = alpha1 * peakTime;
  double k      = findDecayScaler(c);
  double alpha2 = k*alpha1;
  double hp     = exp(-alpha1*peakTime) - exp(-alpha2*peakTime); // peak height
  tau2   = 1/alpha2;
  scaler = 1/hp; 
}


//=================================================================================================
// class ModalSynthesizer:

//-------------------------------------------------------------------------------------------------
// construction/destruction:

ModalSynthesizer::ModalSynthesizer()
{
  sampleRate = 44100.0;
  numModes   = maxNumModes;
  modalFilters.reserve(maxNumModes);
  for(int m=0; m<maxNumModes; m++)
    modalFilters.push_back(ModalFilter());
}

ModalSynthesizer::~ModalSynthesizer()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void ModalSynthesizer::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  calculateModalFilterCoefficients();
}

void ModalSynthesizer::setModalParameters(rosic::Vector newFrequencies, 
  rosic::Vector newAmplitudes, rosic::Vector newDecayTimes, rosic::Vector newStartPhases)
{
  frequencies = newFrequencies;
  amplitudes  = newAmplitudes;
  decayTimes  = newDecayTimes;
  startPhases = newStartPhases;
}

//-------------------------------------------------------------------------------------------------
// others:

void ModalSynthesizer::resetModalFilters()
{
  for(unsigned int m = 0; m < modalFilters.size(); m++)
    modalFilters[m].reset();
}

void ModalSynthesizer::calculateModalFilterCoefficients()
{
  int nm = rmin(numModes, frequencies.dim, amplitudes.dim, decayTimes.dim);
  nm = rmin(nm, startPhases.dim);
  for(int m=0; m<nm; m++)
  {
    modalFilters[m].setModalParameters(frequencies.v[m], amplitudes.v[m], decayTimes.v[m], 
      startPhases.v[m], sampleRate); 
  }
}
