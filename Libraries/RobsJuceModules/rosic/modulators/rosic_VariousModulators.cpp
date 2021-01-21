//-------------------------------------------------------------------------------------------------
// construction/destruction:

AttackDecayEnvelope::AttackDecayEnvelope()
{  
  ca = 0.0;
  cd = 0.0;    
  ya = 0.0;
  yd = 0.0;
  n  = 0;  
  np = 0;
  ta = 0.1;
  td = 1.5;     
  tp = 0.5;         
  fs = 44100.0;

  calculateAttackTimeConstant();
  calculateCoeffsAndInitValue();
  reset();
}

AttackDecayEnvelope::~AttackDecayEnvelope()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void AttackDecayEnvelope::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    fs = newSampleRate;
}

void AttackDecayEnvelope::setPeakTime(double newPeakTime)
{
  if( newPeakTime >= 0.0 )
  {
    tp = 0.001*newPeakTime;
    calculateAttackTimeConstant();
    calculateCoeffsAndInitValue();
  }
}

void AttackDecayEnvelope::setDecayTimeConstant(double newTimeConstant)
{
  if( newTimeConstant >= 0.0 )
  {
    td = 0.001*newTimeConstant;
    calculateAttackTimeConstant();
    calculateCoeffsAndInitValue();
  }
}

//-------------------------------------------------------------------------------------------------
// others:

void AttackDecayEnvelope::reset()
{
  ya = 0.0;
  yd = ydInit;
}

void AttackDecayEnvelope::trigger(bool startFromCurrentLevel)
{
  yd = ydInit;
  if( startFromCurrentLevel == false )
    ya = 0.0;
}

bool AttackDecayEnvelope::endIsReached(double threshold)
{
  if( n >= np && ya < threshold )
    return true;
  else
    return false;
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void AttackDecayEnvelope::calculateAttackTimeConstant()
{
  // compute an initial guess for ta (the formulas were found by curve-fitting)
  double rp = tp/td;
  double ra;
  if( rp > 1.0 )
    ra = 1.128*exp(0.9837*rp) - 1.355*(0.42*rp);  
  else
    ra = 0.2073*rp*rp*rp + 0.5512*rp*rp + 0.2398*rp;
  ta = ra*td;

  // refine the initial guess via one or another fixed point iteration:
  double taOld;
  if( tp > 0.0 && tp <= td )
  {
    double k = -log(td)-tp/td; 
    for(int i=0; i<10000; i++)
    {
      taOld = ta;
      ta    = -tp / (k+log(ta));
      if( abs(taOld-ta) < 1.e-12 )
        break;
    }
  }
  else if( tp > 0.0 && tp > td )
  {
    double k  = tp/td + log(td);
    for(int i=0; i<10000; i++)
    {
      taOld = ta;
      ta    = exp(k-tp/ta);
      if( abs(taOld-ta) < 1.e-12 )
        break;
    }
  }
  else
    ta = 0.0;
}

void AttackDecayEnvelope::calculateCoeffsAndInitValue()
{
  // compute the filter coefficients:
  double x  = exp( -1.0 / (fs*td)  );
  double bd = 1-x;
  double ad = -x;
  x         = exp( -1.0 / (fs*ta)  );
  double ba = 1-x;
  double aa = -x;

  // compute the normalizer:
  double normalizer = 1, xp;
  if( td == 0 )
  {
    // zero decay time not allowed
  }
  else if( ta == 0 )
  {
    np         = 0;
    normalizer = 1/bd;
  }
  else if( ta == td )
  {
    np = roundToInt(tp*fs);  // maybe we should use a double for np
    xp = (np+1)*ba*ba*pow(aa, (double)np);
    normalizer = 1/xp; 
  }
  else
  {
    double s   = 1 / (aa-ad);
    double b01 = s * aa*ba*bd;
    double b02 = s * ad*ba*bd;
    double a01 = s * (ad-aa)*aa;
    double a02 = s * (ad-aa)*ad;
    np         = roundToInt(tp*fs);
    xp         = b01*pow(a01, (double)np) - b02*pow(a02, (double)np);
    normalizer = 1/xp;
  }

  cd     = -ad;
  ca     = -aa;
  ydInit = bd*normalizer/cd;
}

//=================================================================================================

void rsTriSawModulator::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  updateOscParameters();
}

void rsTriSawModulator::setAttackTime(double newAttack)
{
  attack = 0.001*newAttack;
  updateOscParameters();
}

void rsTriSawModulator::setDecayTime(double newDecay)
{
  decay = 0.001*newDecay;
  updateOscParameters();
}

void rsTriSawModulator::setTimeScaler(double newScaler)
{
  timeScale = newScaler;
  updateOscParameters();
}

void rsTriSawModulator::updateOscParameters()
{
  double period = timeScale*(attack+decay);
  double freq   = 1/period;
  setPhaseIncrement(freq / sampleRate);
  setAsymmetry(1 - 2*(timeScale*decay)/period);
}