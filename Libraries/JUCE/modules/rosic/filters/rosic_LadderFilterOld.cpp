namespace rosic // temporary - for as long as there's a RAPT::LadderFilter class, too
{               // maybe at some point, rename this class into rsLadderFilterDD as template
                // instantiation of RAPT::rsLadderFilter<double, double>

//-------------------------------------------------------------------------------------------------
// construction/destruction:

LadderFilterOld::LadderFilterOld()
{
  parameters = new LadderFilterParameters;
  isMaster   = true;

  cutoff = 1000.0;
  b0 = 1.0;
  a1 = 0.0;
  k  = 0.0;
  c0 = 0.0;
  c1 = 0.0;
  c2 = 0.0;
  c3 = 0.0;
  c4 = 1.0;

  makeupGain    = 1.0;
  //makeupGainRec = 1.0;
  //makeupGainSq  = 1.0;

  allpass.setSampleRate(parameters->sampleRate);
  allpass.setMode(OnePoleFilterStereo::ALLPASS);
  //allpass.setCutoff(0.25*cutoff);
  allpass.setCutoff(20000.0);
  allpass.setShelvingGain(1.0);

  calculateCoefficients();
  reset();
}

LadderFilterOld::~LadderFilterOld()
{
  if(isMaster && parameters != NULL)
    delete parameters;
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void LadderFilterOld::setSampleRate(double newSampleRate)
{
  if(newSampleRate > 0.0)
    parameters->sampleRate = newSampleRate;
  allpass.setSampleRate(parameters->sampleRate);

  calculateCoefficients();
  for(unsigned int s=0; s<slaves.size(); s++)
    slaves[s]->calculateCoefficients();
}

void LadderFilterOld::setAllpassFreq(double newAllpassFreq)
{
  parameters->allpassFreq = newAllpassFreq;
  allpass.setCutoff(newAllpassFreq);
  for(unsigned int s=0; s<slaves.size(); s++)
    slaves[s]->setAllpassFreq(newAllpassFreq);
}

void LadderFilterOld::setMakeUp(double newMakeUp, bool updateCoefficients)
{
  parameters->makeUp = newMakeUp;

  if(updateCoefficients == true)
  {
    calculateCoefficients();
    for(unsigned int s=0; s<slaves.size(); s++)
      slaves[s]->calculateCoefficients();
  }
}

void LadderFilterOld::setDrive(double newDrive)
{
  parameters->driveFactor = RAPT::rsDbToAmp(newDrive);
}

void LadderFilterOld::setDcOffset(double newDcOffset)
{
  parameters->dcOffset = newDcOffset;
}

void LadderFilterOld::setOutputStage(int newOutputStage)
{
  if(newOutputStage >= 0 && newOutputStage <= 4)
    parameters->outputStage = newOutputStage;
  else
    DEBUG_BREAK; // ouput stage must be in 0...4
}

void LadderFilterOld::setMode(int newMode)
{
  if(newMode >= 0 && newMode < LadderFilterParameters::NUM_MODES)
    parameters->mode = newMode;
}

void LadderFilterOld::setMorph(double newMorph)
{
  parameters->morph = newMorph;

/*
c0 = [ 1  0  0  0  0];
c1 = [-4  1  0  0  0];
c2 = [ 6 -3  1  0  0];
c3 = [-4  3 -2  1  0];
c4 = [ 1 -1  1 -1  1];
*/

  // may be streamlined
  double m = newMorph;
  if(m < 0.25)
  {
    c0 =  1.0 + 4.0*m* (0.0 -  1.0);
    c1 = -4.0 + 4.0*m* (1.0 - -4.0);
    c2 =  6.0 + 4.0*m* (-3.0 -  6.0);
    c3 = -4.0 + 4.0*m* (3.0 - -4.0);
    c4 =  1.0 + 4.0*m* (-1.0 -  1.0);
  }
  else if(m < 0.5)
  {
    m  -= 0.25;
    c0  =  0.0 + 4.0*m* (0.0 -  0.0);
    c1  =  1.0 + 4.0*m* (0.0 -  1.0);
    c2  = -3.0 + 4.0*m* (1.0 - -3.0);
    c3  =  3.0 + 4.0*m* (-2.0 -  3.0);
    c4  = -1.0 + 4.0*m* (1.0 - -1.0);
  }
  else if(m < 0.75)
  {
    m  -= 0.5;
    c0  =  0.0 + 4.0*m* (0.0 -  0.0);
    c1  =  0.0 + 4.0*m* (0.0 -  0.0);
    c2  =  1.0 + 4.0*m* (0.0 -  1.0);
    c3  = -2.0 + 4.0*m* (1.0 - -2.0);
    c4  =  1.0 + 4.0*m* (-1.0 -  1.0);
  }
  else if(m < 1.0)
  {
    m  -= 0.75;
    c0  =  0.0 + 4.0*m* (0.0 -  0.0);
    c1  =  0.0 + 4.0*m* (0.0 -  0.0);
    c2  =  0.0 + 4.0*m* (0.0 -  0.0);
    c3  =  1.0 + 4.0*m* (0.0 -  1.0);
    c4  = -1.0 + 4.0*m* (1.0 - -1.0);
  }
  else
  {
    c0 = c1 = c2 = c3 = 0.0;
    c4 = 1.0;
  }

  // todo: try a polynomial between the breakpoints for more natural transition
}

//-------------------------------------------------------------------------------------------------
// inquiry:

Complex LadderFilterOld::getTransferFunctionAt(Complex z, bool withFeedback, bool /*withMakeUpBoost*/,
  bool withMakeUpGain, int stage)
{
  Complex G1, G4, G, H;

  G1 = b0 / (1 + a1*(1/z));             // G1: response of one stage
  G4 = G1*G1*G1*G1;                     // G4: response of four stages (without feedback)
  switch(stage)
  {
  case 0: G = 1.0;       break;
  case 1: G = G1;        break;
  case 2: G = G1*G1;     break;
  case 3: G = G1*G1*G1;  break;
  case 4: G = G4;        break;
  }                                     // G: response of the requested stage (without feedback)

  if(withFeedback == true)
    H  = G / (1 + k*(1/z)*G4);         // response of four stages with feedback
  else
    H = G;

  //if( withMakeUpBoost == true )
  //  H *= lowBooster.getTransferFunctionAt(z);

  if(withMakeUpGain == true)
    H *= makeupGain;

  return H;
}

double LadderFilterOld::getMagnitudeAt(double frequency, bool withFeedback, bool withMakeUpBoost,
  bool withMakeUpGain, int stage)
{
  double omega = 2*PI*frequency / parameters->sampleRate;

  Complex z = expC(Complex(0.0, omega));
  Complex H = getTransferFunctionAt(z, withFeedback, withMakeUpBoost, withMakeUpGain, stage);

  return H.getRadius();
}

void LadderFilterOld::getMagnitudeResponse(double *frequencies, double *magnitudes, int numBins,
  bool inDecibels, bool accumulate)
{
  int k;
  int s = parameters->outputStage;
  if(inDecibels == false)
  {
    if(accumulate == false)
    {
      for(k=0; k<numBins; k++)
        magnitudes[k] = getMagnitudeAt(frequencies[k], true, true, true, s);
    }
    else
    {
      for(k=0; k<numBins; k++)
        magnitudes[k] *= getMagnitudeAt(frequencies[k], true, true, true, s);
    }
  }
  else
  {
    if(accumulate == false)
    {
      for(k=0; k<numBins; k++)
        magnitudes[k] = RAPT::rsAmpToDbWithCheck(getMagnitudeAt(frequencies[k], true, true, true, s), 0.000001);
    }
    else
    {
      for(k=0; k<numBins; k++)
        magnitudes[k] += RAPT::rsAmpToDbWithCheck(getMagnitudeAt(frequencies[k], true, true, true, s), 0.000001);
    }
  }
}

double LadderFilterOld::getCutoff()
{
  return cutoff;
}

double LadderFilterOld::getResonance()
{
  return parameters->resonanceRaw;
}

double LadderFilterOld::getDrive()
{
  return RAPT::rsAmpToDb(parameters->driveFactor);
}

int LadderFilterOld::getOutputStage()
{
  return parameters->outputStage;
}

double LadderFilterOld::getAllpassFreq()
{
  return parameters->allpassFreq;
}

double LadderFilterOld::getMakeUp()
{
  return parameters->makeUp;
}

//-------------------------------------------------------------------------------------------------
// master/slave config:

void LadderFilterOld::addSlave(LadderFilterOld* newSlave)
{
  // add the new slave to the vector of slaves:
  slaves.push_back(newSlave);

  // delete the original parameter-set of the new slave and redirect it to ours (with some safety
  // checks):
  if(newSlave->parameters != NULL && newSlave->parameters != this->parameters)
  {
    delete newSlave->parameters;
    newSlave->parameters = this->parameters;
  }
  else
  {
    DEBUG_BREAK;
    // the object to be added as slave did not contain a valid parameter-pointer - maybe it has
    // been already added as slave to another master?
  }

  // set the isMaster-flag of the new slave to false:
  newSlave->isMaster = false;

  // this flag will prevent the destructor of the slave from trying to delete the parameter-set
  // which is now shared - only masters delete their parameter-set on destruction
}

//-------------------------------------------------------------------------------------------------
// others:

void LadderFilterOld::reset()
{
  allpass.resetBuffers();
  y1L   = 0.0;
  y2L   = 0.0;
  y3L   = 0.0;
  y4L   = 0.0;
  yOldL = 0.0;
  y1R   = 0.0;
  y2R   = 0.0;
  y3R   = 0.0;
  y4R   = 0.0;
  yOldR = 0.0;

  kOld  = k;
  a1Old = a1;
}


double LadderFilterOld::getSampleTest(double in)
{
  // hmm...maybe we must scale the internal states of the integrators also
  // maybe we need to take b into account

  // apply drive, feedback and DC-offset:
  double y0L = parameters->driveFactor*in - k*yOldL + parameters->dcOffset;
  //double y0L = in;

  if(a1 != a1Old)
  {
    //double b0Old = 1.0+a1Old;

    // cutoff was changed - update internal states to preserve energy
    double scaler = 1.0;
    double offset = 0.0;

    /*
    double eOld, eNew, eG; // energy represented by integrator states
    eOld = y1L*y1L * (1/(1-a1Old*a1Old) - 1);
    eNew = y1L*y1L * (1/(1-a1   *a1   ) - 1);
    */
    //scaler = sqrt( (1/(1-a1Old*a1Old)-1) /  (1/(1-a1*a1)-1+TINY) );
    //scaler = sqrt( (1/(1-a1*a1)-1) /(1/(1-a1Old*a1Old)-1+TINY) );

    //scaler = a1Old/a1;

    //scaler = a1/a1Old;
    //offset = in*(b0-b0Old)/a1;

    //scaler *= scaler;

    y1L = scaler*y1L + offset;
    y2L = scaler*y2L + offset;
    y3L = scaler*y3L + offset;
    y4L = scaler*y4L + offset;
  }

  //double b0 = 1.0;  // only during development
  //y1L = y0L - a1*y1L;

  //int dummy = 0;

  // cascade of 4 1st order sections:
  y1L = b0*y0L - a1*y1L;
  y2L = b0*y1L - a1*y2L;
  y3L = b0*y2L - a1*y3L;
  y4L = b0*y3L - a1*y4L;

  // test some scalings of the feedback value in order to counteract energy changes under
  // time-varying conditions:
  double s = 1.0;
  /*
  if( k != kOld )
  {
    //s  = ((kOld)/(k+TINY));
    //s *= s*s*s;
    s = 0.0; // test
  }
  */
  yOldL = s * y4L;

  kOld  = k;
  a1Old = a1;

  //return y1L;  // preliminary
  return y4L;
}


} // end namespace rosic
