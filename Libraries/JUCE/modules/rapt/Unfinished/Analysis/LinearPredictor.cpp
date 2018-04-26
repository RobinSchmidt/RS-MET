using namespace RSLib;

// Construction/Destruction:

rsLinearPredictor::rsLinearPredictor(int newMaximumOrder)
{
  if( newMaximumOrder >= 1 )
    maxOrder = newMaximumOrder;
  else 
    maxOrder = 128;

  weightVector = new double[maxOrder];
  pastInputs   = new double[maxOrder];
  updateVector = new double[maxOrder];

  order        = maxOrder;
  learnRate    = 0.01;
  forgetRate   = 0.001;
  forgetFactor = 1.0-forgetRate;
  momentum     = 0.5;

  reset();
}

rsLinearPredictor::~rsLinearPredictor()
{
  if( weightVector != NULL )
    delete[] weightVector;
  if( pastInputs != NULL )
    delete[] pastInputs;
  if( updateVector != NULL )
    delete[] updateVector;
}

// Setup:

void rsLinearPredictor::setOrder(int newOrder)
{
  if( newOrder >= 1 && newOrder <= maxOrder )
    order = newOrder;
  else
    order = maxOrder;

  reset();
}

void rsLinearPredictor::setLearnRate(double newLearnRate)
{
  if( newLearnRate >= 0.0 && newLearnRate < 1.0 )
    learnRate = newLearnRate;
}

void rsLinearPredictor::setForgetRate(double newForgetRate)
{
  if( newForgetRate >= 0.0 && newForgetRate < 1.0 )
    forgetRate = newForgetRate;
  forgetFactor = 1.0-forgetRate;
}

void rsLinearPredictor::setMomentum(double newMomentum)
{
  if( newMomentum >= 0.0 && newMomentum < 1.0 )
    momentum = newMomentum;
}

// Misc:

void rsLinearPredictor::reset()
{
  for(int i = 0; i < maxOrder; i++)
  {
    weightVector[i] = 0.0;
    pastInputs[i]   = 0.0;
    updateVector[i] = 0.0;
  }
}
