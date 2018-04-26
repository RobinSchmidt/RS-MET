// Construction/Destruction:

template<class TSig, class TPar>
rsLinearPredictor<TSig, TPar>::rsLinearPredictor(int newMaximumOrder)
{
  if( newMaximumOrder >= 1 )
    maxOrder = newMaximumOrder;
  else 
    maxOrder = 128;

  weightVector = new TSig[maxOrder];
  pastInputs   = new TSig[maxOrder];
  updateVector = new TSig[maxOrder];

  order        = maxOrder;
  learnRate    = 0.01;
  forgetRate   = 0.001;
  forgetFactor = 1.0-forgetRate;
  momentum     = 0.5;

  reset();
}

template<class TSig, class TPar>
rsLinearPredictor<TSig, TPar>::~rsLinearPredictor()
{
  if( weightVector != NULL )
    delete[] weightVector;
  if( pastInputs != NULL )
    delete[] pastInputs;
  if( updateVector != NULL )
    delete[] updateVector;
}

// Setup:

template<class TSig, class TPar>
void rsLinearPredictor<TSig, TPar>::setOrder(int newOrder)
{
  if( newOrder >= 1 && newOrder <= maxOrder )
    order = newOrder;
  else
    order = maxOrder;

  reset();
}

template<class TSig, class TPar>
void rsLinearPredictor<TSig, TPar>::setLearnRate(TPar newLearnRate)
{
  if( newLearnRate >= 0.0 && newLearnRate < 1.0 )
    learnRate = newLearnRate;
}

template<class TSig, class TPar>
void rsLinearPredictor<TSig, TPar>::setForgetRate(TPar newForgetRate)
{
  if( newForgetRate >= 0.0 && newForgetRate < 1.0 )
    forgetRate = newForgetRate;
  forgetFactor = 1.0-forgetRate;
}

template<class TSig, class TPar>
void rsLinearPredictor<TSig, TPar>::setMomentum(TPar newMomentum)
{
  if( newMomentum >= 0.0 && newMomentum < 1.0 )
    momentum = newMomentum;
}

// Misc:

template<class TSig, class TPar>
void rsLinearPredictor<TSig, TPar>::reset()
{
  for(int i = 0; i < maxOrder; i++)
  {
    weightVector[i] = 0.0;
    pastInputs[i]   = 0.0;
    updateVector[i] = 0.0;
  }
}
