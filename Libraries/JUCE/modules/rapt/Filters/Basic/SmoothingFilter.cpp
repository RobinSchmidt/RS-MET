template<class TSig, class TPar>
rsSmoothingFilter<TSig, TPar>::rsSmoothingFilter()
{
  reset();
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::setTimeConstantAndSampleRate(TPar timeConstant, TPar sampleRate)
{
  decay = sampleRate * timeConstant;
  updateCoeff();
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::setOrder(int newOrder)
{
  order = rsMin(newOrder, maxOrder);
  updateCoeff();
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::reset()
{
  for(int i = 0; i < maxOrder; i++)
    y1[i] = 0;
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::updateCoeff()
{
  coeff = exp(-order/decay); // amounts to divide the time-constant by the order
}