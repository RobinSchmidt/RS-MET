template<class TSig, class TPar>
rsSmoothingFilter<TSig, TPar>::rsSmoothingFilter()
{
  y1.resize(1);
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
  order = rsMax(1, newOrder);

  y1.resize(order);

  reset();
  // todo: if newOrder > oldOrder, init only the vector values in y1 above oldOrder-1 to 0
  // not all of them

  updateCoeff();
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::reset()
{
  for(int i = 0; i < order; i++)
    y1[i] = 0;
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::updateCoeff()
{
  coeff = exp(-order/decay); // amounts to divide the time-constant by the order
}


/*

ToDo:
-make an envelope follower based on this filter for use in dynamics processors
-the env-follower should be availbale as modulation source in chainer

*/