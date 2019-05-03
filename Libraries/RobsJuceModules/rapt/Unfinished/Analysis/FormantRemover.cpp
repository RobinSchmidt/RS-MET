template<class TSig, class TPar>
rsFormantRemover<TSig, TPar>::rsFormantRemover(int newMaxOrder)
  : rsLinearPredictor<TSig, TPar>(newMaxOrder)
{
  reset();
}

template<class TSig, class TPar>
void rsFormantRemover<TSig, TPar>::reset()
{
  coeff   = 0.0;
  pastIn  = 0.0;
  pastOut = 0.0;
  rsLinearPredictor<TSig, TPar>::reset();
}
