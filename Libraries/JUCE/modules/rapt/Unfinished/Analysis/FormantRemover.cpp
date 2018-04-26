using namespace RSLib;

rsFormantRemover::rsFormantRemover(int newMaxOrder) 
  : rsLinearPredictor(newMaxOrder)
{
  reset();
}

void rsFormantRemover::reset()
{
  coeff   = 0.0;
  pastIn  = 0.0;
  pastOut = 0.0;
  rsLinearPredictor::reset();
}
