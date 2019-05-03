//#include "rosic_FormantRemover.h"
//using namespace rosic;

FormantRemover::FormantRemover(int newMaxOrder) : LinearPredictor(newMaxOrder)
{
  reset();
}

void FormantRemover::reset()
{
  coeff   = 0.0;
  pastIn  = 0.0;
  pastOut = 0.0;
  LinearPredictor::reset();
}
