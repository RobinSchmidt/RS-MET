#include "rosic_FormantPreserver.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

FormantPreserver::FormantPreserver(int newMaxOrder) : FormantRemover(newMaxOrder)
{
  pastInputsL  = new double[maxOrder];
  pastInputsR  = new double[maxOrder];
  pastOutputsL = new double[maxOrder];
  pastOutputsR = new double[maxOrder];

  reset();
}

FormantPreserver::~FormantPreserver()
{
  delete[] pastInputsL;
  delete[] pastInputsR;
  delete[] pastOutputsL;
  delete[] pastOutputsR;
}

//-------------------------------------------------------------------------------------------------
// others:

void FormantPreserver::reset()
{
  for(int i=0; i<maxOrder; i++)
  {
    pastInputsL[i]  = 0.0;
    pastInputsR[i]  = 0.0;
    pastOutputsL[i] = 0.0;
    pastOutputsR[i] = 0.0;
  }
  preEmphPastInL = 0.0;
  preEmphPastInR = 0.0;
  deEmphPastOutL = 0.0;
  deEmphPastOutR = 0.0;
  FormantRemover::reset();
}

