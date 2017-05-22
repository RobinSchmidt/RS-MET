//#include "rosic_EllipticSubBandFilterDirectForm.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

EllipticSubBandFilterDirectForm::EllipticSubBandFilterDirectForm()
{
  reset();  
  setSubDivision(2.0); 
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void EllipticSubBandFilterDirectForm::setSubDivision(double newSubDivision)
{
  EllipticSubBandFilter tmpBiquad;
  tmpBiquad.setSubDivision(newSubDivision);
  FilterCoefficientConverter::biquadCascadeToDirectForm(6, tmpBiquad.getAddressB0(), tmpBiquad.getAddressB1(), 
    tmpBiquad.getAddressB2(), tmpBiquad.getAddressA1(), tmpBiquad.getAddressA2(), b, a);
}

void EllipticSubBandFilterDirectForm::reset()
{
  for(int i=0; i<12; i++)
    w[i] = 0.0;
}

