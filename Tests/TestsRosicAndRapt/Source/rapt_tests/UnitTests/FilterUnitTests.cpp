#include "FilterUnitTests.h"

bool prototypeDesignUnitTest()
{
	bool r = true; // test result

  // shorthands:
  typedef rsPrototypeDesignerF PD;
  typedef complex<float> CF;
  float inf = std::numeric_limits<float>::infinity();

  CF p[5], z[5]; // arrays for retrieved poles
  PD pd;         // prototype designer object

  // 5th order Papoulis lowpass:
  pd.setApproximationMethod(PD::PAPOULIS);
  pd.setPrototypeMode(PD::LOWPASS_PROTOTYPE);
  pd.setOrder(5);
  pd.getPolesAndZeros(&p[0], &z[0]);                 // non-redundant upper halfplane poles/zeros
  r &= p[0] == CF(-0.153586745f, 0.968145967f);
  r &= p[1] == CF(-0.388139844f, 0.588632464f);
  r &= p[2] == CF(-0.468089849f, 0.f);
  r &= z[0] == CF(inf, 0.f);
  r &= z[1] == CF(inf, 0.f);
  r &= z[2] == CF(inf, 0.f);
  r &= pd.getNumFinitePoles() == 5;
  r &= pd.getNumFiniteZeros() == 0;




  //pd.setPrototypeMode(PD::LOWSHELV_PROTOTYPE);
  //pd.setGain(6.02f);
	
	return r;
}
