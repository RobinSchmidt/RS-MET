#include "FilterUnitTests.h"

// maybe move somewhere else for sharing
bool isCloseTo(complex<float> x, complex<float> y, float tol)
{
  return abs(x-y) <= tol;
}

bool prototypeDesignUnitTest()
{
  // shorthands:
  typedef rsPrototypeDesignerF PD;
  typedef complex<float> CF;
  float inf = std::numeric_limits<float>::infinity();

	bool r = true;   // test result
  float tol = 0.f; // zero tolerance - float comparisons are currently exact
  CF p[3], z[3];   // arrays for retrieved poles
  PD pd;           // prototype designer object

  // Papoulis filters:
  pd.setApproximationMethod(PD::PAPOULIS);
  pd.setPrototypeMode(PD::LOWPASS_PROTOTYPE);
  pd.setOrder(5);                                               // 5th order lowpass
  pd.getPolesAndZeros(&p[0], &z[0]);              
  r &= isCloseTo(p[0], CF(-0.153586745f, 0.968145967f), tol);
  r &= isCloseTo(p[1], CF(-0.388139844f, 0.588632464f), tol);
  r &= isCloseTo(p[2], CF(-0.468089849f, 0.f),          tol);
  r &= z[0] == CF(inf, 0.f);
  r &= z[1] == CF(inf, 0.f);
  r &= z[2] == CF(inf, 0.f);
  r &= pd.getNumFinitePoles() == 5;
  r &= pd.getNumFiniteZeros() == 0;
  pd.setPrototypeMode(PD::LOWSHELV_PROTOTYPE);                  // 5th order low-boost
  pd.setGain(+6);
  pd.getPolesAndZeros(&p[0], &z[0]);              
  r &= isCloseTo(p[0], CF(-0.147757247f, 0.931399286f), tol);
  r &= isCloseTo(p[1], CF(-0.373407722f, 0.566290498f), tol);
  r &= isCloseTo(p[2], CF(-0.450323164f, 0.f),          tol);
  r &= isCloseTo(z[0], CF(-0.185926169f, 0.996051848f), tol);
  r &= isCloseTo(z[1], CF(-0.481426001f, 0.611809969f), tol);
  r &= isCloseTo(z[2], CF(-0.590867400f, 0.f),          tol);
  r &= pd.getNumFinitePoles() == 5;
  r &= pd.getNumFiniteZeros() == 5;
  pd.setGain(-6);                                               // 5th order low-cut
  pd.getPolesAndZeros(&p[0], &z[0]); 
  r &= isCloseTo(p[0], CF(-0.185926169f, 0.996051848f), tol);
  r &= isCloseTo(p[1], CF(-0.481426001f, 0.611809969f), tol);
  r &= isCloseTo(p[2], CF(-0.590867400f, 0.f),          tol);
  r &= isCloseTo(z[0], CF(-0.147757247f, 0.931399286f), tol);
  r &= isCloseTo(z[1], CF(-0.373407722f, 0.566290498f), tol);
  r &= isCloseTo(z[2], CF(-0.450323164f, 0.f),          tol);
  r &= pd.getNumFinitePoles() == 5;
  r &= pd.getNumFiniteZeros() == 5;
  pd.setPrototypeMode(PD::LOWPASS_PROTOTYPE);
  pd.setOrder(6);                                               // 6th order lowpass
  pd.getPolesAndZeros(&p[0], &z[0]);
  r &= isCloseTo(p[0], CF(-0.115192428f, 0.977922320f), tol);
  r &= isCloseTo(p[1], CF(-0.308961034f, 0.698167443f), tol);
  r &= isCloseTo(p[2], CF(-0.438901573f, 0.239981338f), tol);
  r &= z[0] == CF(inf, 0.f);
  r &= z[1] == CF(inf, 0.f);
  r &= z[2] == CF(inf, 0.f);
  r &= pd.getNumFinitePoles() == 6;
  r &= pd.getNumFiniteZeros() == 0;
  pd.setPrototypeMode(PD::LOWSHELV_PROTOTYPE);                  // 6th order low-boost
  pd.setGain(+6);
  pd.getPolesAndZeros(&p[0], &z[0]);  
  r &= isCloseTo(p[0], CF(-0.111895919f, 0.949936688f), tol);
  r &= isCloseTo(p[1], CF(-0.300119370f, 0.678187668f), tol);
  r &= isCloseTo(p[2], CF(-0.426341325f, 0.233113691f), tol);
  r &= isCloseTo(z[0], CF(-0.139246285f, 0.999132454f), tol);
  r &= isCloseTo(z[1], CF(-0.380611271f, 0.720496416f), tol);
  r &= isCloseTo(z[2], CF(-0.534630060f, 0.254959404f), tol);

  // 6th order Bessel low-cut:
  pd.setApproximationMethod(PD::BESSEL);
  pd.setPrototypeMode(PD::LOWSHELV_PROTOTYPE);
  pd.setGain(-6);
  pd.getPolesAndZeros(&p[0], &z[0]);
  r &= isCloseTo(p[0], CF(-0.659143269f, 1.40492857f),  tol);
  r &= isCloseTo(p[1], CF(-1.15471625f,  1.03821945f),  tol);
  r &= isCloseTo(p[2], CF(-1.53899157f,  0.387341291f), tol);
  r &= isCloseTo(z[0], CF(-0.750605762f, 1.34034932f),  tol);
  r &= isCloseTo(z[1], CF(-1.11451340f,  0.783524334f), tol);
  r &= isCloseTo(z[2], CF(-1.26746202f,  0.258809835f), tol);
	
	return r;
}
