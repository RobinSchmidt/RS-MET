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

	bool r = true;      // test result
  //float tol = 0.f;  // zero tolerance - float comparisons are currently exact
  float tol = 1.e-7f;
  CF p[3], z[3];      // arrays for retrieved poles
  PD pd;              // prototype designer object

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

bool filterSpecUnitTest()
{
  bool r = true;      // test result

  typedef RAPT::rsFilterSpecificationBA<double>  BA;
  typedef RAPT::rsFilterSpecificationZPK<double> ZPK;
  typedef std::complex<double> Complex;

  // create example ZPK-filter with 3 zeros, 2 poles and gain = 3
  Complex q1(-1.0, +2.0), q2(-1.0, -2.0), q3(-4.0, 0.0);
  Complex p1(-0.5, +0.5), p2(-0.5, -0.5);
  Complex k = 4.0;
  double inf = RS_INF(double);
  ZPK zpk32({ q1, q2, q3 }, { p1, p2}, k, inf); // sampleRate = inf -> analog filter

  // Analog case:
  //             (s-q1)*(s-q2)*(s-q3)     B0 + B1*s + B2*s^2 + B3*s^3
  // H(s) = k * ---------------------- = -----------------------------
  //             (s-p1)*(s-p2)            A0 + A1*s + A2*s^2
  //
  // multiplying out the zpk representation gives:
  Complex B0 = -k*q1*q2*q3, B1 = k*(q1*q2+q1*q3+q2*q3), B2 = -k*(q1+q2+q3), B3 = k;
  Complex A0 = p1*p2, A1 = -(p1+p2), A2 = 1;

  // we check now, if the built-in conversion-function gives the desired result:
  BA ba32 = zpk32.toBA();
  r &= ba32.b[0] == B0;
  r &= ba32.b[1] == B1;
  r &= ba32.b[2] == B2;
  r &= ba32.b[3] == B3;
  r &= ba32.a[0] == A0;
  r &= ba32.a[1] == A1;
  r &= ba32.a[2] == A2;

  // now, we convert back from ba to zpk and check, if we get our original zpk specifiction
  // properly reconstructed:
  ZPK zpkTmp = ba32.toZPK();
  r &= zpkTmp.equals(zpk32); 

  // Digital case:
  //             (1-q1/z)*(1-q2/z)*(1-q3/z)     b0 + b1/z + b2/z^2 + b3/z^3
  // H(z) = k * ---------------------------- = -----------------------------
  //             (1-p1/z)*(1-p2/z)              a0 + a1/z + a2/z^2
  // multiplying out the zpk representation gives:
  Complex b0 = k, b1 = -k*(q1+q2+q3), b2 = k*(q1*q2+q1*q3+q2*q3), b3 = -k*q1*q2*q3;
  Complex a0 = 1, a1 = -(p1+p2), a2 = p1*p2;
  // The coeffs are the same as in the analog case but in reverse order because in the digital 
  // domain, we multiply inverse powers of z (instead of regular powers of s in the analog domain)

  // We re-interpret the zpk32 specification as a digital one by setting the sample-rate to 
  // some finite value (1 in this case):
  zpk32.sampleRate = 1;
  // ...and now do the same tests as we did in the analog case:

  // ZPK -> BA:
  ba32 = zpk32.toBA();
  r &= ba32.b[0] == b0;
  r &= ba32.b[1] == b1;
  r &= ba32.b[2] == b2;
  r &= ba32.b[3] == b3;
  r &= ba32.a[0] == a0;
  r &= ba32.a[1] == a1;
  r &= ba32.a[2] == a2;

  // BA -> ZPK:
  zpkTmp = ba32.toZPK();
  r &= zpkTmp.equals(zpk32); 

  // test of conversions is done - now we evaluate the transfer-function at a couple of randomly
  // selected values for s or z an see, if both representations (ZPK and BA) give the same results:
  int numValues = 50;
  RAPT::rsNoiseGenerator<double> prng;
  prng.setRange(-2.0, +2.0);
  Complex z, s, H_zpk, H_ba, d;
  double tol = 1.e-10; // seems like we need a quite high tolerance - check for numeric issues
  int i;
  for(i = 0; i < numValues; i++) {
    Complex z = Complex(prng.getSample(), prng.getSample());
    H_zpk = zpk32.transferFunctionAt(z);
    H_ba  = ba32.transferFunctionAt( z);
    d     = H_zpk - H_ba;
    r &= abs(d) <= tol;
    //rsAssert(r);
  }




  return r;
}