//#include "FilterTests.h"

typedef std::complex<double> rsComplexDbl;  // get rid

bool testFilterPolynomials()
{
  bool testResult = true;

  double a[30], at[30];
  for(int n = 1; n <= 10; n++)
  {
    halpernT2(at, n);
    rsPrototypeDesignerD::halpernPolynomial(a, n);
    testResult &= RAPT::rsArrayTools::equal(a, at, 2*n+1);

    papoulisL2(at, n);
    rsPrototypeDesignerD::papoulisPolynomial(a, n);
    testResult &= RAPT::rsArrayTools::equal(a, at, 2*n+1);
  }

  return testResult;
}

bool testHighOrderFilter1()
{
  bool testResult = true;

  static const int N = 6;    // prototype filter order
  double fs    = 44100.0;  // samplerate
  double fc    =  1000.0;  // cutoff frequency
  double Ap    =     1.0;  // passband ripple in dB
  double As    =    50.0;  // stopband rejection in dB
  double bw    =     1.0;  // bandwidth in octaves
  int mode   = rsInfiniteImpulseResponseDesignerD::BANDPASS;
  int method = rsPrototypeDesignerD::ELLIPTIC;

  double tol = 1.e-10;

  // create and set up the filter designer object:
  rsInfiniteImpulseResponseDesignerD designer;
  designer.setPrototypeOrder(N);
  designer.setSampleRate(fs);
  designer.setFrequency(fc);
  designer.setRipple(Ap);
  designer.setStopbandRejection(As);
  designer.setBandwidth(bw);
  designer.setMode(mode);
  designer.setApproximationMethod(method);


  // compute poles and zeros:
  rsComplexDbl z[2*N], p[2*N];        // arrays of poles and zeros
  designer.getPolesAndZeros(p, z);

  // target polse and zeros:
  rsComplexDbl zt[2*N], pt[2*N];
  zt[0]  = rsComplexDbl(0.97661958515038094,  0.21497484946081000);
  zt[1]  = conj(zt[0]);
  zt[2]  = rsComplexDbl(0.99561798226291653, -0.093513813924568437);
  zt[3]  = conj(zt[2]);
  zt[4]  = rsComplexDbl(0.97224253719572762,  0.23397531679049874);
  zt[5]  = conj(zt[4]);
  zt[6]  = rsComplexDbl(0.99631586399951144, -0.085759542576362605);
  zt[7]  = conj(zt[6]);
  zt[8]  = rsComplexDbl(0.91804188410822718,  0.39648341582343183);
  zt[9]  = conj(zt[8]);
  zt[10] = rsComplexDbl(0.99878505098211623, -0.049279021242830033);
  zt[11] = conj(zt[10]);

  pt[0]  = rsComplexDbl(0.99381379769518985, -0.10051120051999071);
  pt[1]  = conj(pt[0]);
  pt[2]  = rsComplexDbl(0.97762023598210512,  0.19956873870621017);
  pt[3]  = conj(pt[2]);
  pt[4]  = rsComplexDbl(0.98945304748144991, -0.10564891612703108);
  pt[5]  = conj(pt[4]);
  pt[6]  = rsComplexDbl(0.97331129154065377,  0.18766906880048043);
  pt[7]  = conj(pt[6]);
  pt[8]  = rsComplexDbl(0.97964553629083029, -0.12289201190055950);
  pt[9]  = conj(pt[8]); 
  pt[10] = rsComplexDbl(0.97095305270689836,  0.15780232434708155);
  pt[11] = conj(pt[10]);

  for(int i = 0; i < 2*N; i++)
  {
    testResult &= abs(p[i]-pt[i]) <= tol;
    testResult &= abs(z[i]-zt[i]) <= tol;
    //testResult &= rsIsCloseTo(p[i], pt[i], tol);
    //testResult &= rsIsCloseTo(z[i], zt[i], tol);
  }

  // the zeros are in the wrong order
  // rsPrototypeDesigner::getPolesAndZeros seems OK

  //testResult &= z[0] == 

  return testResult;
}

// move to test utilities or something:
bool equal(double x, double y, double tol)
{
  if(x == RS_INF(double) && y == RS_INF(double))
    return true;
  return abs(x-y) <= tol;
}
bool equal(rosic::Complex z1, std::complex<double> z2, double tol)
{
  return equal(z1.re, z2.real(), tol) && equal(z1.im, z2.imag(), tol); 
}
bool equal(const std::vector<rosic::Complex>& array1,
  const std::vector<std::complex<double>>& array2)
{
  if(array1.size() != array2.size())
    return false;
  double tol = 1.e-12;
  for(size_t i = 0; i < array1.size(); i++) {
    if(!equal(array1[i], array2[i], tol))
      return false; }
  return true;
}

bool testIIRDesign(int method, int mode, int protoOrder)
{
  // method: 0..5, mode: 0..7, order: 1..20
  // Compare results from the old rosic and new rapt implementation
  // ...obsolete now ...delete soon

  bool r = true;

  double fs = 44100;  // sample rate
  double fc = 1000;   // cutoff/center frequency
  double bw = 2;      // bandwidth in octaves
  double g  = 6;      // peak/shelf gain
  double rp = 3;      // passband ripple in dB
  double rj = 60;     // stopband rejection in dB
  size_t N  = protoOrder;

  // typedefs for convenience:
  typedef rosic::rsPrototypeDesignerD PTD1;
  typedef rsPrototypeDesignerD   PTD2;
  //typedef rosic::rsPoleZeroMapper PZM1;
  //typedef rsPoleZeroMapperD   PZM2;
  //typedef rosic::rsFilterCoefficientConverter FCC1;
  //typedef rsFilterCoefficientConverterD FCC2;
  typedef rosic::rsInfiniteImpulseResponseDesignerD IIRD1;
  typedef rsInfiniteImpulseResponseDesignerD   IIRD2;

  // maybe check here, if the enums in PTD1/PTD2, etc. match

  // create prototype poles/zeros:
  std::vector<std::complex<double>> pp1(N), pz1(N);
  std::vector<std::complex<double>> pp2(N), pz2(N);
  int protoMode = PTD1::LOWPASS_PROTOTYPE;
  if(mode == IIRD1::LOW_SHELV || mode == IIRD1::HIGH_SHELV || mode == IIRD1::PEAK )
    protoMode = PTD1::LOWSHELV_PROTOTYPE;
  PTD1 ptd1;
  PTD2 ptd2;
  ptd1.setPrototypeMode(protoMode);
  ptd2.setPrototypeMode(protoMode);
  ptd1.setApproximationMethod(method);
  ptd2.setApproximationMethod(method);
  ptd1.setPassbandRipple(rp);
  ptd2.setPassbandRipple(rp);
  ptd1.setStopbandRejection(rj);
  ptd2.setStopbandRejection(rj);
  ptd1.setGain(g);
  ptd2.setGain(g);
  ptd1.setOrder(protoOrder);
  ptd2.setOrder(protoOrder);
  ptd1.getPolesAndZeros(&pp1[0], &pz1[0]);
  ptd2.getPolesAndZeros(&pp2[0], &pz2[0]);
  //r &= equal(pp1, pp2);
  //r &= equal(pz1, pz2);

  // s-domain frequency transform:
  //size_t N2 = N;
  //if(mode == IIRD1::BANDPASS || mode == IIRD1::BANDREJECT || mode == IIRD1::PEAK )
  //  N2 *= 2; // order doubling modes

  // s-to-z bilinear transform:


  // convert to biquad coeffs:


  // the whole design encapsulated:
  IIRD1 iird1;
  IIRD2 iird2;
  iird1.setMode(mode);
  iird2.setMode(mode);
  iird1.setApproximationMethod(method);
  iird2.setApproximationMethod(method);
  iird1.setRipple(rp);  // rename to setPassbandRipple
  iird2.setRipple(rp);
  iird1.setStopbandRejection(rj);
  iird2.setStopbandRejection(rj);
  iird1.setGain(g);
  iird2.setGain(g);
  iird1.setPrototypeOrder(protoOrder);
  iird2.setPrototypeOrder(protoOrder);
  iird1.setSampleRate(fs);
  iird2.setSampleRate(fs);
  iird1.setFrequency(fc);
  iird2.setFrequency(fc);
  iird1.setBandwidth(bw);
  iird2.setBandwidth(bw);
  size_t N2 = (size_t) iird1.getFinalFilterOrder();
  //std::vector<rosic::Complex> p1(N2), z1(N2);
  std::vector<std::complex<double>> p1(N2), z1(N2);
  std::vector<std::complex<double>> p2(N2), z2(N2);
  iird1.getPolesAndZeros(&p1[0], &z1[0]);
  iird2.getPolesAndZeros(&p2[0], &z2[0]);
  //r &= equal(p1, p2);
  //r &= equal(z1, z2);

  return r;
}

bool testHighOrderFilter()
{
  bool r = true;

  //r &= testHighOrderFilter1();
  //r &= testIIRDesign(1, 1, 2);
  //r &= testIIRDesign(0, 3, 1);  // 1st order butterworth highpass -> fails because the old one
  r &= testIIRDesign(4, 4, 4); 

  for(int method = 0; method <= 5; method++)
    for(int mode = 0; mode <= 7; mode++)
      for(int order = 1; order <= 20; order++)
      {
        r &= testIIRDesign(method, mode, order);
        rsAssert(r == true);
      }

  return r;
}





