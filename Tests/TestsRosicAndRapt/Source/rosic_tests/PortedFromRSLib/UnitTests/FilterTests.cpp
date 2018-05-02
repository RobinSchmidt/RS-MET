#include "FilterTests.h"

bool testFilterPolynomials(std::string &reportString)
{
  std::string testName = "FilterPolynomials";
  bool testResult = true;

  double a[30], at[30];
  for(int n = 1; n <= 10; n++)
  {
    halpernT2(at, n);
    rsHalpernPolynomial(a, n);
    testResult &= rsAreBuffersEqual(a, at, 2*n+1);

    papoulisL2(at, n);
    rsPapoulisPolynomial(a, n);
    testResult &= rsAreBuffersEqual(a, at, 2*n+1);
  }

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testHighOrderFilter(std::string &reportString)
{
  std::string testName = "HighOrderFilter";
  bool testResult = true;


  static const int N = 6;    // prototype filter order
  double fs    = 44100.0;  // samplerate
  double fc    =  1000.0;  // cutoff frequency
  double Ap    =     1.0;  // passband ripple in dB
  double As    =    50.0;  // stopband rejection in dB
  double bw    =     1.0;  // bandwidth in octaves
  int mode   = rsInfiniteImpulseResponseDesigner::BANDPASS;
  int method = rsPrototypeDesigner::ELLIPTIC;

  double tol = 1.e-10;

  // create and set up the filter designer object:
  rsInfiniteImpulseResponseDesigner designer;
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
  zt[1]  = zt[0].getConjugate();
  zt[2]  = rsComplexDbl(0.99561798226291653, -0.093513813924568437);
  zt[3]  = zt[2].getConjugate();
  zt[4]  = rsComplexDbl(0.97224253719572762,  0.23397531679049874);
  zt[5]  = zt[4].getConjugate();
  zt[6]  = rsComplexDbl(0.99631586399951144, -0.085759542576362605);
  zt[7]  = zt[6].getConjugate();
  zt[8]  = rsComplexDbl(0.91804188410822718,  0.39648341582343183);
  zt[9]  = zt[8].getConjugate();
  zt[10] = rsComplexDbl(0.99878505098211623, -0.049279021242830033);
  zt[11] = zt[10].getConjugate();

  pt[0]  = rsComplexDbl(0.99381379769518985, -0.10051120051999071);
  pt[1]  = pt[0].getConjugate();
  pt[2]  = rsComplexDbl(0.97762023598210512,  0.19956873870621017);
  pt[3]  = pt[2].getConjugate();
  pt[4]  = rsComplexDbl(0.98945304748144991, -0.10564891612703108);
  pt[5]  = pt[4].getConjugate();
  pt[6]  = rsComplexDbl(0.97331129154065377,  0.18766906880048043);
  pt[7]  = pt[6].getConjugate();
  pt[8]  = rsComplexDbl(0.97964553629083029, -0.12289201190055950);
  pt[9]  = pt[8].getConjugate(); 
  pt[10] = rsComplexDbl(0.97095305270689836,  0.15780232434708155);
  pt[11] = pt[10].getConjugate();

  for(int i = 0; i < 2*N; i++)
  {
    testResult &= rsIsCloseTo(p[i], pt[i], tol);
    testResult &= rsIsCloseTo(z[i], zt[i], tol);
  }

  // the zeros are in the wrong order
  // rsPrototypeDesigner::getPolesAndZeros seems OK

  //testResult &= z[0] == 


  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}



