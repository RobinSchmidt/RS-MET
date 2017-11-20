#include "ComplexTests.h"

bool testComplex(std::string &reportString)
{
  std::string testName = "rsComplex";
  bool testResult = true;

  /*
  // note: we first assign the sizes to variables (instead of driectly using sizeof inside the comparison) to be able to inspect the actual
  // sizes in the debugger when something goes wrong

  */


  // test, if the sizes are right:
  int complexDoubleSize = sizeof( rsComplex<double> );
  testResult &= ( complexDoubleSize  ==  16 );

  int complexFloatSize = sizeof( rsComplex<float> );
  testResult &= ( complexFloatSize  ==  8 );

  // test, if the memory layout is right:
  rsComplex<double> zd(3.0, 4.0);    // a Complex<double> variable
  double *pzd = (double*) &zd;       // a pointer to the first data-field (casted to double)
  double zdre = *pzd;                // the content of the pointer should be the real part
  double zdim = *(pzd+1);            // the imaginary part should be one size-of-double further
  testResult &= ( zdre == zd.re );
  testResult &= ( zdim == zd.im );

  rsComplex<float> zf(3.0, 4.0);
  float  *pzf = (float*) &zf;
  float  zfre = *pzf;
  float  zfim = *(pzf+1);
  testResult &= ( zfre == zf.re );
  testResult &= ( zfim == zf.im );

  //double r = zd.getRadius();
  //Complex<double> z2(4.0, 3.0);
  //double r = z1.getRadius();
  // ...move to other function ...and the content of this function into testComplexMemoryLayout


  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}





