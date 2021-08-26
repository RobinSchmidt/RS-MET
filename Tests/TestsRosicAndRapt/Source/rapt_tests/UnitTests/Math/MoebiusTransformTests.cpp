
bool testMoebiusTransform()
{
  // This test fails on mac. I guess that's due to the exact float comparisons - we may need a
  // tolerance. But before implementing that, it may make sense to switch to using rsComplex
  // instead of std::complex.
  
  std::string testName = "rsMoebiusTransform";
  bool testResult = true;

  //typedef std::complex<double> rsComplex;

  std::complex<double> i(0, 1);     // imaginary unit
  std::complex<double> a1, b1, c1, d1;
  std::complex<double> z, w;
  std::complex<double> det, fp1, fp2;

  // test the rsSqrtC function (move this later to testComplexFunctions):
  testResult &= (sqrt(+3.0 + 4.0*i) == (2.0 + 1.0*i));
  testResult &= (sqrt(+3.0 - 4.0*i) == (2.0 - 1.0*i));
  testResult &= (sqrt(-3.0 + 4.0*i) == (1.0 + 2.0*i));
  testResult &= (sqrt(-3.0 - 4.0*i) == (1.0 - 2.0*i));

  testResult &= (sqrt(+8.0 + 6.0*i) == (3.0 + 1.0*i));
  testResult &= (sqrt(+8.0 - 6.0*i) == (3.0 - 1.0*i));
  testResult &= (sqrt(-8.0 + 6.0*i) == (1.0 + 3.0*i));
  testResult &= (sqrt(-8.0 - 6.0*i) == (1.0 - 3.0*i));

  //testResult &= (rsSqrtC( 0.0 + 2.0*i) == (1.0 + 1.0*i)); // numerically close, but not exactly
  //w = rsSqrtC( 0.0 + 2.0*i); 

  a1 = 1.0 + 2.0*i;
  b1 = 3.0 + 4.0*i;
  c1 = 5.0 + 6.0*i;
  d1 = 7.0 + 8.0*i;
  rsMoebiusTransform<double> M1(a1, b1, c1, d1);

  z = 3.0 + 4.0*i;
  w = M1.getMappedNumber(z);

  // test inversion (the inverse of M1 applied to w = M1(z) should recover z):
  rsMoebiusTransform<double> invM1 = M1.getInverse();
  w = invM1.getMappedNumber(w);
  testResult &= rsIsCloseTo(abs(z-w), 0., 1.e-14);

  // test composition (M1 composed with its inverse should give the identity):
  rsMoebiusTransform<double> idM = M1.followedBy(invM1);
  testResult &= idM.isIdentity();

  det = idM.getDeterminant();

  // test normalize:
  idM.normalize();
  testResult &= idM.getDeterminant() == 1.0;

  // test fixpoint computation:
  M1.getFixPoints(fp1, fp2);
  w = M1.getMappedNumber(fp1);
  testResult &= rsIsCloseTo( abs(fp1-w), 0., 1.e-14);
  w = M1.getMappedNumber(fp2);
  testResult &= rsIsCloseTo( abs(fp2-w), 0., 1.e-14);

  /*
  M1.normalize();
  det = M1.getDeterminant();
  */

  // test construction from 3-points z1, z2, z3 and their images w1, w2, w3:
  std::complex<double> z1, z2, z3, w1, w2, w3;
  z1 = 1.0 + 2.0*i;
  z2 = 2.0 - 3.0*i;
  z3 = 3.0 + 4.0*i;
  w1 = 4.0 - 5.0*i;
  w2 = 5.0 + 6.0*i;
  w3 = 6.0 + 7.0*i;
  rsMoebiusTransform<double> M2 = 
    rsMoebiusTransform<double>::from3PointsAndImages(z1, z2, z3, w1, w2, w3);

  w = M2.getMappedNumber(z1); testResult &= (w == w1); 
  w = M2.getMappedNumber(z2); testResult &= (w == w2); 
  w = M2.getMappedNumber(z3); testResult &= (w == w3); 

  //....stuff to come

  //appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}
