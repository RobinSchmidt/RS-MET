  
rsVectorDbl dummy(const rsVectorDbl &v)
{
  rsVectorDbl r(v.dim);
  for(int n = 0; n < r.dim; n++)
    r[n] = 2.0*v[n];
  return r;
}

bool testVector()
{
  std::string testName = "rsVector";
  bool testResult = true;


  double a1[3] = {2.0,  3.0,  5.0};
  double a2[3] = {7.0, 11.0, 13.0};

  rsVectorDbl v1(3, a1);
  rsVectorDbl v2(3, a2);

  rsVectorDbl r; 

  // operators that take a scalar and a vector:
  r = 1.0 + v1; testResult &= r[0] ==  3.0; testResult &= r[1] ==  4.0;     testResult &= r[2] ==  6.0;
  r = v1 + 1.0; testResult &= r[0] ==  3.0; testResult &= r[1] ==  4.0;     testResult &= r[2] ==  6.0;
  r = 1.0 - v1; testResult &= r[0] == -1.0; testResult &= r[1] == -2.0;     testResult &= r[2] == -4.0;
  r = v1 - 1.0; testResult &= r[0] ==  1.0; testResult &= r[1] ==  2.0;     testResult &= r[2] ==  4.0;
  r = 2.0 * v1; testResult &= r[0] ==  4.0; testResult &= r[1] ==  6.0;     testResult &= r[2] == 10.0;
  r = v1 * 2.0; testResult &= r[0] ==  4.0; testResult &= r[1] ==  6.0;     testResult &= r[2] == 10.0;
  r = v1 / 2.0; testResult &= r[0] ==  1.0; testResult &= r[1] ==  1.5;     testResult &= r[2] == 2.5;
  r = 2.0 / v1; testResult &= r[0] ==  1.0; testResult &= r[1] ==  2.0/3.0; testResult &= r[2] == 2.0/5.0;

  // operators that take 2 vectors:
  r = v1 + v2; testResult &= r[0] ==  9.0; testResult &= r[1] == 14.0; testResult &= r[2] == 18.0;
  r = v2 - v1; testResult &= r[0] ==  5.0; testResult &= r[1] ==  8.0; testResult &= r[2] ==  8.0;

  // ==, !=, =,  +=, -=, *=, /=

  // element-wise application of functions to a vector:
  r = rsApplyFunction(v1, &exp);
  testResult &= r[0] == exp(2.0);
  testResult &= r[1] == exp(3.0);
  testResult &= r[2] == exp(5.0);

  r = rsApplyFunction(v1, 2.0, &pow);
  testResult &= r[0] == pow(2.0, 2.0);
  testResult &= r[1] == pow(3.0, 2.0);
  testResult &= r[2] == pow(5.0, 2.0);

  r = rsApplyFunction(2.0, v1, &pow);
  testResult &= r[0] == pow(2.0, 2.0);
  testResult &= r[1] == pow(2.0, 3.0);
  testResult &= r[2] == pow(2.0, 5.0);

  r = rsApplyFunction(v1, v2, &pow);
  testResult &= r[0] == pow(2.0,  7.0);
  testResult &= r[1] == pow(3.0, 11.0);
  testResult &= r[2] == pow(5.0, 13.0);

  //rsVectorDbl v3(1, a1);
  //r = dummy(v3);

  return testResult;
}





