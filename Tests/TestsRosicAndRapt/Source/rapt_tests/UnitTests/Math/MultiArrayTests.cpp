
bool testMultiArrayInit(std::string &reportString)
{
  std::string testName = "rsMultiArrayInit";
  bool testResult = true;


  rsUint32 indices[1];
  indices[0] = 3;    // that value should be irrelevant
  float data[10];
  rsArray::fillWithRangeLinear(data, 10, 1.f, 10.f);

  // create array with 0 indices (a scalar) and see, if it initializes correctly:
  rsMultiArrayOld<float> a0;
  testResult &= a0.getNumElements()    == 1;
  testResult &= a0.getElement(indices) == 0.f;
  a0.setDataFromFlatArray(&data[5]);
  testResult &= a0.getElement(indices) == 6.f;


  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testMultiArrayOffsetComputation(std::string &reportString)
{
  std::string testName = "rsMultiArrayOffsetComputation";
  bool testResult = true;

  // create an 7x11x5 array (3 indices):
  rsUint32 indices[3];
  indices[0] = 7;
  indices[1] = 11;
  indices[2] = 5;
  rsMultiArrayOld<float> a3(3, indices);
  testResult &= a3.getNumElements() == 385;

  // test offset/index conversion with this array:
  int i, j, k;
  rsUint32 targetOffset = 0;
  for(i = 0; i < 7; i++)
  {
    indices[0] = i;
    for(j = 0; j < 11; j++)
    {
      indices[1] = j;
      for(k = 0; k < 5; k++)
      {
        indices[2] = k;

        // convert the set of indices to an offset:
        rsUint32 offset = a3.offsetFromIndices(indices);
        testResult &= (offset == targetOffset);
        targetOffset++;

        // convert back from the offset to indices:
        a3.indicesFromOffset(offset, indices);
        testResult &= ((int)indices[0] == i);
        testResult &= ((int)indices[1] == j);
        testResult &= ((int)indices[2] == k);
      }
    }
  }

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testMultiArrayOuterProduct(std::string &reportString)
{
  std::string testName = "rsMultiArrayOuterProduct";
  bool testResult = true;

  // create an 4x2x3 array (3 indices):
  rsUint32 indices1[3];
  indices1[0] = 4;
  indices1[1] = 2;
  indices1[2] = 3;

  rsUint32 indices2[2];
  indices2[0] = 6;
  indices2[1] = 5;

  rsUint32 indicesResult[5];

  float data[30];
  rsArray::fillWithRangeLinear(data, 30, 1.f, 30.f);
  data[23] = 32; // these 2 values will be used in the divisions, so we want powers of 2 to make
  data[29] = 64; // numerically exact division possible (for easier comparison witn target results)

  rsMultiArrayOld<float> a(3, indices1);
  rsMultiArrayOld<float> b(2, indices2);
  rsMultiArrayOld<float> c;
  a.setDataFromFlatArray(data);
  b.setDataFromFlatArray(data);

  rsMultiArrayOld<float>::outerProduct(a, b, c);

  int i, j, k, p, q;
  for(i = 0; i < 4; i++)
  {
    indicesResult[0] = i;
    indices1[0] = i;
    for(j = 0; j < 2; j++)
    {
      indicesResult[1] = j;
      indices1[1] = j;
      for(k = 0; k < 3; k++)
      {
        indicesResult[2] = k;
        indices1[2] = k;
        for(p = 0; p < 6; p++)
        {
          indicesResult[3] = p;
          indices2[0] = p;
          for(q = 0; q < 5; q++)
          {
            indicesResult[4] = q;
            indices2[1] = q;
            testResult &= (c.getElement(indicesResult)
              == a.getElement(indices1) * b.getElement(indices2));
          }
        }
      }
    }
  }

  // test tail-division (retrieve the left factor):
  indices1[0] = 4;
  indices1[1] = 2;
  indices1[2] = 3;
  rsMultiArrayOld<float> a2(3, indices1);
  rsMultiArrayOld<float>::leftFactor(c, b, a2);
  testResult &= (a == a2);

  // test head-division (retrieve the right factor):
  indices2[0] = 6;
  indices2[1] = 5;
  rsMultiArrayOld<float> b2(2, indices2);
  rsMultiArrayOld<float>::rightFactor(c, a, b2);
  testResult &= (b == b2);

  // write tests to verify that division-compatibility check works

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool testMultiArrayContraction(std::string &reportString)
{
  std::string testName = "rsMultiArrayContraction";
  bool testResult = true;

  // create an 2x3x4x3x5 array (5 indices):
  rsUint32 subjectIndices[5];
  subjectIndices[0] = 2;
  subjectIndices[1] = 3;  // summation index 1
  subjectIndices[2] = 4;
  subjectIndices[3] = 3;  // summation index 2
  subjectIndices[4] = 5;

  float data[360];
  rsArray::fillWithRangeLinear(data, 360, 1.f, 360.f);

  rsMultiArrayOld<float> a(5, subjectIndices);
  a.setDataFromFlatArray(data);

  // check contraction with respect to index-pair (1, 3):
  rsMultiArrayOld<float> a13 = rsMultiArrayOld<float>::contract(a, 1, 3);
  rsUint32 resultIndices[3];
  for(int i = 0; i < 2; i++)
  {
    resultIndices [0] = i;
    subjectIndices[0] = i;
    for(int j = 0; j < 4; j++)
    {
      resultIndices [1] = j;
      subjectIndices[2] = j;
      for(int k = 0; k < 5; k++)
      {
        resultIndices [2] = k;
        subjectIndices[4] = k;
        float sum = 0.f;
        for(int s = 0; s < 3; s++)
        {
          subjectIndices[1] = subjectIndices[3] = s;
          sum += a.getElement(subjectIndices);
        }
        testResult &= (sum == a13.getElement(resultIndices));
      }
    }
  }

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}


bool testMultiArray1()
{
  // maybe split the function int 1D,2D,3D tests

  bool r = true;

  typedef std::vector<int> VecI;
  typedef std::vector<float> VecF;
  typedef rsMultiArray<float> MA;


  // let's see, if the example code in the documentation of rsMultiArrayView works:
  float data[24];                             // flat C-array of data
  rsMultiArrayView<float> A({2,4,3}, data);   // we want to interpret the data as 2x4x3 3D array
  A(0,0,0) = 111.f;                           // set element at position 0,0,0 (the first one)
  A(1,3,2) = 243.f;                           // set element at position 1,3,2 (the last one)
  // ...yep - it does :-)





  // test memory layout:

  float* p;  // pointer to the data 

  // 3D vector:
  MA a1 = MA(VecI{3});
  a1(0) = 1; a1(1) = 2; a1(2) = 3;
  p = a1.getDataPointer();
  r &= p[0] == 1 && p[1] == 2 && p[2] == 3;

  //a1 = a1 + a1;


  // 3x2 matrix:
  MA a2 = MA(VecI{3,2});  
  a2(0,0) = 11; a2(0,1) = 12;
  a2(1,0) = 21; a2(1,1) = 22;
  a2(2,0) = 31; a2(2,1) = 32;
  p = a2.getDataPointer();
  r &= p[0] == 11 && p[1] == 12 && p[2] == 21 && p[3] == 22 && p[4] == 31 && p[5] == 32;

  // reshape into 2x3 matrix:
  a2.setShape(VecI({2,3}));
  r &= a2(0,0) == 11; r &= a2(0,1) == 12; r &= a2(0,2) == 21; 
  r &= a2(1,0) == 22; r &= a2(1,1) == 31; r &= a2(1,2) == 32;

  // 2x4x3 block/cuboid:
  MA a3 = MA(VecI{2,4,3}); 
  p = a3.getDataPointer();
  a3(0,0,0) = 111; a3(0,0,1) = 112; a3(0,0,2) = 113; r &= p[0] ==111 && p[1] ==112 && p[2] ==113;
  a3(0,1,0) = 121; a3(0,1,1) = 122; a3(0,1,2) = 123; r &= p[3] ==121 && p[4] ==122 && p[5] ==123;
  a3(0,2,0) = 131; a3(0,2,1) = 132; a3(0,2,2) = 133; r &= p[6] ==131 && p[7] ==132 && p[8] ==133;
  a3(0,3,0) = 141; a3(0,3,1) = 142; a3(0,3,2) = 143; r &= p[9] ==141 && p[10]==142 && p[11]==143;

  a3(1,0,0) = 211; a3(1,0,1) = 212; a3(1,0,2) = 213; r &= p[12]==211 && p[13]==212 && p[14]==213;
  a3(1,1,0) = 221; a3(1,1,1) = 222; a3(1,1,2) = 223; r &= p[15]==221 && p[16]==222 && p[17]==223;
  a3(1,2,0) = 231; a3(1,2,1) = 232; a3(1,2,2) = 233; r &= p[18]==231 && p[19]==232 && p[20]==233;
  a3(1,3,0) = 241; a3(1,3,1) = 242; a3(1,3,2) = 243; r &= p[21]==241 && p[22]==242 && p[23]==243;





  // allow the user to specify an allocator so we can unit-test the memory allocation avoidance
  // (copy elision for return values, for example) - explorin allocators should be done the research 
  // repo

  return r;
}

bool testMultiArray()
{
  std::string testName = "rsMultiArray (old and new)";
  std::string dummy;
  bool testResult = true;

  // these are the tests that work on the odl implementation -> adapt them!
  testResult &= testMultiArrayInit(dummy);
  testResult &= testMultiArrayOffsetComputation(dummy);
  testResult &= testMultiArrayOuterProduct(dummy);
  testResult &= testMultiArrayContraction(dummy);

  // tests with the new class:
  testResult &= testMultiArray1();



  //appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}