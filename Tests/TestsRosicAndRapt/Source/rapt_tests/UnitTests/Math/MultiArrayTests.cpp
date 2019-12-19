
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

bool testMultiArray()
{
  std::string testName = "rsMultiArrayOld";
  std::string dummy;
  bool testResult = true;

  testResult &= testMultiArrayInit(dummy);
  testResult &= testMultiArrayOffsetComputation(dummy);
  testResult &= testMultiArrayOuterProduct(dummy);
  testResult &= testMultiArrayContraction(dummy);

  //appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}