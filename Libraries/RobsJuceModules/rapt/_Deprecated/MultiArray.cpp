// construction/destruction:

template<class ElementType>
rsMultiArrayOld<ElementType>::rsMultiArrayOld()
{
  numIndices  = 0;
  indexRanges = NULL;
  allocateDataArray(true);
}

template<class ElementType>
rsMultiArrayOld<ElementType>::rsMultiArrayOld(const rsUint32 numIndices,
                                        const rsUint32 indexRanges[],
                                        bool initWithZeros)
{
  this->indexRanges = NULL;
  allocateAndInitIndexArray(numIndices, indexRanges);
  allocateDataArray(initWithZeros);
}

template<class ElementType>
rsMultiArrayOld<ElementType>::rsMultiArrayOld(const rsMultiArrayOld& other)
{
  allocateAndInitIndexArray(other.numIndices, other.indexRanges);
  allocateDataArray(false);
  memcpy(data, other.data, getNumElements()*sizeof(ElementType));
}

template<class ElementType>
rsMultiArrayOld<ElementType>::~rsMultiArrayOld()
{
  delete[] indexRanges;
  delete[] data;
}

// operators:

template<class ElementType>
bool rsMultiArrayOld<ElementType>::operator==(const rsMultiArrayOld<ElementType> &other) const
{
  if( !this->isOfSameTypeAs(other) )
    return false;
  else
    return rsArrayTools::equal(data, other.data, getNumElements());
}

template<class ElementType>
bool rsMultiArrayOld<ElementType>::operator!=(const rsMultiArrayOld<ElementType> &other) const
{
  return !(*this == other);
}

/** Adds a scalar and a multi-array. */
/*
template<class ElementType>
RS_INLINE rsMultiArray<ElementType> operator+(const ElementType &x, const rsMultiArray<ElementType> &v)
{
  rsMultiArray<ElementType> result(v.dim);
  for(int i=0; i<v.dim; i++)
    result.v[i] = v.v[i] + x;
  return result;
}
*/

/** Subtracts a multi-array from a scalar. */
/*
template<class ElementType>
RS_INLINE rsMultiArray<ElementType> operator-(const ElementType &x, const rsMultiArray<ElementType> &v)
{
  rsMultiArray<ElementType> result(v.dim);
  for(int i=0; i<v.dim; i++)
    result.v[i] = x - v.v[i];
  return result;
}
*/

/** Multiplies a scalar and a multi-array. */
/*
template<class ElementType>
RS_INLINE rsMultiArray<ElementType> operator*(const ElementType &x, const rsMultiArray<ElementType> &v)
{
  rsMultiArray<ElementType> result(v.dim);
  for(int i=0; i<v.dim; i++)
    result.v[i] = x * v.v[i];
  return result;
}
*/

// manipulations:

template<class ElementType>
void rsMultiArrayOld<ElementType>::setSize(const rsUint32 newNumIndices, const rsUint32 newIndexRanges[])
{
  if( getNumElements() != rsArrayTools::product(newIndexRanges, newNumIndices) )
  {
    freeIndexArray();
    freeDataArray();
    allocateAndInitIndexArray(newNumIndices, newIndexRanges);
    allocateDataArray(false); // maybe have a boolean parameter to pass through
  }
  else
  {
    freeIndexArray();
    allocateAndInitIndexArray(newNumIndices, newIndexRanges);
  }
}

template<class ElementType>
void rsMultiArrayOld<ElementType>::applyFunction(ElementType (*f) (ElementType))
{
  for(rsUint32 i = 0; i < getNumElements(); i++)
    data[i] = f(data[i]);
}

template<class ElementType>
void rsMultiArrayOld<ElementType>::fillWithValue(const ElementType &value)
{
  for(rsUint32 i = 0; i < getNumElements(); i++)
    data[i] = value;
}

template<class ElementType>
void rsMultiArrayOld<ElementType>::setElement(const rsUint32 indices[], const ElementType &value)
{
  data[offsetFromIndices(indices)] = value;
}

template<class ElementType>
void rsMultiArrayOld<ElementType>::setDataFromFlatArray(const ElementType flatArray[])
{
  memcpy(data, flatArray, getNumElements()*sizeof(ElementType));
}

// inquiry:

template<class ElementType>
rsUint32 rsMultiArrayOld<ElementType>::offsetFromIndices(const rsUint32 indices[]) const
{
  rsUint32 offset = 0;
  rsUint32 scaler = 1;
  for(rsInt32 i = numIndices-1; i >= 0; i--)
  {
    offset += scaler * indices[i];
    scaler *= indexRanges[i];
  }
  return offset;
}

template<class ElementType>
void rsMultiArrayOld<ElementType>::indicesFromOffset(rsUint32 offset, rsUint32 indices[]) const
{
  rsUint32 divider = getNumElements();
  for(rsUint32 k = 0; k < numIndices; k++)
  {
    divider    /= indexRanges[k];
    indices[k]  = offset / divider;
    offset     -= divider * indices[k];
  }
  // maybe a more efficient algorithm can be found (using one mod and one mul instead of the 2
  // divs)
}

template<class ElementType>
ElementType rsMultiArrayOld<ElementType>::getElement(const rsUint32 indices[]) const
{
  return data[offsetFromIndices(indices)];
}

template<class ElementType>
rsUint32 rsMultiArrayOld<ElementType>::getNumElements() const
{
  return rsArrayTools::product(indexRanges, numIndices);
}

template<class ElementType>
void rsMultiArrayOld<ElementType>::getDataAsFlatArray(ElementType flatArray[]) const
{
  memcpy(flatArray, data, getNumElements()*sizeof(ElementType));
}

template<class ElementType>
bool rsMultiArrayOld<ElementType>::isOfSameTypeAs(const rsMultiArrayOld<ElementType> &other) const
{
  if( numIndices != other.numIndices )
    return false;
  else
    return rsArrayTools::equal(indexRanges, other.indexRanges, numIndices);
}

// internal functions:

template<class ElementType>
void rsMultiArrayOld<ElementType>::allocateAndInitIndexArray(const rsUint32 newNumIndices,
                                                          const rsUint32 newIndexRanges[])
{
  numIndices = newNumIndices;
  if( numIndices > 0 )
  {
    indexRanges = new rsUint32[numIndices];
    memcpy(indexRanges, newIndexRanges, numIndices*sizeof(rsUint32));
  }
}

template<class ElementType>
void rsMultiArrayOld<ElementType>::allocateDataArray(bool initWithZeros)
{
  data = new ElementType[getNumElements()];
  if( initWithZeros )
    fillWithValue(ElementType(0));
}

template<class ElementType>
void rsMultiArrayOld<ElementType>::freeIndexArray()
{
  delete[] indexRanges;
  indexRanges = NULL;
  numIndices = 0;
}

template<class ElementType>
void rsMultiArrayOld<ElementType>::freeDataArray()
{
  delete[] data;
  data = NULL;
}

template<class ElementType>
void rsMultiArrayOld<ElementType>::add(const rsMultiArrayOld<ElementType> &left,
                                    const rsMultiArrayOld<ElementType> &right,
                                    rsMultiArrayOld<ElementType> &result)
{
  rsAssert(left.isOfSameTypeAs(right));
  result.setSize(left.numIndices, left.indexRanges);
  rsArrayTools::add(left.data, right.data, result.data, left.getNumElements());
}

template<class ElementType>
void rsMultiArrayOld<ElementType>::subtract(const rsMultiArrayOld<ElementType> &left,
                                         const rsMultiArrayOld<ElementType> &right,
                                         rsMultiArrayOld<ElementType> &result)
{
  rsAssert(left.isOfSameTypeAs(right));
  result.setSize(left.numIndices, left.indexRanges);
  rsArrayTools::subtract(left.data, right.data, result.data, left.getNumElements());
}

template<class ElementType>
void rsMultiArrayOld<ElementType>::outerProduct(const rsMultiArrayOld<ElementType> &left,
                                             const rsMultiArrayOld<ElementType> &right,
                                             rsMultiArrayOld<ElementType> &result)
{
  rsUint32 *resultIndices = new rsUint32[left.numIndices + right.numIndices];
  memcpy( resultIndices,                  left.indexRanges,  left.numIndices*sizeof(rsUint32));
  memcpy(&resultIndices[left.numIndices], right.indexRanges, right.numIndices*sizeof(rsUint32));
  result.setSize(left.numIndices + right.numIndices, resultIndices);
  for(rsUint32 i = 0; i < result.getNumElements(); i++)
  {
    result.indicesFromOffset(i, resultIndices);
    result.data[i] = left.getElement(resultIndices)
      * right.getElement(&resultIndices[left.numIndices]);
  }
  delete[] resultIndices;
  // maybe a more efficient algorithm can be derived that avoids the repeated back-and-forth
  // conversions between flat and hierarchical indices
}

template<class ElementType>
void rsMultiArrayOld<ElementType>::leftFactor(const rsMultiArrayOld<ElementType> &outerProduct,
                                           const rsMultiArrayOld<ElementType> &rightFactor,
                                           rsMultiArrayOld<ElementType> &result)
{
  rsUint32 numResultIndices = outerProduct.numIndices-rightFactor.numIndices;
  rsUint32 *resultIndices   = new rsUint32[numResultIndices];

  rsAssert(rsArrayTools::equal(&outerProduct.indexRanges[numResultIndices],
    rightFactor.indexRanges, rightFactor.numIndices)); // multiarrays incompatible

  memcpy(resultIndices, outerProduct.indexRanges, numResultIndices*sizeof(rsUint32));
  result.setSize(numResultIndices, resultIndices);

  // For each element in the result, we have several (equal to the number of elements in
  // rightFactor) options how to compute it. We should use the one which gives the numerically
  // most precise result - for the time being, we pick the element with maximum absolute
  // value (maybe write function getBestDivisorIndex):
  rsUint32      stepSize = rightFactor.getNumElements();
  rsUint32      offset   = (rsUint32) rsArrayTools::maxAbsIndex(rightFactor.data, (int) stepSize);
  ElementType scaler   = ElementType(1) / rightFactor.data[offset];

  for(rsUint32 i = 0; i < result.getNumElements(); i++)
  {
    result.data[i] = scaler * outerProduct.data[offset];
    offset += stepSize;
  }

  delete[] resultIndices;
}

template<class ElementType>
void rsMultiArrayOld<ElementType>::rightFactor(const rsMultiArrayOld<ElementType> &outerProduct,
                                            const rsMultiArrayOld<ElementType> &leftFactor,
                                            rsMultiArrayOld<ElementType> &result)
{
  rsUint32 numResultIndices = outerProduct.numIndices-leftFactor.numIndices;
  rsUint32 *resultIndices   = new rsUint32[numResultIndices];

  rsAssert(rsArrayTools::equal(outerProduct.indexRanges,
    leftFactor.indexRanges, leftFactor.numIndices)); // multiarrays incompatible for this operation

  memcpy(resultIndices, &outerProduct.indexRanges[leftFactor.numIndices],
         numResultIndices*sizeof(rsUint32));
  result.setSize(numResultIndices, resultIndices);

  rsUint32      offset  = (rsUint32) rsArrayTools::maxAbsIndex(leftFactor.data,
                                                   (int) leftFactor.getNumElements());
  ElementType scaler  = ElementType(1) / leftFactor.data[offset];
  offset             *= result.getNumElements();

  for(rsUint32 i = 0; i < result.getNumElements(); i++)
    result.data[i] = scaler * outerProduct.data[offset+i];

  delete[] resultIndices;
}

template<class ElementType>
rsMultiArrayOld<ElementType> rsMultiArrayOld<ElementType>::contract(
  const rsMultiArrayOld<ElementType> &subject, rsUint32 index1, rsUint32 index2)
{
  rsUint32 i, j, k;

  // check some conditions that must hold:
  rsAssert(index1 >= 0 && index1 < subject.numIndices);
  rsAssert(index2 >= 0 && index2 < subject.numIndices);
  rsAssert(index1 != index2);
  rsAssert(subject.indexRanges[index1] == subject.indexRanges[index2]);
  rsUint32 sumRange = subject.indexRanges[index1];
  if( index2 < index1 )
    rsSwap(index1, index2);
    //rsSwap(&index1, &index2);

  // allocate and init object for the result:
  rsMultiArrayOld<ElementType> result;
  rsUint32 numResultIndices = subject.numIndices-2;
  rsUint32 *resultIndices = new rsUint32[numResultIndices];
  j = 0;
  for(i = 0; i < subject.numIndices; i++)
  {
    if( i != index1 && i != index2 )
    {
      resultIndices[j] = subject.indexRanges[i];
      j++;
    }
  }
  result.setSize(numResultIndices, resultIndices);

  // the actual summation loop:
  rsUint32 *subjectIndices = new rsUint32[subject.numIndices];
  for(i = 0; i < result.getNumElements(); i++)
  {
    result.indicesFromOffset(i, resultIndices);

    // copy the result-indices array into the subject-indices array, leaving space for the two
    // summation-indices:
    k = 0;
    for(j = 0; j < subject.numIndices; j++)
    {
      if( j != index1 && j != index2 )
      {
        subjectIndices[j] = resultIndices[k];
        k++;
      }
    }

    // in this inner loop, the summation indices (for which we left the spaces during copying in
    // the previous block) traverse all possible values and the result is accumulated into the
    // appropriate element in the result:
    result.data[i] = ElementType(0);
    for(j = 0; j < sumRange; j++)
    {
      subjectIndices[index1] = subjectIndices[index2] = j;
      result.data[i] += subject.getElement(subjectIndices);
    }
  }

  delete[] resultIndices;
  delete[] subjectIndices;
  return result;
}
