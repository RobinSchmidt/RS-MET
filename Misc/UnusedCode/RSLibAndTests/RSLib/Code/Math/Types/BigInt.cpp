using namespace RSLib;

// construction/destruction:

rsBigInt::rsBigInt(rsInt64 initialValue, rsUint64 baseToUse)
{
  base      = baseToUse;
  numDigits = 0;
  digits    = nullptr;
  setValue(initialValue);
}

rsBigInt::rsBigInt(rsUint32 *digits, rsUint32 numDigits, bool shouldBeNegative, rsUint64 baseToUse)
{
  base = baseToUse;
  int firstNonZeroIndex = rsFirstIndexWithNonZeroValue(digits, numDigits);
  if( firstNonZeroIndex == -1 )
  {
    this->numDigits   = 1;
    firstNonZeroIndex = numDigits-1;
  }
  else
  {
    this->numDigits = numDigits-firstNonZeroIndex;
  }
  allocateMemory();
  rsCopyBuffer(&digits[firstNonZeroIndex], this->digits, this->numDigits);
  negative = shouldBeNegative;
}

rsBigInt::rsBigInt(const rsBigInt& other)
{
  base = other.base;
  this->numDigits = other.numDigits;
  allocateMemory();
  rsCopyBuffer(other.digits, this->digits, numDigits);
  negative = other.negative;
}

// operators:

rsBigInt& rsBigInt::operator=(const rsBigInt& other)
{
  base = other.base;
  if( other.numDigits != numDigits )
  {
    freeMemory();
    numDigits = other.numDigits;
    allocateMemory();
  }
  rsCopyBuffer(other.digits, digits, numDigits);
  negative = other.negative;
  return *this;
}

rsBigInt rsBigInt::operator-() const
{
  rsBigInt result = *this;
  result.negative = !this->negative;
  return result;
}

bool rsBigInt::operator==(const rsBigInt& other) const
{
  rsAssert(this->base == other.base);
  if( other.numDigits != numDigits )
    return false;
  else
    return rsAreBuffersEqual(digits, other.digits, numDigits) && negative == other.negative;
}

bool rsBigInt::operator!=(const rsBigInt& other) const
{
  rsAssert(this->base == other.base);
  return !(*this == other);
}

bool rsBigInt::operator>=(const rsBigInt& other) const
{
  rsAssert(this->base == other.base);
  return (*this > other) || (*this == other);
}

bool rsBigInt::operator<=(const rsBigInt& other) const
{
  rsAssert(this->base == other.base);
  return !(*this > other);
}

bool rsBigInt::operator>(const rsBigInt& other) const
{
  rsAssert(this->base == other.base);
  if( this->isNegative() && !other.isNegative() )
    return false;
  else if( !this->isNegative() && other.isNegative())
    return true;
  else if( !this->isNegative() && !other.isNegative() )
  {
    if( this->numDigits > other.numDigits )
      return true;
    else if( this->numDigits < other.numDigits )
      return false;
    else
      return rsCompare(this->digits, other.digits, numDigits) == 1;
  }
  else
  {
    if( this->numDigits > other.numDigits )
      return false;
    else if( this->numDigits < other.numDigits )
      return true;
    else
      return rsCompare(this->digits, other.digits, numDigits) == -1;
  }
}

bool rsBigInt::operator<(const rsBigInt& other) const
{
  rsAssert(this->base == other.base);
  return !(*this >= other);
}

rsBigInt rsBigInt::operator+(const rsBigInt &other)
{
  rsAssert(this->base == other.base);
  int  operandLength   = rsMax(numDigits, other.numDigits);
  int  maxResultLength = operandLength + 1;
  int  tempLength      = 2*operandLength + maxResultLength;
  rsUint32 *tempDigits = new rsUint32[tempLength];
  rsFillWithZeros(tempDigits, tempLength);
  rsCopyBuffer(this->digits, &tempDigits[1*operandLength - this->numDigits], this->numDigits);
  rsCopyBuffer(other.digits, &tempDigits[2*operandLength - other.numDigits], other.numDigits);

  bool resultIsNegative;
  if( this->hasSameSignAs(other) )
  {
    //resultIsNegative = flags.isFlagTrue(negative);
    resultIsNegative = negative;
    tempDigits[2*operandLength] = add(tempDigits, &tempDigits[operandLength], operandLength,
                                      &tempDigits[2*operandLength+1], base);
    rsBigInt result(&tempDigits[2*operandLength], maxResultLength, resultIsNegative, base);
    delete[] tempDigits;
    return result;
  }
  else
  {
    // operands have different signs, subtract the smaller absolute value from the larger:
    int compareResult = rsCompare(tempDigits, &tempDigits[operandLength], operandLength);
    if( compareResult == 0 )     // both operands have equal absolute value
    {
      delete[] tempDigits;
      return rsBigInt(0, base);
    }
    else if( compareResult > 0 ) // left operand has greater absolute value
    {
      // |L| > |R| -> result negative, iff L < 0:
      subtract(tempDigits, &tempDigits[operandLength], operandLength,
               &tempDigits[2*operandLength+1], base);
      resultIsNegative = this->isNegative();
    }
    else if( compareResult < 0 ) // left operand has smaller absolute value
    {
      // |L| < |R| -> result negative, iff R < 0:
      subtract(&tempDigits[operandLength], tempDigits, operandLength,
               &tempDigits[2*operandLength+1], base);
      resultIsNegative = other.isNegative();
    }
    rsBigInt result(&tempDigits[2*operandLength+1], maxResultLength-1, resultIsNegative, base);
    delete[] tempDigits;
    return result;
  }
}

rsBigInt rsBigInt::operator-(const rsBigInt &other)
{
  rsAssert(this->base == other.base);
  rsBigInt tmp = other;
  tmp.negative = !tmp.negative;
  return *this + tmp;
}

rsBigInt rsBigInt::operator*(const rsBigInt &other)
{
  rsAssert(this->base == other.base);
  int resultLength = numDigits + other.numDigits;
  rsUint32 *tempResult = new rsUint32[resultLength];
  multiply(digits, numDigits, other.digits, other.numDigits, tempResult, base);
  bool resultIsNegative;
  if( this->hasSameSignAs(other) )
    resultIsNegative = false;
  else
    resultIsNegative = true;
  rsBigInt result(tempResult, resultLength, resultIsNegative, base);
  delete[] tempResult;
  return result;
}

rsBigInt rsBigInt::operator/(const rsBigInt &other)
{
  rsAssert(this->base == other.base);
  rsBigInt q, r;
  divMod(other, q, r);
  return q;
}

rsBigInt rsBigInt::operator%(const rsBigInt &other)
{
  rsAssert(this->base == other.base);
  rsBigInt q, r;
  divMod(other, q, r);
  return r;
}

rsBigInt& rsBigInt::operator+=(const rsBigInt &other)
{
  *this = *this + other;
  return *this;
}
rsBigInt& rsBigInt::operator-=(const rsBigInt &other)
{
  *this = *this - other;
  return *this;
}
rsBigInt& rsBigInt::operator*=(const rsBigInt &other)
{
  *this = *this * other;
  return *this;
}
rsBigInt& rsBigInt::operator/=(const rsBigInt &other)
{
  *this = *this / other;
  return *this;
}
rsBigInt& rsBigInt::operator%=(const rsBigInt &other)
{
  *this = *this % other;
  return *this;
}

rsBigInt& rsBigInt::operator++()
{
  *this = *this + 1;
  return *this;
}

rsBigInt& rsBigInt::operator--()
{
  *this = *this - 1;
  return *this;
}

// setup:

void rsBigInt::setValue(rsInt64 newValue)
{
  if( newValue == 0 )
  {
    // maybe wrap into a setZero function:
    negative = false;
    if( numDigits != 1 )
    {
      freeMemory();
      numDigits = 1;
      allocateMemory();
    }
    digits[0] = 0;
  }

  negative = newValue < 0;
  newValue = rsAbs(newValue);

  rsUint64 r = newValue;
  rsUint64 q;

  // new:
  rsUint32 a[64]; // will contain the digits in reverse order
  rsUint32 newNumDigits = 0;
  while( r > 0 )
  {
    q = r / base;
    r = r % base;
    a[newNumDigits] = (rsUint32)r;
    newNumDigits++;
    r = q;
  }

  if( newNumDigits != numDigits )
  {
    freeMemory();
    numDigits = newNumDigits;
    allocateMemory();
  }
  for(rsUint32 i = 0; i < numDigits; i++)
    digits[i] = a[numDigits-1-i];

  /*
  // old:
  // using a dynamically growing array here sucks. it would be better to allocate enough memory
  // at once - for that, we would need an upper bound for the number of digits needed.
  // 64 digits should be such an upper bound:
  rsArray<rsUint32> a;  // will contain the digits in reverse order
  while( r > 0 )
  {
    q = r / base;
    r = r % base;
    a.appendElement((rsUint32)r);
    r = q;
  }

  if( a.getNumElements() != numDigits )
  {
    freeMemory();
    numDigits = a.getNumElements();
    allocateMemory();
  }
  for(rsUint32 i = 0; i < numDigits; i++)
    digits[i] = a.getElement(numDigits-1-i);
  */
}

// inquiry:

void rsBigInt::divMod(const rsBigInt& divisor, rsBigInt& quotient, rsBigInt& remainder) const
{
  rsAssert(this->base == divisor.base);
  rsAssert(divisor != rsBigInt(0), "rsBigInt division by zero.");

  int resultLength = rsMax(this->numDigits, divisor.numDigits);
  int tempLength   = 4*resultLength;
  rsUint32 *tempDigits = new rsUint32[tempLength];
  rsFillWithZeros(tempDigits, tempLength);
  rsCopyBuffer(this->digits,   &tempDigits[1*resultLength - this->numDigits],   this->numDigits);
  rsCopyBuffer(divisor.digits, &tempDigits[2*resultLength - divisor.numDigits], divisor.numDigits);

  divide(tempDigits, &tempDigits[resultLength], resultLength, &tempDigits[2*resultLength],
         &tempDigits[3*resultLength], base);
  bool quotientIsNegative  = !this->hasSameSignAs(divisor);
  bool remainderIsNegative =  this->isNegative();

  quotient  = rsBigInt(&tempDigits[2*resultLength], resultLength, quotientIsNegative,  base);
  remainder = rsBigInt(&tempDigits[3*resultLength], resultLength, remainderIsNegative, base);
  delete[] tempDigits;
}

rsString rsBigInt::toString(rsUint64 stringBase) const
{
  rsArray<rsUint32> a = changeBaseForInteger(digits, numDigits, base, stringBase);
  rsString result;
  result.ensureAllocatedSize(a.getNumElements() + negative);
  if( negative )
    result.appendElement('-');
  for(int i = 0; i < a.getNumElements(); i++)
    result.appendElement(rsUint32ToChar(a[i]));
  return result;
}

// misc:

rsBigInt rsBigInt::fromString(const rsString& s, rsUint64 stringBase, rsUint64 baseToUse)
{
  rsUint32 numDigits = s.getNumElements();
  bool negative     = false;
  if( s.getElement(0) == '-' )
  {
    negative   = true;
    numDigits -= 1;
  }
  rsUint32 *tmp = new rsUint32[numDigits];
  for(rsUint32 i = 0; i < numDigits; i++)
    tmp[i] = rsCharToUint32(s.getElement(i+(int)negative));
  rsArray<rsUint32> a = changeBaseForInteger(tmp, numDigits, stringBase, baseToUse);
  delete[] tmp;
  return rsBigInt(a.getRawData(), a.getNumElements(), negative, baseToUse);
}

rsString RSLib::rsConvertIntegerIntoBase(rsString value, rsUint64 oldBase, rsUint64 newBase)
{
  return rsBigInt::fromString(value, oldBase).toString(newBase);
    // it seems silly to convert back and forth through a string. the code must be old, from a
    // point where no other option were available, replace this
}

rsBigInt RSLib::rsUnityValue(rsBigInt value)
{
  return rsBigInt(1, value.getBase());
}
