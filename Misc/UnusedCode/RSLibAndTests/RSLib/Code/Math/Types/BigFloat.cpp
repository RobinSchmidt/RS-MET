using namespace RSLib;

rsBigFloat::rsBigFloat(rsUint32 precision, rsUint64 baseToUse)
{
  numDigits = precision;
  base      = baseToUse;
  allocateMemory();
  setZero();
}

rsBigFloat::rsBigFloat(const char s[] /*, rsUint64 baseToUse*/)
{
  numDigits = 0;
  digits  = nullptr;
  //fromString(s, baseToUse);
  fromString(s /*, baseToUse*/);
}

rsBigFloat::rsBigFloat(double value, rsUint32 precisionToUse, rsUint64 baseToUse)
{
  numDigits = precisionToUse;
  base      = baseToUse;
  allocateMemory();
  fromDouble(value, precisionToUse, baseToUse);
}

rsBigFloat::rsBigFloat(const rsBigInt& bigInt)
{
  numDigits = bigInt.getNumDigits();
  negative  = bigInt.isNegative();
  base      = bigInt.getBase();
  exponent  = bigInt.getNumDigits();
  allocateMemory();
  bigInt.copyDigitsToBuffer(digits);
}

rsBigFloat::rsBigFloat(const rsBigFloat& other)
{
  negative  = other.negative;
  base      = other.base;
  exponent  = other.exponent;
  numDigits = other.numDigits;
  allocateMemory();
  rsCopyBuffer(other.digits, this->digits, numDigits);
}

// operators:

rsBigFloat& rsBigFloat::operator=(const rsBigFloat& other)
{
  if( other.numDigits != numDigits )
  {
    freeMemory();
    numDigits = other.numDigits;
    allocateMemory();
  }
  rsCopyBuffer(other.digits, digits, numDigits);
  exponent = other.exponent;
  negative = other.negative;
  base     = other.base;
  return *this;
}

rsBigFloat rsBigFloat::operator-()
{
  rsBigFloat result = *this;
  result.negative = !this->negative;
  return result;
}

// maybe we should have a macro that generates comparison-operators like that for any class:
bool rsBigFloat::operator==(const rsBigFloat& other) const
{
  return compare(*this, other) == 0;
}
bool rsBigFloat::operator!=(const rsBigFloat& other) const
{
  return compare(*this, other) != 0;
}
bool rsBigFloat::operator>=(const rsBigFloat& other) const
{
  return compare(*this, other) >= 0;
}
bool rsBigFloat::operator<=(const rsBigFloat& other) const
{
  return compare(*this, other) <= 0;
}
bool rsBigFloat::operator>(const rsBigFloat& other) const
{
  return compare(*this, other) > 0;
}
bool rsBigFloat::operator<(const rsBigFloat& other) const
{
  return compare(*this, other) < 0;
}

void rsBigFloat::subtractZeroExtended(const rsBigFloat& a, const rsBigFloat& b, rsBigFloat& d)
{
  rsAssert(a.numDigits == b.numDigits); // the algorithm assumes this

  int N         = a.numDigits;
  int shift     = rsMin(N, abs(rsClippedDifference(a.exponent, b.exponent)));
  int opLength  = N + shift;
  rsUint32 *opL = new rsUint32[opLength];
  rsUint32 *opR = new rsUint32[opLength];
  rsUint32 *sum = new rsUint32[opLength+2]; // 2 more for carry and shifting
  rsCopyBuffer(a.digits, opL, a.numDigits);
  rsCopyBuffer(b.digits, opR, b.numDigits);
  rsFillWithZeros(&opL[a.numDigits], opLength-a.numDigits);
  rsFillWithZeros(&opR[b.numDigits], opLength-b.numDigits);

  rsRightShift(opR, opLength, shift);

  sum[0]  = rsBigNumber::subtract(opL, opR, opLength, &sum[1], a.base);
  d.negative = false;

  // this should be common with the other version except for the cleanup:
  sum[opLength+1] = 0;
  d.exponent = rsMax(a.exponent, b.exponent);
  if( sum[0] > 0 )
  {
    rsRightShift(sum, opLength+2, 1);
    d.incrementExponent(1);
  }
  sum[0] = roundHeadAccordingToTail(&sum[1], N, a.base, &sum[N+1], shift, a.base);
  if( sum[0] > 0 )
  {
    rsRightShift(sum, opLength+1, 1); // last is garbage
    d.incrementExponent(1);
  }
  rsCopyBuffer(&sum[1], d.digits, N);
  d.normalize();

  delete[] opL;
  delete[] opR;
  delete[] sum;
}

void rsBigFloat::addOrSubtractWithTailCopy(const rsBigFloat& a, const rsBigFloat& b, rsBigFloat& s)
{
  rsAssert(a.numDigits == b.numDigits); // the algorithm assumes this

  int N         = a.numDigits; // we can assume that they are the same
  int shift     = rsMin(N, abs(rsClippedDifference(a.exponent, b.exponent))); // rename to shiftAmount
  int opLength  = N + shift;

  rsUint32 *sum       = new rsUint32[opLength+2];
  rsUint32 *shiftedOp = nullptr;
  rsUint32 *tmpOp     = nullptr;
  rsUint32 *opL, *opR;
  if( a.exponent > b.exponent )
  {
    tmpOp     = new rsUint32[N+shift];
    opL       = a.digits;
    opR       = tmpOp;
    shiftedOp = opR;
    rsFillWithZeros(tmpOp, shift);
    rsCopyBuffer(b.digits, &tmpOp[shift], a.numDigits);
  }
  else if( b.exponent > a.exponent )
  {
    tmpOp     = new rsUint32[N+shift];
    opL       = tmpOp;
    opR       = b.digits;
    shiftedOp = opL;
    rsFillWithZeros(tmpOp, shift);
    rsCopyBuffer(a.digits, &tmpOp[shift], a.numDigits);
  }
  else
  {
    opL = a.digits;
    opR = b.digits;
    shiftedOp = nullptr;
  }

  if( a.hasSameSignAs(b) )
  {
    sum[0] = rsBigNumber::add(opL, opR, N, &sum[1], a.base);
    s.negative = a.negative;
  }
  else
  {
    bool rightGreater = rsCompare(opL, opR, N) == -1;
    if( rightGreater )
      sum[0] = rsBigNumber::subtract(opR, opL, N, &sum[1], a.base);
    else
      sum[0] = rsBigNumber::subtract(opL, opR, N, &sum[1], a.base);

    if( !rightGreater && a.negative ||  rightGreater && b.negative  )
      s.negative = true;
  }
  rsCopyBuffer(&shiftedOp[N], &sum[N+1], shift);
  sum[opLength+1] = 0;

  s.exponent = rsMax(a.exponent, b.exponent);
  if( sum[0] > 0 )
  {
    rsRightShift(sum, opLength+2, 1);
    s.incrementExponent(1);
  }
  sum[0] = roundHeadAccordingToTail(&sum[1], N, a.base, &sum[N+1], shift, a.base);
  if( sum[0] > 0 )
  {
    rsRightShift(sum, opLength+1, 1); // last is garbage
    s.incrementExponent(1);
  }
  rsCopyBuffer(&sum[1], s.digits, N);
  s.normalize();

  delete[] sum;
  delete[] tmpOp;
}

void rsBigFloat::addOrSubtractAllCases(const rsBigFloat& a, const rsBigFloat& b, rsBigFloat& r)
{
  int N         = rsMax(a.numDigits, b.numDigits);
  int shift     = rsMin(N, abs(rsClippedDifference(a.exponent, b.exponent)));
  int opLength  = N + shift;
  rsUint32 *opL = new rsUint32[opLength];
  rsUint32 *opR = new rsUint32[opLength];
  rsUint32 *sum = new rsUint32[opLength+2]; // 2 more for carry and shifting
  rsCopyBuffer(a.digits, opL, a.numDigits);
  rsCopyBuffer(b.digits, opR, b.numDigits);
  rsFillWithZeros(&opL[a.numDigits], opLength-a.numDigits);
  rsFillWithZeros(&opR[b.numDigits], opLength-b.numDigits);

  if( a.exponent > b.exponent )
    rsRightShift(opR, opLength, shift);
  else if( b.exponent > a.exponent )
    rsRightShift(opL, opLength, shift);

  if( a.hasSameSignAs(b) )
  {
    sum[0] = rsBigNumber::add(opL, opR, opLength, &sum[1], a.base);
    r.negative = a.negative;
  }
  else
  {
    bool rightGreater = rsCompare(opL, opR, N) == -1;
    if( rightGreater )
      sum[0] = rsBigNumber::subtract(opR, opL, opLength, &sum[1], a.base);
    else
      sum[0] = rsBigNumber::subtract(opL, opR, opLength, &sum[1], a.base);

    if( !rightGreater && a.negative ||  rightGreater && b.negative  )
      r.negative = true;
  }
  sum[opLength+1] = 0;

  r.exponent = rsMax(a.exponent, b.exponent);
  if( sum[0] > 0 )
  {
    rsRightShift(sum, opLength+2, 1);
    r.incrementExponent(1);
  }
  sum[0] = roundHeadAccordingToTail(&sum[1], N, a.base, &sum[N+1], shift, a.base);


  if( sum[0] > 0 )
  {
    rsRightShift(sum, opLength+1, 1); // last is garbage
    r.incrementExponent(1);
  }
  rsCopyBuffer(&sum[1], r.digits, N);
    // we may get rid of the rightshift by either copying from sum[0] (and incrementing) or
    // copying from sum[1] without incrementing - then we also need one digit less for the sum


  r.normalize();

  delete[] opL;
  delete[] opR;
  delete[] sum;
}

void rsBigFloat::add(const rsBigFloat& a, const rsBigFloat& b, rsBigFloat& s)
{
  rsAssert(a.base == b.base);

  // the function below handles all cases uniformly and may even deal with operands a, b which
  // have different numbers of digits, but that comes at the cost of more memory allocation and
  // doing always the full zero-extended addition/subtraction.
  addOrSubtractAllCases(a, b, s);

  // below is the switch that distiguishes the 2 cases: (1): we need to borrow-from-zero due
  // right-shifting an operator that is subtracted - in this case, we must do the full
  // zero-extended version. (2) no borrow-from zero occurs and the tail of the sum is equal to
  // one of the operand's tails - in this case, we can copy the tail.  i don't know, if it's better to treat the cases seperately
  // treating them uniformly certainly needs less code. maybe do performance tests
  /*
  bool borrowFromZeroInTail = a.exponent > b.exponent && !a.negative && b.negative;
  if( borrowFromZeroInTail )
  {
    subtractZeroExtended(a, b, s);
    return;
  }
  else
  {
    addOrSubtractWithTailCopy(a, b, s);
    return;
  }
  */
}

void rsBigFloat::divide(const rsBigFloat& a, const rsBigFloat& b, rsBigFloat& result)
{
  rsAssert(a.base == b.base);
  rsAssert(!b.isZero(),"rsBigFloat: division by zero");

  rsUint32 Nq = rsMax(a.numDigits, result.numDigits) + b.numDigits;
  rsUint32 *q = new rsUint32[Nq];

  bool overflow = divideFractions(a.digits, a.numDigits, b.digits, b.numDigits, q, Nq, a.base);
  if( overflow )
  {
    rsFillWithZeros(q, Nq);
    q[0] = 1;
  }
  // i think, we don't need to handle overflow - it should be impossible

  result.setValue(q, Nq, a.exponent - b.exponent + b.numDigits + overflow,
                  rsExOr(a.negative, b.negative));
  delete[] q;
}

void rsBigFloat::multiply(const rsBigFloat& a, const rsBigFloat& b, rsBigFloat& result)
{
  rsAssert(a.base == b.base);

  rsUint32 Np = a.numDigits + b.numDigits;
  rsUint32 *p = new rsUint32[Np];
  rsBigNumber::multiply(a.digits, a.numDigits, b.digits, b.numDigits, p, a.base);
  result.setValue(p, Np, a.exponent + b.exponent, rsExOr(a.negative, b.negative));
  delete[] p;
}

rsBigFloat rsBigFloat::operator+(const rsBigFloat &other)
{
  rsBigFloat result(rsMax(numDigits, other.numDigits), base);
  if( this->numDigits > other.numDigits )
    add(*this, other.withPrecision(this->numDigits), result);
  else if( this->numDigits < other.numDigits )
    add(this->withPrecision(other.numDigits), other, result);
  else
    add(*this, other, result);
  return result;
}

rsBigFloat rsBigFloat::operator-(const rsBigFloat &other)
{
  rsBigFloat tmp = other;
  tmp.negative = !tmp.negative;
  return *this + tmp;
}

rsBigFloat rsBigFloat::operator*(const rsBigFloat &other)
{
  rsBigFloat result(rsMax(numDigits, other.numDigits), base);
  multiply(*this, other, result);
  return result;
}

rsBigFloat rsBigFloat::operator/(const rsBigFloat &other)
{
  rsBigFloat result(rsMax(numDigits, other.numDigits), base);
  divide(*this, other, result);
  return result;
}

// these operators look exactly the same as in rsBigInt and in fact, it's likely that they look
// like that in just about any class - maybe write a macro that generates them to avoid the
// code-duplication, add ++ and -- operator
rsBigFloat& rsBigFloat::operator+=(const rsBigFloat &other)
{
  *this = *this + other;
  return *this;
}
rsBigFloat& rsBigFloat::operator-=(const rsBigFloat &other)
{
  *this = *this - other;
  return *this;
}
rsBigFloat& rsBigFloat::operator*=(const rsBigFloat &other)
{
  *this = *this * other;
  return *this;
}
rsBigFloat& rsBigFloat::operator/=(const rsBigFloat &other)
{
  *this = *this / other;
  return *this;
}

void rsBigFloat::setZero()
{
  rsFillWithZeros(digits, numDigits);
  exponent = 0;
  negative = false;
}

void rsBigFloat::setData(rsUint32 *digits, rsUint32 numDigits, rsUint64 base, rsInt32 exponent,
                         bool negative, bool shouldNormalize)
{
  if( this->numDigits != numDigits )
  {
    freeMemory();
    this->numDigits = numDigits;
    allocateMemory();
  }
  this->base      = base;
  this->exponent  = exponent;
  this->negative  = negative;
  rsCopyBuffer(digits, this->digits, numDigits);
  if( shouldNormalize )
    normalize();
}

void rsBigFloat::setValue(rsUint32 *newDigits, rsUint32 numNewDigits, rsInt32 exponent,
                          bool negative)
{
  int start = rsFirstIndexWithNonZeroValue(newDigits, numNewDigits);
  if( start == -1 )
  {
    setZero();
    return;
  }
  newDigits    += start;
  numNewDigits -= start;

  this->negative = negative;
  setExponent(exponent-start); // must be done before handleRoundingOverflow()

  if( numNewDigits > numDigits )
  {
    rsCopyBuffer(newDigits, digits, numDigits);
    bool overflow = roundHeadAccordingToTail( digits,               numDigits,              base,
                                             &newDigits[numDigits], numNewDigits-numDigits, base);
    if( overflow )
      handleRoundingOverflow();
  }
  else
  {
    rsCopyBuffer(newDigits, digits, numNewDigits);
    rsFillWithZeros(&digits[numNewDigits], numDigits-numNewDigits);
  }
}

// inquiry:

void rsBigFloat::fromString(const char s[] /*, rsUint64 baseToUse*/)
{
  rsUint32 length = (rsUint32)strlen(s);
  rsUint32 *tmp   = new rsUint32[length];
  rsFillWithZeros(tmp, length);

  base                  = 10;  /*baseToUse;*/
  exponent              = 0;
  negative              = s[0] == '-';
  bool nonZeroDigitSeen = false;
  bool exponentSeen     = false;
  bool dotSeen          = false;
  bool hasSign          = s[0] == '+' || s[0] == '-';
  rsUint32 newNumDigits = 0;
  rsUint32 i;
  for(i = hasSign; i < length; i++)
  {
    if( exponentSeen == false )
    {
      // copy digits:
      if( s[i] >= '0' && s[i] <= '9' ) // more generally: doesCharacterRepresentDigit(s[i])
      {                                // for 10 < base <= 36
        tmp[newNumDigits] = rsBigInt::rsCharToUint32(s[i]);
        nonZeroDigitSeen |= s[i] != '0';
        newNumDigits++;
      }
      else if( s[i] == '.' )
        dotSeen = true;
      else if( s[i] == 'e' || s[i] == 'E' )
        exponentSeen = true;

      // adjust exponent:
      if( !exponentSeen )
      {
        if( !dotSeen && nonZeroDigitSeen )
          exponent++;
        if( dotSeen && !nonZeroDigitSeen && s[i] != '.'  )
          exponent--;
      }
    }
    else
    {
      // add/subtract the exponent given in the string:
      if( s[i] == '-' )
        exponent -= strtol(&s[i+1], NULL, (int)base);
      else if( s[i] == '+' )
        exponent += strtol(&s[i+1], NULL, (int)base);
      else
        exponent += strtol(&s[i],   NULL, (int)base);
      break;
    }
  }

  // catch all-zeros case:
  if( rsIsAllZeros(tmp, newNumDigits) )
  {
    if( numDigits != 1 )
    {
      freeMemory();
      numDigits = 1;
      allocateMemory();
    }
    setZero();
    delete[] tmp;
    return;
  }

  // count leading and trailing zeros:
  int numLeadingZeros = 0;
  for(i = 0; i < length; i++)
  {
    if( tmp[i] == 0 )
      numLeadingZeros++;
    else
      break;
  }
  int numTrailingZeros = 0;
  for(i = length-1; i > 0; i--)
  {
    if( tmp[i] == 0 )
      numTrailingZeros++;
    else
      break;
  }

  // copy digits from tmp-buffer into member-variable (scrapping leading and trailing zeros):
  newNumDigits = length - numLeadingZeros - numTrailingZeros;
  if( newNumDigits != numDigits )
  {
    freeMemory();
    numDigits = newNumDigits;
    allocateMemory();
  }
  rsCopyBuffer(&tmp[numLeadingZeros], digits, numDigits);
  delete[] tmp;
}

int rsBigFloat::compare(const rsBigFloat& a, const rsBigFloat& b)
{
  rsAssert(a.base == b.base);     // arguments must have same "type"

  if( a.isZero() && b.isZero() )  // different representations of zero are considered equal
    return 0;

  if( a.negative && !b.negative )
    return -1;
  else if( !a.negative && b.negative )
    return +1;
  else
  {
    bool bothNegative = a.negative;
    if( a.exponent > b.exponent )
    {
      if( !bothNegative )
        return +1;
      else
        return -1;
    }
    else if( a.exponent < b.exponent )
    {
      if( !bothNegative )
        return -1;
      else
        return +1;
    }
    else
    {
      if( !bothNegative )
        return  rsCompare(a.digits, a.numDigits, b.digits, b.numDigits);
      else
        return -rsCompare(a.digits, a.numDigits, b.digits, b.numDigits);
    }
  }
}

rsString rsBigFloat::toString(rsUint32 numDigits, rsUint32 stringBase) const
{
  rsBigFloat tmp = this->toBase(stringBase, numDigits);
  rsString result;
  result.ensureAllocatedSize(numDigits + negative);
  if( negative )
    result.appendElement('-');
  result += rsString("0.");
  for(rsUint32 i = 0; i < numDigits; i++)
    result.appendElement(rsBigInt::rsUint32ToChar(tmp.digits[i]));
  result += rsString("e");
  result += rsString(tmp.exponent);
  return result;
}

rsBigFloat rsBigFloat::withPrecision(rsUint32 precisionToUse) const
{
  rsBigFloat result(precisionToUse, base);

  result.negative = negative;
  result.exponent = exponent;
  rsUint32 copyLength = rsMin(numDigits, precisionToUse);
  rsCopyBuffer(digits, result.digits, copyLength);
  if( numDigits > precisionToUse )
  {
    bool overflow = roundHeadAccordingToTail(result.digits, precisionToUse, base,
      &digits[precisionToUse], numDigits-precisionToUse, base);
    if( overflow )
      result.handleRoundingOverflow();
  }
  else if( numDigits > precisionToUse )
    rsFillWithZeros(&result.digits[precisionToUse], numDigits-precisionToUse);

  return result;
}

rsBigFloat rsBigFloat::toBase(rsUint64 newBase, rsUint32 numDigitsToUse, bool exact) const
{
  // re-interpret fraction as integer (amounts to multiplying it by base^numDigits) and divide this
  // intermediate result by base^numDigits/base^exponent = base^(numDigits-exponent):

  // convert scaled fraction exactly (by re-interpreting it as integer):
  rsArray<rsUint32> a;
  rsUint32 tmpNumDigits = numDigits;
  if( exact == true )
    a = changeBaseForInteger(digits, tmpNumDigits, base, newBase);
  else
  {
    rsAssert(false);
      // fast (and approximate) conversion not yet implemented. the idea is to assign a possibly
      // lower value to tmpNumDigits than numDigits - we need a rule to determine the required
      // number of digits for mantissa conversion, such that in the  base-converted result only the
      // last digit may be wrong

    a = changeBaseForInteger(digits, tmpNumDigits, base, newBase);
  }
  rsBigInt tmp(a.getRawData(), a.getNumElements(), false, newBase);

  // un-scale the scaled converted fraction by dividing by an appropriate power of the old base
  // (represented exactly as to-float-converted-integer in the new base):
  rsInt64 ex = (rsInt64)tmpNumDigits - (rsInt64)exponent;
  rsBigFloat scaler(rsPow(rsBigInt(base, newBase), abs(ex))); // maybe use rsAbs(ex)
  rsBigFloat result(numDigitsToUse, newBase);
  if( ex > 0 )
    divide(  rsBigFloat(tmp), scaler, result);
  else
    multiply(rsBigFloat(tmp), scaler, result);

  return result;
}
// edit: i just found this article - may be worth a read in this context:
// https://www.exploringbinary.com/number-of-digits-required-for-round-trip-conversions/


/*
// the apporach of converting the fractionand muliplying the result by an appropriate power of the
// old base (as determined by the (old) exponent) introduces large rounding errors - the rounding
// error cuased by the conversion of the fraction may be magnified by the subsequent multiplication
rsBigFloat rsBigFloat::toBase(rsUint64 newBase, rsUint32 numDigitsToUse) const
{
  static const int numGuardDigits = 32; // to eliminate rounding errors - hmm - maybe this is not good?
  rsBigFloat result(numDigitsToUse + numGuardDigits, newBase);

  // perhaps, it's not good to change the base for the fraction first and then multiply because
  // the rounding error will be magnified by the multiplier

  // convert normalized fraction digits:
  bool overflow = rsBigNumber::changeBaseForFraction(digits, numDigits, base,
                    result.digits, result.numDigits, result.base);
  if( overflow )
    result.handleRoundingOverflow();

  rsUint32 tmp[1000];
  result.copyDigitsToBuffer(tmp);


  // multiply result by appropriate power of our base, represented as rsBigFloat in the new base:

  rsBigInt test = rsBigInt(base, result.base);

  rsBigInt multiplier = rsBigInt::pow(rsBigInt(base, result.base), abs(exponent));
     // when converting 9 from base 10 to 2, the multiplier has a leading 2 (in base 2) which
     // should be impossible

  if( exponent > 0 )
    result *= rsBigFloat(multiplier);
  else if( exponent < 0 )
    result /= rsBigFloat(multiplier);

  result.copyDigitsToBuffer(tmp);

  result.normalize();

  result.copyDigitsToBuffer(tmp);

  return result.withPrecision(numDigitsToUse);
    // .withPrecision is necessary, because of guard-digits and the multiplication with
    // "multiplier" may have increased the precision of the result. the precision of the multiplier
    // is not known beforehand - it equals the number of digits, required to represent is exactly,
    // because the rsBigInt class uses always as many digits as necessarry to represent a number
}
*/


rsBigFloat rsBigFloat::mimicBuiltInType(double value) const
{
  rsUint32 bits[64];
  rsGetBits64(value, bits);
  rsBigFloat result((rsUint32)53, (rsUint64)2);
  result.negative = value < 0.0;
  result.exponent = rsExtractExponentFromDouble(value) + 1; // +1 because of implied leading 1?
  if( !rsIsDenormal(value) )  // test this...
    result.digits[0] = 1;                                 // set implied leading 1
  for(int i = 1; i <= 52; i++)
    result.digits[i] = bits[i+11];

  rsUint32 tmp[1000];
  result.copyDigitsToBuffer(tmp);


  return result;
}

void rsBigFloat::fromDouble(double value, rsUint32 precisionToUse, rsUint64 baseToUse)
{
  *this = mimicBuiltInType(value).toBase(baseToUse, precisionToUse);
}

double rsBigFloat::toDouble() const
{
  rsBigFloat tmp = this->toBase(2, 53);
  rsUint32 bits[64];

  rsFillWithZeros(bits, 64);
  if( !tmp.isNormalized() )
  {
    //rsAssert(false); // denormal handling not yet implemented
  }
  for(int i = 1; i <= 52; i++)
    bits[i+11] = tmp.digits[i];

  double result;
  rsSetBits64(&result, bits);
  rsSetExponentOfDouble(&result, tmp.exponent-1);
  if( negative )
    result = -result;

  return result;
}


rsBigFloat rsBigFloat::computePi(rsUint32 precisionToUse, rsUint64 baseToUse)
{
  rsAssert(false); // not yet implemented
  return rsBigFloat(precisionToUse, baseToUse); // preliminary
}

// misc:

void rsBigFloat::normalize()
{
  int shiftAmount = rsFirstIndexWithNonZeroValue(digits, numDigits);

  if( (rsInt64)exponent - shiftAmount < minExponent )
    shiftAmount = exponent - minExponent; // unnormalized numbers

  if( shiftAmount > 0 )
  {
    rsLeftShift(digits, numDigits, shiftAmount);
    exponent -= shiftAmount;
  }
}

void rsBigFloat::incrementExponent(rsUint32 amount)
{
  exponent = rsClippedSum(exponent, amount, minExponent, maxExponent);
}

void rsBigFloat::decrementExponent(rsUint32 amount)
{
  exponent = rsClippedDifference(exponent, amount, minExponent, maxExponent);
}

void rsBigFloat::handleRoundingOverflow()
{
  incrementExponent(1);
  digits[0] = 1;
  rsFillWithZeros(&digits[1], numDigits-1);
}

//-------------------------------------------------------------------------------------------------
// mathematical functions for rsBigFloat:

rsBigFloat RSLib::rsSqrt(const rsBigFloat& x)
{
  rsUint32 yt[100], dt[100]; // for debug
  x.copyDigitsToBuffer(yt);



  // initialize:
  double xd = x.toDouble();
  rsBigFloat y(x.getNumDigits(), x.getBase());
  y.fromDouble(1.0/rsSqrt(xd), x.getNumDigits(), x.getBase());




  // iterate:
  rsBigFloat one( 1.0, x.getNumDigits(), x.getBase());
  rsBigFloat half(0.5, x.getNumDigits(), x.getBase());
  rsBigFloat delta;
  while( true )
  {
    delta = y*(half*(one-y*y*x));

    // for debug:
    y.copyDigitsToBuffer(yt);
    delta.copyDigitsToBuffer(dt);

    if( delta.isSmallComparedTo(y) )
      break;

    //if( delta.isZero() )
    //  break; // this may lead to limit-cycles - we should instead use
             // delta.isSmallComparedTo(y) which returns true, if
             // y.exponent - delta.exponent >= y.numDigits-1 or delta == zero
             // it's small, if it's with one 2 ulps of y

    y += delta;
  }


  return y*x;
}


double rsAGM(double a, double b) // just for testing the algorithm
{
  double tmp;
  while( a != b )
  {
    tmp = 0.5*(a+b);
    b   = rsSqrt(a*b);
    a   = tmp;
  }
  return a;
}
rsBigFloat RSLib::rsAGM(const rsBigFloat& a, const rsBigFloat& b)
{
  double agm = ::rsAGM(4.0, 2.0);

  rsUint32 N = rsMax(a.getNumDigits(), b.getNumDigits());
  rsUint64 B = a.getBase();


  rsUint32 ab[100], bb[100]; // for debug

  rsBigFloat at(a), bt(b), tmp, half;
  half.fromDouble(0.5, N, B);
  while( at != bt )
  {
    at.copyDigitsToBuffer(ab);
    bt.copyDigitsToBuffer(bb);

    tmp = half *(at+bt);
    bt  = rsSqrt(at*bt);
    at  = tmp;

    at.copyDigitsToBuffer(ab);
    bt.copyDigitsToBuffer(bb);
    int dummy = 0;

  }
  return at;
}

// Notes / stuff to do:
// epsilon = 1 / base^numDigits (difference between 1 and closest lower number)
// unnormalized: digits[numDigits-1] = 1, all other 0, exponent = 0
//   normalized: digits[0] = 1, all others 0, exponent = -(numDigits-1) -> check this
// hmm: http://www.cs.mcgill.ca/~chang/teaching/cs350/slides/ieee.pdf, says epsilon is the
// difference bewteen 1 and the next *larger* number, so it should be 1 / base^(numDigits-1)?
// here are yet 2 other another definitions: http://en.wikipedia.org/wiki/Machine_epsilon
// "the maximum error that can occur when rounding to the unit value".
// "Machine epsilon is defined as the smallest number that, when added to one, yields a result
//  different from one."
// however, std::numeric_limits<double>::epsilon() returns 2.2204460492503131e-16 which is
// 2^(-52) = 2^(-(N-1))

// 1/x is representable exactly (i.e. with terminatig sequence of digits) if base / x is integer,
// i.e. base % x == 0, so in order to have all fractions with small denominator (and
// their multiples) represented exactly, we should choose a base 2*3*4*5*6*7..., so
// 1/2,1/3,1/4,1/5,1/6,1/7, ... an their multiples have exact representations -> choose the
// smallest factorial <= 2^32 which is 479001600 = 12! - so fractions with numerators up to 12
// would have exact representations (verify this)

// for base = 2, numDigits = 24, the precision should be equivalent to the built-in C++ single
// precision float datatype, numDigits = 53 would correspond to double -> maybe we should test
// this in a unit-test

// write a function that returns the largest integer that is exactly representable with all smaller
// integers also exactly representable. this is exactly base^numDigits. for example, with 2 digits
// in base 10, we can represent integers up to 100. above 100, the stepsize between integers is 10,
// above 1000, it's 100 etc.

// have a static int variable roundingMode which can take values TOWARD_ZERO, TOWARD_INFINITY,
// TO_NEAREST, TO_FLOOR, TO_CEIL - maybe these should be known already in rsBigNumber and a
// parameter roundingModeshould be passed to roundHeadAccordingToTail with default value TO_NEAREST

// implement floor function: the floor of a number can be obtained by setting all digits above
// numDigits - exponent to zero. for negative numbers, one would have to add 1 to the last retained
// digit. if a carry digit occurs, we need to rightshift, prepend the carry digit (which is 1) and
// decrement the exponent

// implement fmod function fmod(x, m): k = floor(x/m); return x - k*m;
// test with examples, like fmod(7.5, 2.3) = 0.6

// for example code for arbitrary precision math, see:
// http://www.mathematik.uni-muenchen.de/~forster/sw/aribas.html
