using namespace RSLib;

// construction/destruction:

rsBigNumber::~rsBigNumber()
{
  freeMemory();
}

// operators:


// setup:


// inquiry:


// misc:

rsUint32 rsBigNumber::rsCharToUint32(char c)
{
  if( c >= 48 && c <= 57 )
    return c - 48; // '0'...'9' map to 0...9
  else if( c >= 65 && c <= 90 )
    return c - 55; // 'A'...'Z' map to 10...35
  else if( c >= 97 && c <= 122 )
    return c - 87; // 'a'...'z' map to 10...35
  rsError("Character does not represent a digit.");
  return 0;
}

char rsBigNumber::rsUint32ToChar(rsUint32 x)
{
  if( x <= 9 )
    return x + 48;
  else if( x <= 35 )
    return x + 55;
  rsError("Number not representable by single character");
  return 'X';
}

void rsBigNumber::allocateMemory()
{
  digits = new rsUint32[numDigits];
}

void rsBigNumber::freeMemory()
{
  delete[] digits;
}

// arithmetic:

rsUint32 rsBigNumber::add(rsUint32 *a, rsUint32 *b, rsUint32 N, rsUint32 *s, rsUint64 base)
{
  rsUint64 tmp;
  rsUint64 carry = 0;
  for(int i = N-1; i >= 0; i--)
  {
    tmp   = (rsUint64)a[i] + (rsUint64)b[i] + carry;
    carry = tmp >= base;                     // == tmp / base
    s[i]  = (rsUint32) (tmp - carry * base); // == tmp % base
  }
  return (rsUint32)carry; 
}

rsUint32 rsBigNumber::add(rsUint32 *a, rsUint64 b, rsUint32 N, rsUint32 *s, rsUint64 base)
{
  rsUint64 tmp;
  rsUint64 carry = b;
  for(int i = N-1; i >= 0; i--)
  {
    tmp   = (rsUint64)a[i] + carry;
    carry = tmp >= base;                     // == tmp / base
    s[i]  = (rsUint32) (tmp - carry * base); // == tmp % base
  }
  return (rsUint32)carry; 
}

rsUint32 rsBigNumber::subtract(rsUint32 *a, rsUint32 *b, rsUint32 N, rsUint32 *d, rsUint64 base)
{
  rsInt64 tmp;
  rsInt64 borrow = 0;
  for(int i = N-1; i >= 0; i--)
  {
    tmp    = (rsInt64)a[i] - (rsInt64)b[i] - borrow; 
    borrow = tmp < 0;
    d[i]   = (rsUint32)(tmp + borrow * base);
  }
  return (rsUint32)borrow;
}

rsUint32 rsBigNumber::subtract(rsUint32 *a, rsUint32 b, rsUint32 N, rsUint32 *d, rsUint64 base)
{
  rsInt64 tmp;
  rsInt64 borrow = b;
  for(int i = N-1; i >= 0; i--)
  {
    tmp    = (rsInt64)a[i] - borrow; 
    borrow = tmp < 0;
    d[i]   = (rsUint32)(tmp + borrow * base);
  }
  return (rsUint32)borrow;
}
   
void rsBigNumber::multiply(rsUint32 *a, rsUint32 *b, rsUint32 N, rsUint32 *p, rsUint64 base)
{
  multiply(a, N, b, N, p, base);
}

void rsBigNumber::multiply(rsUint32 *a, rsUint32 Na, rsUint32 *b, rsUint32 Nb, rsUint32 *p, 
                           rsUint64 base)
{
  rsUint64 tmp;
  rsUint32 carry;
  rsFillWithZeros(&p[Nb], Na); 
  for(int j = Nb - 1; j >= 0; j--)
  {
    carry = 0;
    for(int i = Na - 1; i >= 0; i--)
    {
      tmp      = (rsUint64)a[i] * (rsUint64)b[j] + (rsUint64)p[i+j+1] + (rsUint64)carry;
      p[i+j+1] = (rsUint32) (tmp % base);
      carry    = (rsUint32) (tmp / base);
    }
    p[j] = carry;
  }
}

rsUint32 rsBigNumber::multiply(rsUint32 *a, rsUint32 N, rsUint64 b, rsUint32 *p, rsUint64 base)
{  
  rsUint64 tmp;
  rsUint32 carry = 0;
  for(int i = N-1; i >= 0; i--)
  {
    tmp   = b * (rsUint64)a[i] + (rsUint64)carry; 
    carry = (rsUint32) (tmp / base);
    p[i]  = (rsUint32) (tmp % base);
  }
  return carry;
}

void rsBigNumber::divide(rsUint32 *a, rsUint32 *b, rsUint32 N, rsUint32 *q, rsUint32 *r, 
                         rsUint64 base)
{
  // the algorithm is from (2), page 257-258 with a modified indexing scheme for u[]: instead of 
  // using u[j], we use u[0] and leftshift the u-array after each completion of the loop such that 
  // we don't need to allocate twice as much memory for the array.

  // the general algorithm fails in the special case N == 1, handle it separately:
  if( N == 1 )
  {
    rsUint32 tmp = a[0] / b[0];
    r[0]         = a[0] % b[0];
    q[0]         = tmp;
    return;
  }

  // if a < b, the quotient is zero and the remainder is a:
  if( rsCompare(a, b, N) == -1 )
  {    
    rsCopyBuffer(a, r, N);
    rsFillWithZeros(q, N);
    return;
  }

  // arrays for normalized dividend, divisor and weighted (normalized) divisor
  rsUint32 *u = new rsUint32[N+1]; rsFillWithZeros(u, N+1);
  rsUint32 *v = new rsUint32[N+1]; rsFillWithZeros(v, N+1);
  rsUint32 *w = new rsUint32[N+1]; rsFillWithZeros(w, N+1);
    // \todo: optimization: allocate a single "workArea" for size 3*(N+1) - but measure, if that's
    // really faster

  // find the number of non-leading-zero places of the remainder and quotient:
  int ia = rsFirstIndexWithNonZeroValue(a, N);
  int ib = rsFirstIndexWithNonZeroValue(b, N);
  int  n = N - ib;
  int  m = N - n - ia;

  // D1: "normalize" - multiply dividend and divisor by a 1-digit number d that ensures that the 
  // divisor's 1st digit is reasonably large (this leaves the overall quotient unchanged and 
  // multiplies the remainder by the same number):
  rsUint32 d = (rsUint32) (base / (b[ib] + 1));
  u[0] = multiply(&a[ia], N-ia, d, &u[1], base);
  v[0] = multiply(&b[ib], N-ib, d, &v[1], base); // v[0] should always be 0

  // initialize pointer into the quotient, where the actual nonzero-digits go and prepend an
  // appropriate number of zeros:
  rsUint32 *pq = &q[N-m-1];
  rsFillWithZeros(q, N-m-1);

  // main loop - compute m+1 digits for quotient:
  for(int j = 0; j <= m; j++) // D2, D7 - init and check loop conditions
  {
    // D3: calculate trial value for quotient digit pq[j] - if the divisor's 1st digit is 
    // reasonably large (which we have assured above), this heuristic is pretty good:
    if( u[0] == v[1] )
      pq[j] = (rsUint32) (base - 1);
    else
      pq[j] = (rsUint32) ((u[0]*base + u[1]) / v[1]);

    // eliminate all cases where trial value pq[j] is two too large and most cases where it is one 
    // too large:
    while( (rsUint64)v[2]*pq[j] > (u[0]*base + u[1] - (rsUint64)v[1]*pq[j]) * base + u[2] )
      pq[j]--;

    // remark: when both operands of a multiplication are of type rsUint32 (like v[2]*pq[j]), we 
    // use an explicit typecast to rsUint64 on one of the operands to avoid overflow of the
    // multiplication result. the other operand will be casted automatically. it seems to work 
    // without, though - maybe there are some conditions which ensure that such overflow never
    // occurs. but i'm not sure. but then, i tested only base 10 in which case we don't expect 
    // overflow anyway

    // D4: multiply v with trial quotient digit pq[j] and subtract result w from u:
    w[0]            = multiply(&v[1], n, pq[j], &w[1], base); 
    rsUint32 borrow = subtract(&u[0], w, n+1,   &u[0], base);

    // D6: trial value pq[j] was one too high - add back one v to u and decrease pq[j] by one:
    if( borrow != 0 )
    {
      pq[j]--;
      add(&u[0], v, n+1, &u[0], base);
    }

    rsLeftShift(u, N+1, 1); // compensates for our different indexing
  }

  // obtain actual remainder by dividing the "normalized" remainder (which is what remains in u) by 
  // the normalization factor d, write it into the appropriate places in r and prepend zeros:
  divide(u, n, d, &r[N-n], base);
  rsFillWithZeros(r, N-n);
    // todo: wrap this into a conditional if( r != nullptr ), so we may not touch the buffer at all
    // -> useful when no remainder is desired and we don't want to allocate memory for it

  delete[] u;
  delete[] v;
  delete[] w;
}

rsUint32 rsBigNumber::divide(rsUint32 *a, rsUint32 N, rsUint64 b, rsUint32 *q, rsUint64 base)
{
  rsUint64 tmp;
  rsUint32 rem = 0;
  for(rsUint32 j = 0; j < N; j++)
  {
    tmp  = (rsUint64)a[j] + (rsUint64)rem * base;
    q[j] = (rsUint32) (tmp / b);
    rem  = (rsUint32) (tmp % b);
  }
  return rem;
}

bool rsBigNumber::divideFractions(rsUint32 *a, rsUint32 Na, rsUint32 *b, rsUint32 Nb, 
                                  rsUint32 *q, rsUint32 Nq, rsUint64 base)
{
  rsAssert( Nq >= Na+Nb );

  // (in-place) integer division on zero-extended operand arrays:
  rsUint32 *r = new rsUint32[Nq];
  rsCopyBuffer(a,  q,        Na); rsFillWithZeros(&q[Na], Nq-Na);
  rsCopyBuffer(b, &r[Nq-Nb], Nb); rsFillWithZeros( r,     Nq-Nb);
  divide(q, r, Nq, q, r, base);

  // round upward, if 2*remainder > b:  
  bool roundUp = false;
  r[Nq-Nb-1] = multiply(&r[Nq-Nb], Nb, 2, &r[Nq-Nb], base);
  if( r[Nq-Nb-1] > 0 )
    roundUp = true;
  else
  {
    int d = rsCompare(&r[Nq-Nb], b, Nb);
    if( d == 1 )
      roundUp = true;
    else if( d == 0 )
      roundUp = shouldRoundUpAtMidway(q[Nq-1], base);
  }

  delete[] r;

  if( roundUp )
    return add(q, 1, Nq, q, base) > 0;
  else
    return false;
}

bool rsBigNumber::divideFractions(rsUint32 *a, rsUint32 *b, rsUint32 *q, rsUint32 N, 
                                  rsUint64 base)
{
  return divideFractions(a, N, b, N, q, 2*N, base);
}

// conversion:

void rsBigNumber::digitsFromInt(rsUint64 value, rsUint32 *digits, rsUint32 N, rsUint64 base)
{
  for(int i = N-1; i >= 0; i--)
  {
    digits[i]  = (rsUint32) (value % base);
    value     /= base;
  }
}

rsUint64 rsBigNumber::digitsToInt(rsUint32 *digits, rsUint32 N, rsUint64 base)
{
  rsUint64 m = 1;
  rsUint64 r = 0;
  for(int i = N-1; i >= 0; i--)
  {
    r += m * digits[i];
    m *= base;
  }
  return r;
}

rsArray<rsUint32> rsBigNumber::changeBaseForInteger(rsUint32 *digits, rsUint32 numDigits, 
                                                    rsUint64 oldBase, rsUint64 newBase)
{
  rsUint32 *q = new rsUint32[numDigits];
  rsCopyBuffer(digits, q, numDigits);
  rsArray<rsUint32> result;
  while( !rsIsAllZeros(q, numDigits) )
    result.appendElement(divide(q, numDigits, newBase, q, oldBase));
  result.reverse();
  delete[] q;
  return result;

  // \todo if we could determine an upper bound for the required number of digits in the new base,
  // we could directly allocate enough memory, maybe something like:
  // newNumDigits <= oldNumDigits * ceil(oldBase/newBase)? ..or:
  // newNumDigits <= oldNumDigits * numDigitsForOldBaseInNewBase? ...no idea yet
}

// rounding:

int rsBigNumber::compareToOneHalf(rsUint32 *digits, rsUint32 numDigits, rsUint64 base)
{
  rsUint32 thresh = (rsUint32) (base/2); 
  if( rsIsEven(base) )
  {
    if( digits[0] > thresh )
      return +1;
    else if( digits[0] < thresh )
      return -1;
    else
    {
      for(rsUint32 i = 1; i < numDigits; i++)
      {
        if( digits[i] > 0 )
          return +1;
      }  
      return 0;
    }
  }
  else
  {
    for(rsUint32 i = 0; i < numDigits; i++)
    {
      if( digits[i] > thresh )
        return +1;
      else if( digits[i] < thresh )
        return -1;
    }
    return -1;
  }
}

bool rsBigNumber::shouldRoundUpAtMidway(rsUint32 lastHeadDigit, rsUint64 base)
{
  if( rsIsEven(base/2) )
  {
    if( rsIsEven(lastHeadDigit) )
      return true; // code not tested
  }
  else
  {
    if( rsIsOdd(lastHeadDigit) )
      return true;
  }
  return false;
}

bool rsBigNumber::upwardRoundRequired(rsUint32 *tail, rsUint32 tailLength, rsUint64 tailBase,
                                      rsUint32 lastHeadDigit, rsUint64 headBase)
{
  int d = compareToOneHalf(tail, tailLength, tailBase);
  if( d == 0 )
    return shouldRoundUpAtMidway(lastHeadDigit, headBase);
  else
    return d == 1;
}

bool rsBigNumber::roundHeadAccordingToTail(rsUint32* head, rsUint32 headLength, rsUint64 headBase, 
                                           rsUint32 *tail, rsUint32 tailLength, rsUint64 tailBase)
{
  if( upwardRoundRequired(tail, tailLength, tailBase, head[headLength-1], headBase) )
    return add(head, 1, headLength, head, headBase) > 0;
  else
    return false;
}

bool rsBigNumber::round(rsUint32 *digits, rsUint32 oldNumDigits, rsUint32 newNumDigits, 
                        rsUint64 base, bool zeroOutTail)
{
  bool result = roundHeadAccordingToTail( digits,               newNumDigits,              base,
                                         &digits[newNumDigits], oldNumDigits-newNumDigits, base);
  if( zeroOutTail )
    rsFillWithZeros(&digits[newNumDigits], oldNumDigits-newNumDigits);
  return result;
}

bool rsBigNumber::changeBaseForFraction(rsUint32 *oldDigits, rsUint32 oldNumDigits, 
  rsUint64 oldBase, rsUint32 *newDigits, rsUint32 newNumDigits, rsUint64 newBase)
{
  rsUint32 *tmpDigits = new rsUint32[oldNumDigits];     
  rsCopyBuffer(oldDigits, tmpDigits, oldNumDigits);
  bool roundOverflow = false;
  rsUint32 i = 0;
  while( i < newNumDigits )
  {
    newDigits[i] = multiply(tmpDigits, oldNumDigits, newBase, tmpDigits, oldBase);
    i++;
    if( rsIsAllZeros(tmpDigits, oldNumDigits) )
      break;
  }
  if( i < newNumDigits ) 
    rsFillWithZeros(&newDigits[i], newNumDigits-i); // sequence terminated early
  else
    roundOverflow = roundHeadAccordingToTail(newDigits, newNumDigits, newBase, 
                                             tmpDigits, oldNumDigits, oldBase);
  delete[] tmpDigits;
  return roundOverflow;
}
