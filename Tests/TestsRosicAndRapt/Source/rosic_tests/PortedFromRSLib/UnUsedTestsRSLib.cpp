template<class T>
bool checkAlgebraicRules(T a, T b, T c)
{
  T r1, r2, r3, r4;

  r1 = a + b;
  r2 = b + a;
  bool additionCommutative = r1 == r2;

  r1 = a + (b + c);
  r2 = (a + b) + c;
  bool additionAssociative = r1 == r2;

  r1 = a * b;
  r2 = b * a;
  bool multiplicationCommutative = r1 == r2;

  r3 = r1 / a;;
  r4 = r1 / b; 
  bool multiplicationInvertible = (r3 == b) && (r4 == a); 

  r1 = a * (b * c);
  r2 = (a * b) * c;
  bool multiplicationAssociative = r1 == r2;

  r1 = a * (b + c);
  r2 = a*b + a*c;
  bool multiplicationDistributive = r1 == r2;

  r1 = (a + b) / c;
  r2 = a/c + b/c;
  bool divisionDistributive = r1 == r2;

  // multiplicationInvertible and divisionDistributive are false, but that's only numerical 
  // error. maybe use weaker comparisons - check only, if they are numerically close.

  return additionCommutative && additionAssociative
    && multiplicationCommutative && multiplicationInvertible && multiplicationAssociative
    && multiplicationDistributive && divisionDistributive;

  // of course, such informal tests do not "proof" anything. there could be contrived special 
  // cases, where the stuff behaves differently (for example, if the multiplication of two nonzero
  // numbers gives a zero) ...but i think, that should be impossible

  // for commutativity of multiplication, the proof could proceed as follows:
  // P := (a + ib)(c + id) = ac + aid + ibc + i^2bd = (ac - bd) + i(ad + bc)
  // Q := (c + id)(a + ib) = ca + cib + ida + i^2db = (ca - db) + i(cb + da)
  // for P = Q, a sufficient criterion would be that multiplication of the underlying real type
  // is communitative, such that ac = ca, bd = db, bc = cb, ad = da. for the imaginary parts to be
  // equal we additionaly need commutativity of addition such that:
  // x + y = y + x, where we identify x = ad = da and y = bc = cb
  // i guess, all the other algebraic rules could similarly be reduced to agebraic rules
  // of the underlying "real" type (which may indeed itself be actually "complex" itself as long as
  // it satifies the desired algebraic rules). one could the conclude, by induction, that the rules
  // hold for any order K.
}
void nestedComplexOrder2()
{
  // investigate algebraic properties of nested complex numbers:
  rsComplex< rsComplex<double> > a, b, c, r1, r2, r3, r4; 

  // questions: are the algebraic operations equivalent to the algebra over quarternions and
  // (as the nesting level increases), octonions, sedenions, ...?
  // hmm...apparently not: quarternion-multiplication is non-commutative, but the nested-complex 
  // multiplication is. is this interesting for anything?
  // maybe define such 2^K-dimensional numbers as "K-complex" numbers or complex numbers of order K
  // that would be a consistent terminology: usual complex numbers would be order-1 complex numbers
  // and real numbers would be order-0 complex numbers

  // use prime numbers for the elements of the operands (such that the products are all unique):
  a.re.re = 2; b.re.re = 11; c.re.re = 23;
  a.re.im = 3; b.re.im = 13; c.re.im = 19;
  a.im.re = 5; b.im.re = 17; c.im.re = 31;
  a.im.im = 7; b.im.im = 19; c.im.im = 37;

  bool algebraicRulesSatisfied = checkAlgebraicRules(a, b, c);

  // multiplication:
  // prd = ((2+3i) + (5+7i)j)  * ((23+19i) + (31+37i)j)
  //     =  (2+3i)*(23+19i) + (2+3i)*(31+37i)j + (5+7i)j*(23+19i) + (5+7i)j)*(31+37i)j
  //     = 2*23 + 2*19i + 3i*23 + 3i*19i 


  // to do:

  // define a mapping between the K-complex numbers and 2^K vectors: as the vector-element index 
  // increases, the innermost nesting level switches re/im between successive vector elements, the
  // next-to innermost nesting level switches re/im between vector elements that are separated by
  // 2, etc.
  // maybe write the stuff up in a paper "An Algebra for 2^K-dimensional Vector Spaces over the 
  // Real Numbers based on Nested Complex Numbers" 

  // write functions void vectorMultiply/Divide(double *a, double *b, double *result, int N);
  // maybe with a typedef double Real and using "Real" for the parameters. by writing out
  // the multiplication/division operations, it should be possible to conceive an O(N*log(N))
  // algorithm for multiplication and division that operates directly on the vectors.

  // having all the common algebraic operations in place, one could proceed to define the 
  // elementary functions via series expansions. next, one could look at differentiation and 
  // integration - maybe stuff like the Cauchy/Riemann equations for complex differentiability
  // will also generalize and simply "propagate" through the levels of nesting? for 
  // differentiation, we would have to find all the partial derivatives of each vector-element
  // to each other vector element, i.e., a 2^K times 2^K matrix and establish relations between
  // the elements of such a matrix which would generalize the Cauchy/Riemann equations

  // check, how all this relates to vector-analysis, in particular, compare the 2-complex 
  // multiplication (with the 4th element set to zero) to the cross-product in 3 space,
  // and maybe also in higher dimensional spaces (there's a generalization of the cross-product
  // for these - look it up)

  // maybe write a class rsNestedComplex that encapsulates such a nesting with a nestingLevel
  // as member/user-parameter

  // applications: division and multiplication in a vector-space gives rise to an N = 2^K 
  // parametric family of linear transformations (the N elements of the "factor-vector" 
  // being the parameters), some of which could turn out to be useful. question: could linear 
  // transformations be obtained as special cases by requiring some relations between elements of 
  // the factor-vector? the transformations could be computed in O(N*log(N)) time (i think) and 
  // are invertible (via division)
  // maybe one could implement digital filters with such vectors as inputs, outputs and 
  // coefficients

  // maybe, all this stuff could be generalized to vectors of any lengths by just assuming
  // zero elements (i.e. zero padding to the next power of 2) -  we need to check, if K-complex
  // multiplication of two vectors with zeros at some index (in both vectors) will also produce
  // a zero at the same index in the product vector. ...nope - it doesn't, so we have to restrict 
  // ourselves to N = 2^K


  // OK - it's nothing new really:
  // http://reedbeta.wordpress.com/2012/10/15/nested-complex-numbers/
  // http://en.wikipedia.org/wiki/Hypercomplex_numbers
  // http://home.comcast.net/~cmdaven/hyprcplx.htm
}









/** Given (zero-based) indices i, j for a matrix recursively constructed from 4 parameters a,b,c,d
(generalized Hadamard matrix), this functions returns the exponents for a,b,c,d for the 
matrix-element at index i,j in pa, pb, pc, pd respectively */
void findExponentsForRecursiveMatrix(int i, int j, int N, int &pa, int &pb, int &pc, int &pd)
{
  pa = pb = pc = pd = 0;
  while( N > 1 )
  {
    N /= 2;
    if( i < N )
    {
      if( j < N )
        pa++;       // top-left quadrant
      else
        pb++;       // top-right quadrant
    }
    else
    {
      if( j < N )
        pc++;       // bottom-left quadrant
      else
        pd++;       // bottom-right quadrant
    }
    i %= N;
    j %= N;
  }
}

/** Given the index into the product i and the index into the 1st factor j, this function 
computes k the index into the 2nd factor, and the sign s for a term in a product of vectors. */
void indexAndSignForSecondFactor(int i, int j, int N, int &k, double &s)
{
  k = 0;
  s = 1.0;
  while( N > 1 )
  {
    N /= 2;
    if( i < N )
    {
      if( j >= N )
      {
        k += N;       // top-right quadrant
        s *= -1.0; 
      }
    }
    else
    {
      if( j < N )
        k += N;       // bottom-left quadrant
    }
    i %= N;
    j %= N;
  }
}
void vectorMultiply(double *a, double *b, double *p, int N)
{
  int k;
  double sign;
  for(int i = 0; i < N; i++)
  {
    p[i] = 0.0;
    for(int j = 0; j < N; j++)
    {
      indexAndSignForSecondFactor(i, j, N, k, sign);
      p[i] += sign * a[j] * b[k];
    }
  }
}
void vectorDivide(double *a, double *b, double *p, int N)
{
  int k;
  double sign;
  double tmp;
  rsCopyBuffer(p, a, N);
  for(int i = 0; i < N; i++)
  {
    for(int j = 0; j < N; j++)
    {
      indexAndSignForSecondFactor(i, j, N, k, sign);
      tmp = p[i] - sign * a[j] * b[k]; 
      if( b[k] != 0.0 )
      {
        a[j] = (p[i] - tmp) / (sign * b[k]);
      }
      p[i] = tmp;
    }
  }
  // for division, run the multiplication algorithm backwards: initialize a with the product p and 
  // and in the inner loop, do: a[j] -= sign * p[i] / b[k]; (if b[k] is zero, just skip the step - 
  // this is justified by observing, that in the corresponding step of the multiplication, a 
  // zero-value was added into the product - which we don't need to subtract)
  // ..hmmm - this still doesn't work

  // the p[i] += sign * a[j] * b[k]; instruction in the multiplication can be rewritten as 
  // equation: p'[i] = p[i] + sign * a[j] * b[k]; where the prime denotes the new, updated value.
  // we may, solve this equation for a[j] but then, the solution would involve a p[i] (the 
  // non-updated value) which is unknown at this point. so, as it seems, simply running 
  // multiplication backwards doesn't work. perhaps, we should write a recursive multiplication
  // procedure and see, if that's easier to invert (or otherwise convert)

  // complex division is actually a multiplication by the complex conjugate followed by a division
  // by the (real) sum of the squares of the elements of the denominator. maybe this generalizes
  // directly: invert some signs (which ones), do the multiplciation and finally normalize by
  // squared norm of the vector

  // define the K-complex conjugate in such a way, that multiplying a K-complex number with its
  // K-complex conjugate gives a real number (sum of the squares of the elements - or maybe not?
  // maybe, there are cross-terms as well?):
  // (a+bi)  (a-bi) = a^2 - abi + abi - b^2i^2 = a^2 + b^2
  //
  // consider the 2-complex number: ((a+bi)+(c+di)j). let's eliminate j by multiplying with a
  // number, where + is replaced by - in the outer level, i.e. with ((a+bi)-(c+di)j):
  //
  // ((a+bi)+(c+di)j) * ((a+bi)-(c+di)j) 
  // = -d^2i^2j^2 - 2cdij^2 - c^2j^2 + b^2i^2 + 2abi + a^2
  // = -d^2 + 2cdi + c^2 - b^2 + 2abi + a^2
  // = a^2 - b^2 + c^2 - d^2 + 2abi + 2cdi
  // = (a^2-b^2+c^2-d^2) + 2(ab+cd)i
  //
  // next, eliminate i by multiplying the result by (a^2-b^2+c^2-d^2) - 2(ab+cd)i:
  //
  // ((a^2-b^2+c^2-d^2) + 2(ab+cd)i) * ((a^2-b^2+c^2-d^2) - 2(ab+cd)i)
  // = -4c^2d^2i^2-8abcdi^2-4a^2b^2i^2+d^4-2c^2d^2+2b^2d^2-2a^2d^2+c^4-2b^2c^2+2a^2c^2+b^4-2a^2b^2+a^4
  // = +4c^2d^2 + 8abcd + 4a^2b^2 + d^4 - 2c^2d^2 + 2b^2d^2 - 2a^2d^2 + c^4 - 2b^2c^2 + 2a^2c^2 + b^4 - 2a^2b^2 + a^4
  // = 2c^2d^2 + 8abcd + 2a^2b^2 + d^4 + 2b^2d^2 - 2a^2d^2 + c^4 - 2b^2c^2 + 2a^2c^2 + b^4 + a^4
  // = a^4 + b^4 + c^4 + d^4 + 2(a^2b^2 + a^2c^2 + b^2d^2 + c^2d^2 - a^2d^2 - b^2c^2) + 8abcd
  //
  // hmm... it seems to be more convenient to multiply with ((a-bi)-(c-di)j), i.e., a number,
  // where all + are turned into -, this eliminates i and j in one go, but introduces ij-terms
  // which have to be eliminated in the next step (but these are simpler):
  //
  // ((a+bi)+(c+di)j) * ((a-bi)-(c-di)j) 
  // = d^2i^2j^2-c^2j^2+2adij-2bcij-b^2i^2+a^2
  // = d^2i^2j^2 - c^2j^2 + 2adij - 2bcij - b^2i^2 + a^2
  // = a^2 + b^2 + c^2 + d^2 + 2adij - 2bcij
  // = a^2+b^2+c^2+d^2 + 2(ad-bc)ij
  // = s + 2(ad-bc)ij  where s := a^2+b^2+c^2+d^2
  //
  // next, eliminate the ij-term by multiplying with s - 2(ad-bc)ij
  //
  // (s + 2(ad-bc)ij) * (s - 2(ad-bc)ij)
  // = s^2 - 4a^2d^2i^2j^2 + 8abcdi^2j^2 - 4b^2c^2i^2j^2
  // = s^2 - 4(a^2d^2 - b^2c^2) + 8abcd
  // = (a^2+b^2+c^2+d^2)^2 - 4(a^2d^2-b^2c^2) + 8abcd
  //
  // so, the 2-complex conjugate to ((a+bi)+(c+di)j) is the product of the two numbers by which
  // we multiplied in the 2 steps:
  //
  // ((a-bi)-(c-di)j) * (a^2+b^2+c^2+d^2 - 2(ad-bc)ij)
  // = -2ad^2i^2j^2+2bcdi^2j^2+2acdij^2-2bc^2ij^2+2abdi^2j-2b^2ci^2j+d^3ij+c^2dij+b^2dij-a^2dij+2abcij-cd^2j-c^3j-b^2cj-a^2cj-bd^2i-bc^2i-b^3i-a^2bi+ad^2+ac^2+ab^2+a^3
  // = -2ad^2 + 2bcd - 2acdi + 2bc^2i - 2abdj + 2b^2cj + d^3ij + c^2dij + b^2dij - a^2dij + 2abcij - cd^2j - c^3j - b^2cj - a^2cj - bd^2i - bc^2i - b^3i - a^2bi + ad^2 + ac^2 + ab^2 + a^3
  //
  // collect terms contain none, i, j, and ij:
  //
  //      : -2ad^2 + 2bcd  + ad^2 + ac^2 + ab^2 + a^3
  //     i: - 2acd + 2bc^2 - bd^2 - bc^2 - b^3 - a^2b
  //    j : - 2abd + 2b^2c - cd^2 - c^3 - b^2c - a^2c
  //    ji: + d^3 + c^2d + b^2d - a^2d + 2abc
  //
  // replace a^2 with aa, etc. and rearranging (for seeing paaterns):
  //
  //      : +aaa +bba +cca -dda   + 2 bcd        + + + -   +
  //     i: -aab -bbb +ccb -ddb   - 2a cd        - - + -   -
  //    j : -aac +bbc -ccc -ddc   - 2ab d        - + - -   -
  //    ji: -aad +bbd +ccd +ddd   + 2abc         - + + +   +
  //
  // pattern: in every column there's only one sign different from the others, and that's
  // the one for the cubed-term - ao, from +aaa, we can conclude -aad, from which we conclude +ddd
  // ...maybe, we can find an algrithm
  //
  // so, at the end of the day, multiplying ((a+bi)+(c+di)j) by
  //   (+aaa+bba+cca-dda+2bcd)  + (-aab-bbb+ccb-ddb-2acd)i 
  // + (-aac+bbc-ccc-ddc-2abd)j + (-aad+bbd+ccd+ddd+2abc)ji
  //
  // should result in a real number, i.e. a number without any i,j tags 
  // -> verify with maxima and/or numerically
  //
  // throwing this at maxima:
  //
  // ((a+b*i)+(c+d*i)*j) * ((+a*a*a+b*b*a+c*c*a-d*d*a+2*b*c*d) + (-a*a*b-b*b*b+c*c*b-d*d*b-2*a*c*d)*i + (-a*a*c+b*b*c-c*c*c-d*d*c-2*a*b*d)*j + (-a*a*d+b*b*d+c*c*d+d*d*d+2*a*b*c)*j*i)
  //
  // results in (already simplified notationally):
  //
  // d^4+c^2d^2+b^2d^2-a^2d^2+2abcd-2abd^2ij^2+2b^2cdij^2-2a^2cdij^2+2abc^2ij^2-c^2d^2j^2-2abcdj^2-c^4j^2+b^2c^2j^2-a^2c^2j^2-2acd^2i^2j+2bc^2di^2j-2a^2bdi^2j+2ab^2ci^2j-2acd^2j+2bc^2dj-2a^2bdj+2ab^2cj-b^2d^2i^2-2abcdi^2+b^2c^2i^2-b^4i^2-a^2b^2i^2-2abd^2i+2b^2cdi-2a^2cdi+2abc^2i-a^2d^2+2abcd+a^2c^2+a^2b^2+a^4
  // d^4+ccdd+bbdd-aadd+2abcd-2abddij^2+2bbcdij^2-2aacdij^2+2abccij^2-ccddj^2-2abcdj^2-c^4j^2+bbccj^2-aaccj^2-2acddi^2j+2bccdi^2j-2aabdi^2j+2abbci^2j-2acddj+2bccdj-2aabdj+2abbcj-bbddi^2-2abcdi^2+bbcci^2-b^4i^2-aabbi^2-2abddi+2bbcdi-2aacdi+2abcci-aadd+2abcd+aacc+aabb+a^4
  // d^4+ccdd+bbdd-aadd+2abcd+2abddi-2bbcdi+2aacdi-2abcci+ccdd+2abcd+c^4-bbcc+aacc+2acddj-2bccdj+2aabdj-2abbcj-2acddj+2bccdj-2aabdj+2abbcj+bbdd+2abcd-bbcc+b^4+aabb-2abddi+2bbcdi-2aacdi+2abcci-aadd+2abcd+aacc+aabb+a^4
  // 
  // collect i,j terms and those without tags:
  //
  //  : d^4+ccdd+bbdd-aadd+2abcd+ccdd+2abcd+c^4-bbcc+aacc+bbdd+2abcd-bbcc+b^4+aabb-aadd+2abcd+aacc+aabb+a^4
  // i: +2abddi-2abddi+2aacdi-2aacdi+2bbcdi-2bbcdi+2abcci-2abcci = 0
  // j: +2acddj-2acddj+2bccdj-2bccdj+2aabdj-2aabdj+2abbcj-2abbcj = 0
  //
  // re-arrange the term without tags:
  //
  // +aaaa+aabb+aabb+aacc+aacc-aadd-aadd +bbbb-bbcc-bbcc+bbdd+bbdd +cccc+ccdd+ccdd + dddd +2abcd+2abcd+2abcd+2abcd
  // +aa(aa+bb+bb+cc+cc-dd-dd) +bb(bb-cc-cc+dd+dd) +cc(cc+dd+dd) + dddd +8abcd
  // +aa(aa+2bb+2cc-2dd) +bb(bb-2cc+2dd) +cc(cc+2dd) +dddd +8abcd
  //
  // the last line should be the real number that results, when multiplying ((a+b*i)+(c+d*i)*j) 
  // with its 2-complex conjugate
  //
  // ((a+b*i)+(c+d*i)*j) = (a + bi + cj + dij)


}
typedef rsComplex<rsComplex<double> > rsComplexOrder2;
typedef rsComplex<rsComplex<rsComplex<double> > > rsComplexOrder3;
typedef rsComplex<rsComplex<rsComplex<rsComplex<double> > > > rsComplexOrder4; // doesn't work
void nestedComplexOrder3()
{
  rsComplexOrder3 ac, bc, cc, r1, r2, r3, r4;  // nested complex numbers
  double av[8], bv[8], cv[8], pv[8], qv[8];    // numbers, mapped into vectors
  double a,b,c,d,e,f,g,h;                      // elements of a
  double p,q,r,s,t,u,v,w;                      // elements of b

                                               // imaginary tags
  a = av[0] = ac.re.re.re =  2;   // 
  b = av[1] = ac.re.re.im =  3;   //   i
  c = av[2] = ac.re.im.re =  5;   //   j
  d = av[3] = ac.re.im.im =  7;   //  ji
  e = av[4] = ac.im.re.re = 11;   // k
  f = av[5] = ac.im.re.im = 13;   // k i
  g = av[6] = ac.im.im.re = 17;   // kj
  h = av[7] = ac.im.im.im = 19;   // kji

  p = bv[0] = bc.re.re.re = 23;
  q = bv[1] = bc.re.re.im = 29;
  r = bv[2] = bc.re.im.re = 31;
  s = bv[3] = bc.re.im.im = 37; 
  t = bv[4] = bc.im.re.re = 41;
  u = bv[5] = bc.im.re.im = 43;
  v = bv[6] = bc.im.im.re = 47;
  w = bv[7] = bc.im.im.im = 53; 

  // we stil have x,y,z free, so maybe rename ac, bc, cc into x,y,z

  cv[0] = cc.re.re.re = 59;
  cv[1] = cc.re.re.im = 61;
  cv[2] = cc.re.im.re = 67;
  cv[3] = cc.re.im.im = 71; 
  cv[4] = cc.im.re.re = 73;
  cv[5] = cc.im.re.im = 79;
  cv[6] = cc.im.im.re = 83;
  cv[7] = cc.im.im.im = 89; 

  bool algebraicRulesSatisfied = checkAlgebraicRules(ac, bc, cc);


  cc = ac * bc;

  cv[0] = cc.re.re.re;
  cv[1] = cc.re.re.im;
  cv[2] = cc.re.im.re;
  cv[3] = cc.re.im.im; 
  cv[4] = cc.im.re.re;
  cv[5] = cc.im.re.im;
  cv[6] = cc.im.im.re;
  cv[7] = cc.im.im.im; 

  vectorMultiply(av, bv, pv, 8);

  bool testResult = true;
  for(int i = 0; i < 8; i++)
    testResult &= cv[i] == pv[i];


  vectorDivide(qv, bv, pv, 8);



  int dummy = 0;



  // ac = (2+3i) + (5+7i)i + ((11+13i) + (17+19i)i)i
  //      (a+bi) + (c+di)i + ((e + fi) + (g + hi)i)i
  // maybe define a dimensionality reduction operation multiplying the inner terms out
  // but maybe it's necessarry to distinguish the "imaginary units" of the different levels
  // if we denote them by i,j,k, we have:
  // ac = (2 + 3i) + (5 + 7i)j + ((11+13i) + (17+19i)j)k
  //      (a + bi) + (c + di)j + ((e + fi) + (g + hi)j)k       and similary:
  // bc = (23+29i) + (31+37i)j + ((41+43i) + (47+53i)j)k
  //      (p + qi) + (r + si)j + ((t + ui) + (v + wi)j)k
  // i is the innermost "imaginary unit" and k is the outermost
  // the product of ac*bc is:
  // prd = ((a+bi)+(c+di)j+((e+fi)+(g+hi)j)k) * ((p+qi)+(r+si)j+((t+ui)+(v+wi)j)k)

  // at the even indices we have terms that do not contain i, at the odd indices, we have terms
  // that do contain i, at indices for which floor(index/2) is even, we don't have j, when index/2
  // is odd, we have a j, at indices for which floor(index/4) is even, we have no k, when it's odd
  // we have a k

  // throwing:
  // ((a+b*i)+(c+d*i)*j+((e+f*i)+(g+h*i)*j)*k) * ((p+q*i)+(r+s*i)*j+((t+u*i)+(v+w*i)*j)*k);
  // at maxima and using factor(%) results in:
  // h*k^2*w+g*i*j^2*k^2*w+f*i^2*j*k^2*w+ ...mayn more terms
  //
  // collect terms an sort by 1st letter:
  //
  //    :   ap - bq - cr + ds - et + fu + gv - hw 
  //   i:   aq + bp - cs - dr - eu - ft + gw + hv
  //  j :   ar - bs + cp - dq - ev + fw - gt + hu
  //  ji:   as + br + cq + dp - ew - fv - gu - ht
  // k  :   at - bu - cv + dw + ep - fq - gr + hs
  // k i:   au + bt - cw - dv + eq + fp - gs - hr
  // kj :   av - bw + ct - du + er - fs + gp - hq
  // kji:   aw + bv + cu + dt + es + fr + gq + hp
  //
  // writing [a b c d e f g h] X [p q r s t u v w]  as matrix/vector product:
  //
  // [a b c d e f g] [+p -q -r +s -t +u +v -w]
  // [a b c d e f g] [+q +p -s -r -u -t +w +v]
  // [a b c d e f g] [+r -s +p -q -v +w -t +u]
  // [a b c d e f g] [+s +r +q +p -w -v -u -t]
  // [a b c d e f g] [+t -u -v +w +p -q -r +s]
  // [a b c d e f g] [+u +t -w -v +q +p -s -r]
  // [a b c d e f g] [+v -w +t -u +r -s +p -q]
  // [a b c d e f g] [+w +v +u +t +s +r +q +p]
  //
  // permutation pattern (let n be the index into the a-vector, m the index into the b-vector):
  // 
  // 0 1 2 3 4 5 6 7  n = 0:      m = n
  // 1 0 3 2 5 4 7 6  n = N/4-1:  m = N/4-1-n (1st quarter), N/4-1-n+N/4 (2nd), N/4-1-n+2N/4 (3rd), ...
  // 2 3 0 1 6 7 4 5
  // 3 2 1 0 7 6 5 4  n = N/2-1:  m = N/2-1-n for 1st half, m =  N/2-1-n+N/2 for 2nd half
  // 4 5 6 7 0 1 2 3  n = N/2:    m = N/2+n   for 1st half, m = -N/2+n for 2nd half
  // 5 4 7 6 1 0 3 2
  // 6 7 4 5 2 3 0 1
  // 7 6 5 4 3 2 1 0  n = N-1:    m = N-1-n    (last row is reversal of 1st)
  //
  // in the top-left and bottom-right quadrants, only the lower values occur (0...N/2-1) and in the
  // top-right and bottom-left quadrants, only the higher values occur (N/2...N-1)
  //
  // sign pattern:  
  //
  // + - - + - + + -  2nd half is negative of 1st
  // + + - - - - + +  same here
  // + - + - - + - +  same here    
  // + + + + - - - -  same here
  // + - - + + - - +  2nd half is same as 1st
  // + + - - + + - -  same here
  // + - + - + - + -  same here
  // + + + + + + + +  same here
  //
  // looks like the +/- pattern can be constructed in a similar way as hadamard matrices:
  //
  // H_n+1 = H_n -H_n
  //         H_n  H_n

  // when multiplying two 3-complex numbers (K=3, N=8), the outermost level spawns 4 calls to
  // the "*" operator of the intermediate level, which each spawn 4 calls of the "*" operator
  // of the lowes level.
  // call:   1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21
  // level:  3 2 1 1 1 1 2 1 1  1  1  2  1  1  1  1  2  1  1  1  1
  // in general, we will see sum_n=k^{K-1} 4^k = (4^K-1)/3 invoacations of the "*" operator
  // the actual real-number multiplications are done in the innermost 4^(K-1) calls and each such
  // call performs 4 real mutliplications, totaling to 4 * 4^(K-1) = 4^K multiplications.
  // with K = log2(N), we see that there are a total number of N^2 real multiplications
  // ...we can forget about looking for an N*log(N) algorithm for multiplications, since these 
  // products all involve possibly distinct factors (in a*b, a (real) product of each 
  // vector-element of a with each vector element of b is computed). when the product 
  // c = a*b = A*b is seen as a matrix/vector product with matrix A, the rows of A would have
  // to be permutations of the vector a (possibly with some elements sign-inverted), and none
  // of these permutations can have any element in the same position as any other (otherwise some
  // of the real products would be computed twice and some other not at all) - we need to find out
  // how exactly these permutations are constructed
  //
  // with only one level of nesting:
  //
  // ((a + b*i) + (c + d*i)*j) * ((p + q*i) + (r + s*i)*j)
  // (a+bi)(b+qi) + (a+bi)(r+si)j + (c+di)(p+qj)j
  //
  // it's managable:
  //
  // prd = d*s+c*i*j^2*s+b*i^2*j*s+a*i*j*s+d*i*j^2*r+c*j^2*r+b*i*j*r+a*j*r+d*i^2*j*q
  //       +c*i*j*q+b*i^2*q+a*i*q+d*i*j*p+c*j*p+b*i*p+a*p
  //
  //
  //   prd =  (ap-bq) + (ds-cr)
  //          (aq+bp) - (dr+cs)   | *i
  //          (ar-bs) + (cp-dq)   |   *j
  //          (as+br) + (cq+dp)   | *i*j
  //
  // i'm not sure, if that's the correct order. it's probably best to just run the multiplication 
  // through the debugger and observe, in which order things are computed. we have 4 different
  // variables, and all their 4^2=16 possible products occur, which suggests, that the computation
  // would require O(N^2) operations. but probably some multiplications can be factored out by 
  // re-arranging the computations appropriately, for 
  //
  // t0=ap, t1=aq, t2=bp, t3=bq, p0r=ap-bq, p0i=aq+bp  // temporaries and product re/im, products
  // t0=cp, t1=cq, t2=dp, t3=dq, p1r=cp-dq, p1i=cq+dp 
  // t0=ar, t1=as, t2=br, t3=bs, p2r=ar-bs, p2i=as+br
  // t0=cr, t1=cs, t2=dr, t3=ds, p3r=cr-ds, p3i=cs+dr
  //                         ... pnr=t0-t3, pni=t1+t2  for general index n
  //
  // we currently have, as intermediate result:
  //   prd = ap-bq  
  //         aq+bp    |  i
  //         ar-bs    | j
  //         as+br    | ji
  //
  // in this intermediate result, p0r, p0i, p1r, etc. represent the re/im parts products of the 
  // inner nesting level for the left operand (of the outer level). to that we must add the re/im
  // parts of the products of the inner level of the right operand (of the oputer level):


  dummy = 0;
}