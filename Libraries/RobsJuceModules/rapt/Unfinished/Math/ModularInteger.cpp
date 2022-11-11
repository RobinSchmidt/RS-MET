// construction/destruction:

template<class T>
rsModularInteger<T>::rsModularInteger(rsUint64 initialValue, rsUint64 modulusToUse)
{
  modulus = modulusToUse;
  value   = initialValue;
  rsAssert( value >= T(0) && value < modulus );
  // Although we received unsigned integers for the parameters, value < 0 could happen when there's
  // wraparound in the assignments...right?
}

template<class T>
rsModularInteger<T>::rsModularInteger(const rsModularInteger<T>& other)
{
  modulus = other.modulus;
  value   = other.value;
}

// operators:

template<class T>
rsModularInteger<T> rsModularInteger<T>::operator-() const
{
  if( value == 0 )
    return *this;
  else
    return rsModularInteger<T>(modulus-value, modulus);
}

template<class T>
bool rsModularInteger<T>::operator==(const rsModularInteger<T>& other) const
{
  return value == other.value && modulus == other.modulus;
}

template<class T>
bool rsModularInteger<T>::operator!=(const rsModularInteger<T>& other) const
{
  return !(*this == other);
}

template<class T>
rsModularInteger<T> rsModularInteger<T>::operator+(const rsModularInteger<T> &other)
{
  rsAssert( modulus == other.modulus );
  T r = this->value + other.value;
  if( r >= modulus )
    r -= modulus;
  return rsModularInteger<T>(r, modulus);
}

template<class T>
rsModularInteger<T> rsModularInteger<T>::operator-(const rsModularInteger<T> &other)
{
  rsAssert( modulus == other.modulus );
  T r;
  if( other.value > this->value )
    r = modulus + this->value - other.value;
  else
    r = this->value - other.value;
  return rsModularInteger<T>(r, modulus);
}

template<class T>
rsModularInteger<T> rsModularInteger<T>::operator*(const rsModularInteger<T> &other)
{
  rsAssert( modulus == other.modulus );
  T r = (this->value * other.value) % modulus;
  return rsModularInteger<T>(r, modulus);
}

template<class T>
rsModularInteger<T> rsModularInteger<T>::operator/(const rsModularInteger<T> &other)
{
  rsWarning("Tests needed for: rsModularInteger<T>::operator/");
  rsAssert( modulus == other.modulus );
  return *this * rsModularInteger<T>(rsModularInverse(other.value, modulus), modulus);
  //return *this * rsModularInverse(other.value, modulus); // old, doesn't compile
}
// what, if there is no modular inverse (i think, it exists only if the value is coprime with the 
// modulus - verify)

template<class T>
rsModularInteger<T>& rsModularInteger<T>::operator+=(const rsModularInteger<T> &other)
{
  *this = *this + other;
  return *this;
}
template<class T>
rsModularInteger<T>& rsModularInteger<T>::operator-=(const rsModularInteger<T> &other)
{
  *this = *this - other;
  return *this;
}
template<class T>
rsModularInteger<T>& rsModularInteger<T>::operator*=(const rsModularInteger<T> &other)
{
  *this = *this * other;
  return *this;
}
template<class T>
rsModularInteger<T>& rsModularInteger<T>::operator/=(const rsModularInteger<T> &other)
{
  rsWarning("Tests needed for: rsModularInteger<T>::operator/=");
  *this = *this / other;
  return *this;
}
template<class T>
rsModularInteger<T>& rsModularInteger<T>::operator++()
{
  *this = *this + rsModularInteger<T>(1, modulus);
  return *this;
}
template<class T>
rsModularInteger<T>& rsModularInteger<T>::operator--()
{
  *this = *this - rsModularInteger<T>(1, modulus);
  return *this;
}

/*

Ideas:
-templatize on the integer type to use - allow for arbitrary size integers (done?)
-make it work also for negative values -> verify, if the % operator works correctly, when the 
 left operand is negative - if not, use a custom function instead, see noiseReverseMode() 
 experiment - it has such a function - maybe drag it out into the RAPT library as rsModulo
-how about negative moduli?
-currently, the arithmetic operations make sense only when the two operands have the same modulus
 -generalize this to a sort of "multi-modular" or "mixed-modular" arithmetic
 -the modulus of the result should be the lowest common multiple of the moduli of the operands 
 -i think, it could make sense because in modular arithmetic, we can either take the remainder after 
  each operation or we can just calculate everything in the integers and take the remainder of the 
  final result - but we could (i think) also do all calculations with a modulus that is a multiple 
  of the original modulus and finally take the remainder with respect to the original modulus and 
  still get the same result (verify this) - so if, in the middle of the computations, we encounter 
  different moduli we could choose (the smallest) one, which is "compatible" with both moduli in 
  this sense...would that make sense?
 -what algebraic structure do we get with so defined multimodular integers? is it still a ring?
 -maybe for defining the equality comparison between such multimodular integers, one should compare
  the remainders modulo the gcd of both moduli? ..but would such a definition actually satisfy the 
  constraints for an equivalence relation? ...i think, it breaks transitivity...
-does the notion of a modular rational number make any sense? i.e. numerator and/or denominator are 
 modular integers?
-What about Galois fields? Maybe we can have a class rsGaloisField where the user can set the base
 and exponent. For the special case of an exponent of 1, it would reduce to modular arithmetic in
 modulus p. For higher exponents, we can't just use modular arithmetic anymore. Instead, a more 
 elaborate implementation is necessary, see:
 https://www.youtube.com/watch?v=4BfCmZgOKP8


ToDo: 
-plot lcm(x,y) / (x*y) ...this should be some sort of measure, how small the lcm of of x and y 
 actually is relative to how big it could be at most - a sort of measure how-much-mutually-prime two 
 numbers are as opposed to a simple boolean yes/no, pairs with a smaller value are more simply 
 related, Mutually prime numbers get a value of 1 whereas other numbers get a smaller number. Plot 
 that as a heat-map...maybe also create decimated versions of it. here, i created a plot of the 
 coprimes:
 https://www.facebook.com/MusicEngineer/posts/2296213557119978
 ...maybe run correlation filters with various patterns over it...especially the 
 [[0,1,0],[1,0,1],[0,1,0]] pattern seems to occur a lot

Finding n-th roots of unity:
https://en.wikipedia.org/wiki/Root_of_unity_modulo_n
https://math.stackexchange.com/questions/158344/how-to-find-the-solutions-for-the-n-th-root-of-unity-in-modular-arithmetic

https://people.cs.ksu.edu//~rhowell/calculator/  there's java source-code
https://people.cs.ksu.edu//~rhowell/calculator/how.html

https://doc.sagemath.org/html/en/prep/Quickstarts/Number-Theory.html
https://doc.sagemath.org/html/en/tutorial/tour_numtheory.html


https://mosullivan.sdsu.edu/Teaching/sdsu-sage-tutorial/mathstruct.html


*/