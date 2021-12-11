template<class T>
T rsDampedSine<T>::getIntegral(const T& t0, const T& t1) const
{

  // Handle special case where t1 == inf:
  if(t1 == RS_INF(T))
  {
    T tmp = getCompleteIntegral();
    if(t0 != T(0))
      tmp -= getIntegral(T(0), t0);
    return tmp;
  }
  // needs test

  // Handle general case where t0,t1 are both finite numbers:
  T c0 = cos(t0*w + p);
  T s0 = sin(t0*w + p);
  T c1 = cos(t1*w + p);
  T s1 = sin(t1*w + p);
  T e0 = exp(d*t0);
  T e1 = exp(d*t1);
  T s  = d*d + w*w;
  return a*((w*c0 + d*s0)/(s*e0) - (w*c1 + d*s1) / (s*e1));

  // Formula was found with sage via:
  //   var("t a d w p t0 t1")
  //   f(t) = a * exp(-d*t) * sin(w*t + p)
  //   integral(f, t, t0, t1)
  // The result was then simplified by hand.
}
// ToDo:
// -assert t0 is finite, t1 is either finite or +inf but not -inf
// -use rsSinCos
// -What happens if t1 = inf? Will the formula give the correct result? If not, make it work. I
//  guess, e1 will become infinite such that the 2nd term will vanish. maybe fall back to using
//  getCompleteIntegral - but then subtract the integral from 0 to t0, if t0 is nonzero
// -implement the function for rsDampedSineSum - this should just sum of the results for the 
//  partials
// -compute energy integral as sum of two components which result from multiplying the signal
//  with itself
// -maybe implement an energy formula that doesn't require to explicitly compute the product
//  because this requires a memory allocation which we may wnat to avoid in certain contexts
// -Get a formula for the center of gravity. I think, that just requires to multiply in a t, 
//  such that f(t) = t * a * exp(-d*t) * sin(w*t + p)

template<class T>
T rsDampedSine<T>::getCompleteIntegral() const
{
  return (a*w*cos(p) + a*d*sin(p))/(d*d + w*w);

  // Formula was found with sage via:
  //   var("t a d w p")
  //   f(t) = a * exp(-d*t) * sin(w*t + p)
  //   assume(d > 0)
  //   I = integral(f, t, 0, oo)
  //   I.simplify_full()
}


template<class T>
T rsDampedSine<T>::getEnvelopeIntegral(const T& t0, const T& t1) const
{
  return a*(exp(-d*t0)/d - exp(-d*t1)/d); // ToDo: use only 1 division

  // Formula was found with sage via:
  //   var("t a d w p t0 t1")
  //   f(t) = a * exp(-d*t) 
  //   integral(f, t, t0, t1)
}

template<class T>
T rsDampedSine<T>::getCompleteEnvelopeIntegral() const
{
  return a / d;

  // Formula was found with sage via:
  //   var("t a d w p")
  //   f(t) = a * exp(-d*t)
  //   assume(d > 0)
  //   integral(f, t, 0, oo)
}

template<class T>
T rsDampedSine<T>::getEnergyIntegral(const T& t0, const T& t1) const
{
  T d2   = d*d;
  T w2   = w*w;
  T d2w2 = d2+w2;
  T D    = 4*d*d2w2;
  T s, c; 
  rsSinCos(2*(p+w*t1), &s, &c);
  T F1 = (d*w*s-d2*c+d2w2) / (D*exp(2*d*t1)); // F(t1), up to scaling
  rsSinCos(2*(p+w*t0), &s, &c);
  T F0 = (d*w*s-d2*c+d2w2) / (D*exp(2*d*t0)); // F(t0), up to scaling
  return -a*a * (F1 - F0);
}
/*
var("t a d w p t0 t1")
f(t) = a * exp(-d*t) * sin(w*t + p)
assume(t1-t0 > 0)
I = integral(f^2, t, t0, t1)
I.simplify_full()
# give names to common subexpressions:
var('c2p s2p e0 e1 cs2')
I = I.subs(cos(2*p)==c2p)
I = I.subs(sin(2*p)==s2p)
I = I.subs(e^(2*d*t0)==e0)
I = I.subs(e^(2*d*t1)==e1)
I = I.subs(c2p^2 + s2p^2==cs2)
I

...this gives quite a messy expression...
perhaps it makes sense to set t0=0
*/


template<class T>
T rsDampedSine<T>::getCenterOfMass() const
{
  T d2 = d*d;
  T w2 = w*w;
  T C  = (2*a*d*w*cos(p) + a*d2*sin(p) - a*w2*sin(p))/(d2*d2 + 2*d2*w2 + w2*w2);
  return C / getCompleteIntegral();  // verify!
  // todo: 
  // -precompute cp, sp as cos(p), sin(p) using rsSinCos
  // -factor out a, then (d2-w2) from sin terms
  // -simplify denominator using binomial formula: its just (d2+w2)^2
  // -can the numerator also be simplified?
}
// Sage code for center of mass:
//   var("t a d w p t0 t1")
//   f(t) = t * a * exp(-d*t) * sin(w*t + p)
//   assume(t1 - t0 > 0)
//   cog = integral(f, t, t0, t1)
//   cog.full_simplify()
// gives a quite monstrous result, even with full simplification. ToDo: try to get sage to 
// simplify it further by manually defining some intermediate variables for subexpressions.
// See here: https://ask.sagemath.org/question/8403/full-simplify-sage-vs-mathematica/
// Oh - but wait - we typically don't need the cog of a partial signal! We can replace t0,t1
// by 0 and inf. Oh yes - that's much better:  
//   var("t a d w p")
//   f(t) = t * a * exp(-d*t) * sin(w*t + p)
//   assume(d > 0)
//   cog = integral(f, t, 0, oo)
//   cog.full_simplify()
// gives:
//   (2*a*d*w*cos(p) + a*d^2*sin(p) - a*w^2*sin(p))/(d^4 + 2*d^2*w^2 + w^4)
//
// ToDo: nomralize the value by dividing by the complete integral, see:
//   https://en.wikipedia.org/wiki/Centroid#By_integral_formula
//   https://en.wikipedia.org/wiki/Center_of_mass
// wait - isn't this just the mean value?
//   https://en.wikipedia.org/wiki/Expected_value#Absolutely_continuous_case
// hmm...there the formula is not normalized by the complete integral - but that may be due to the 
// fact that this integral is assumed to be 1 anyway, because f is assumed to be a probability 
// density

// ToDo: Write functions that compute the same quantities for the envelope. That probably more 
// useful. Just remove the sin(...) factor from the formula that is passed to sage



template<class T>
T rsDampedSine<T>::getEnvelopeCenterOfMass() const
{
  return a / (d*d);
  // var("t a d")
  // f(t) = t * a * exp(-d*t)
  // assume(d > 0)
  // cog = integral(f, t, 0, oo)
  // cog.full_simplify()
}
// we should probably be normalized by getCompleteEnvelopeIntegral() ...then the result would just
// be(a/d^2) / (a/d) = 1/d  ...the time constant - maybe rename to getTimeConstant

//=================================================================================================

template<class T>
T rsDampedSineSum<T>::getIntegral(const T& t0, const T& t1) const
{
  T r(0);
  for(const auto & s : sines)
    r += s.getIntegral(t0, t1);
  return r;
  // maybe use std::accumulate
}

template<class T>
T rsDampedSineSum<T>::getCompleteIntegral() const
{
  T r(0);
  for(const auto & s : sines)
    r += s.getCompleteIntegral();
  return r;
  // maybe use std::accumulate
}

template<class T>
T rsDampedSineSum<T>::getEnergyIntegral(const T& t0, const T& t1) const
{
  T r(0);
  for(const auto & s : sines)
    r += s.getEnergyIntegral(t0, t1);
  return r;
  // maybe use std::accumulate
}



template<class T>
T rsDampedSineSum<T>::getCenterOfMass() const
{
  T weightedSumOfCenters(0);
  T sumOfWeights(0);
  for(const auto& s : sines)
  {
    T center = s.getCenterOfMass();
    T weight = s.getCompleteIntegral();
    weightedSumOfCenters += weight * center;
    sumOfWeights += weight;
  }
  return weightedSumOfCenters / sumOfWeights;
  // ToDo: verify this formula!



  /*
  T r(0);
  for(const auto & s : sines)
    r += s.getCenterOfMass();
  return r / getNumSines();
  */
  // ToDo: verify that the arithmetic mean is actually correct...no, i think, it's not - it should
  // probably be a weighted mean where the weights are given by the respective total mass of the
  // partial..
}


/*

ToDo:
-rsDampedSine: implement functions to compute: total energy, energy up to time t, 
 center of gravity. All of these should be available for the whole signal and for the envelope. 
 Total energy and center of gravity of envelope could be used as sorting criteria in place of
 amplitude and decay time. In a way, total energy and center of gravity can be seen as a decoupling
 of amplitude and decay. If we take the values of the whole signal rather than the envelope, it 
 will also "decouple" the values from freq and phase - but i think, that's not really desirable
 for sorting. Formulas for some of that can be found in MathExperiments.cpp starting a around line
 3000 starting with function rsDampedSineIntegral. Verify the formulas there with sage, drag them
 into class rsDampedSine and add the sage code to the documentation.
-Try to find similar formulas also for a sum of damped sines. Maybe the squaring that occurs should
 first be expressed by turning the product into a sum.
-implement functions for derivative and integral, if possible, i.e. if the result is again a sum of 
 damped sines (i think so, but i'm not sure). The definite integral is just a number, so it can be
 computed anyway. I think, the integral look like a constant minus a decay (RC loading curve) with a 
 sine-wiggle on top...but no - it's monotonic...at least in the case of the energy integral, which
 is a the integral of a damped sine multiplied by itself...we'll see...
-implement efficient computation of an array of values by means of recursion
 -maybe this function should accumulate into the array, i.e. use +=, because in practice, we may 
  want to work with sums of damped sines

rsDampedSineSum:
-Maybe keep the array of components sorted in some way - maybe sort according to frequency first 
 (ascending), then according to total energy of the envelope (descending), then according to 
 amplitude (descending), then according to phase (ascending).
-When the array is known to be sorted, we can also meaningfully imlement the == operator.
...but maybe that is not useful - we'll see...


Math for the multiplication:
  
  f(x) = a1 * exp(-d1*x) * sin(w1*x + p1)  *  a2 * exp(-d2*x) * sin(w2*x + p2)
       = a1*a2 * exp(-(d1+d2)*x) * sin(w1*x + p1) * sin(w2*x + p2)
 
Then, for the product of sines, we use the general result:

  sin(a*x+b)*sin(c*x+d) = ( cos((a-c)*x+(b-d)) - cos((a+c)*x+(b+d)) ) / 2

which can be verified via wolfram alpha. For the resulting frequency and phase, there is some 
freedom: we can always choose a positive or negative frequency and then compensate by adding or 
subtracting pi from the phase. When we flip the sign of the frequency, we need to compensate by a 
phase-shift of 180°. ToDo: maybe we should always enforce positive frequencies in the result?

Questions:
-Can we somehow do division (with remainder)? Or factorization? Maybe first we need to search for 
 partials with equal amplitudes and decay-rates because products alway create such pairs of 
 partials. If P=A*B and a partial in P occurrd from a product partials from A and B, then the 
 product's amplitude must be the product of the amplitudes of the factor partials from A and B. 
-What if we form a product of damped sine filters? Will this still be a filter, i.e. an LTI 
 system? ...probably not...nope...i don't think so..but maybe some sort of nonlinear "filter"
-If we want to create a given target spectrum, say, 12 harmonics at 100,200,...,1200 - can we do
 it by multiplying a sum of 2 and a sum of 3 sines? This will produce 2*2*3=12 partials - but how 
 do we choose the freqs of the partials to achieve a given target product spectrum?
-What about a generalization including a frequency variation term v, like 
   f(t) = a * exp(-d*t) * sin(v*t^2 + w*t + p)
 ...maybe make a class rsDampedSineSweep and rsDampedSineSweepSum. Is a product of such signals
 again of the same type? Would be nice, but even if not, it may be useful for synthesis. The 
 algebraic properties are useful for predicting output spectra analytically but are not necessary
 to build an interesting synthesis scheme. Hmm..it's actually looking good: 
   https://www.wolframalpha.com/input/?i=sin%28a+x%5E2%29+sin%28b+x%5E2%29
   sin(a x^2) sin(b x^2) = (1/2) (cos(a x^2 - b x^2) - cos(a x^2 + b x^2))
 and it generalizes further to arbitray powers:
   sin(a x^n) sin(b x^n) = (1/2) (cos(a x^n - b x^n) - cos(a x^n + b x^n))
-What about exponential frequency variation?


*/