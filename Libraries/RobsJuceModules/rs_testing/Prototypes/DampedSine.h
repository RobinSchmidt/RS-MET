#ifndef RAPT_DAMPEDSINE_H
#define RAPT_DAMPEDSINE_H

/** Class to represent a damped sinusoid of the form: f(x) = a * exp(-d*x) * sin(w*x + p) as a 
function object. ...tbc... */

template<class T>
class rsDampedSine
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  rsDampedSine() {}

  rsDampedSine(T radianFreq, T amplitude, T decay, T phase = T(0))
  { setup(radianFreq, amplitude, decay, phase); }


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets up the parameters for frequency, amplitude, decay and phase for the function:
    f(x) = a * exp(-d*x) * sin(w*x + p)
  The frequency variable w stands for "omega" and is actually the radian frequency, i.e. 
  w = 2*pi*f. The amplitude a is a raw multiplier. The decay d is 1/tau where tau is the time 
  constant, i.e. the time it takes to decay to 1/e. The start phase p is in radians. */
  void setup(T w, T a, T d, T p = T(0)) { this->w = w; this->a = a; this->d = d; this->p = p; }


  void canonicalize();


  //-----------------------------------------------------------------------------------------------
  /** \name Evaluation */

  // ToDo:
  // -consistently use t instead of x as input variable
  // -implement () operator for evaluation
  // -implement getEnergyIntegral, getTotalEnergy, getEnergyCentroid - then also for a sum of sines
  // -the computed values for the energy-stuff, e.g. getEnergyIntegral, should be the same as the 
  //  corresponding getIntegral, etc. values when applied to f^2 = f*f -> this can be tested
  // same functions applied to
  //  non-squared signals


  /** Evaluates the function at the given input x. */
  T evaluate(const T& x) const { return a * exp(-d*x) * sin(w*x + p); }

  /** Computes the definite integral of the signal between time instants t0 and t1. */
  T getIntegral(const T& t0, const T& t1) const;

  /** Computes the definite integral of the signal between time instants zero and infinity. */
  T getCompleteIntegral() const;

  /** Computes the definite integral of the signal's envelope between time instants t0 and t1. */
  T getEnvelopeIntegral(const T& t0, const T& t1) const;

  /** Computes the definite integral of the signal's envelope between time instants zero and 
  infinity. */
  T getCompleteEnvelopeIntegral() const;

  T getEnergyIntegral(const T& t0, const T& t1) const;



  /** Computes the center of mass of the signal, i.e. the time instant where most of the "mass" is
  concentrated, where the function value itself is taken to be a "density". Negative function 
  values contribute "negative density", so the value can actually come out negative and is 
  typically close to zero. The notion of what we compute here is typically not really useful in its
  own right when  applied to a single damped sinuosid. It becomes useful when it's applied to a 
  squared sum of damped sinusoids ...tbc...   */
  T getCenterOfMass() const;

  /** Computes the center of mass of the envelope, i.e. the time instant where most of the "mass" 
  is concentrated where the enevlope value is taken to be a density. */
  T getEnvelopeCenterOfMass() const;
  // maybe rename to get(Envelope)Centroid ...but maybe not - the center of mass takes into account 
  // the density (here taken as the function value) whereas the centroid does not
  // maybe rename to getValueCentroid, getEnergyCentroid, getEnvelopeCentroid




  /** Compares sinusoids for one being less than another. The main criterion here is frequency. 
  When sorting an array of partials, we want the lower frequencies to appear first. ...tbc... */
  bool operator<(const rsDampedSine<T>& q) const;


protected:

  T w = 0, a = 0, d = 0, p = 0; 


  template<class U> friend class rsDampedSineSum;
};

template<class T>
void rsDampedSine<T>::canonicalize()
{
  if(w < 0)
  {
    w = -w;
    a = -a;
  }
  if(a < 0)
  {
    a = -a;
    p = RAPT::rsWrapToInterval(p + PI, -PI, +PI);
  }
}
// todo: test, if this really produces the same sine

template<class T>
bool rsDampedSine<T>::operator<(const rsDampedSine<T>& r) const
{
  if(w < r.w) return true;  // lower freqs appear first
  if(r.w < w) return false;
  if(a > r.a) return true;  // louder partials appear first
  if(r.a > a) return false;
  if(d < r.d) return true;  // longer partials appear first
  if(r.d < d) return false;
  if(p < r.p) return true; 
  return false;
  // ToDo: maybe make the 2nd check based on the total energy, not amplitude alone. I think, it's
  // E = a/d or something -> look it up or derive. Then, if two partials have the same energy, use
  // amplitude alone as next comparison - this puts those partials in front, whose early energy is
  // greater (i think))
}


//=================================================================================================

/** Represents a function given by a weighted sum of damped sinusoids. Just like polynomials, the 
set of such functions actually forms the algebraic structure of a ring: we can add, subtract and 
multiply functions of that type and get again a function of the same type. This is realized here
by implementing the +,-,* arithmetic operators for the function objects. It's not recomended to use
this feature for realtime synthesis, though - its far more efficient to just add, subtract and 
multiply the output signals of the operand functions. But for proof of concept, it has been 
implemented. It can also be useful for analyzing the spectra that would be produced by such 
multiplications. */

template<class T>
class rsDampedSineSum
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void addSine(T w, T a=1, T d=0, T p=0) { sines.push_back(rsDampedSine<T>(w, a, d, p)); }


  void setSine(T w, T a=1, T d=0, T p=0) { clear(); addSine(w, a, d, p); }


  void scaleFrequencies(T scale);
  void scaleAmplitudes(T scale);
  void scaleDecayRates(T scale);
  void shiftPhases(T shift);

  void sortPartials();
  void canonicalize();


  void clear() { sines.clear(); }



  //-----------------------------------------------------------------------------------------------
  /** \name Evaluation */

  /** Evaluates the function at the given input x. */
  T evaluate(const T& x) const
  { 
    T y(0);
    for(size_t i = 0; i < sines.size(); i++)
      y += sines[i].evaluate(x);
    return y;
  }


  /** Computes the definite integral of the signal between time instants t0 and t1. */
  T getIntegral(const T& t0, const T& t1) const;

  /** Computes the definite integral of the signal between time instants zero and infinity. */
  T getCompleteIntegral() const;


  T getEnvelopeIntegral(const T& t0, const T& t1) const;

  T getCenterOfMass() const;

  T getEnvelopeCenterOfMass() const;


  int getNumSines() const { return (int) sines.size(); }

  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Adds two damped sine-sums. */
  rsDampedSineSum<T> operator+(const rsDampedSineSum<T>& q) const;

  /** Subtracts two damped sine-sums. */
  rsDampedSineSum<T> operator-(const rsDampedSineSum<T>& q) const;

  /** Multiplies two damped sine-sums. */
  rsDampedSineSum<T> operator*(const rsDampedSineSum<T>& q) const;

  /** Multiplies a scalar and a damped sine-sum. */
  rsDampedSineSum<T> operator*(const T& s) const;


  // unary minus:
  rsDampedSineSum<T> operator-() const;


  // we can implement this in 3 ways: negate amp, negate freq, shift phase...or wait - is this 
  // correct? ..i think freq and phase are linked somehow. anyway, we should probably implement it
  // by negating the amp. that's most natural



  // ToDo:
  // -implement () operator for evaluation
  // -implement unary and binary minus operator

protected:

  std::vector<rsDampedSine<T>> sines; /**< The damped sinusoidal components of this function. */

};

template<class T>
void rsDampedSineSum<T>::scaleAmplitudes(T scale)
{
  for(size_t i = 0; i < sines.size(); i++)
    sines[i].a *= scale;
}
template<class T>
void rsDampedSineSum<T>::scaleFrequencies(T scale)
{
  for(size_t i = 0; i < sines.size(); i++)
    sines[i].w *= scale;
}
template<class T>
void rsDampedSineSum<T>::scaleDecayRates(T scale)
{
  for(size_t i = 0; i < sines.size(); i++)
    sines[i].d *= scale;
}
template<class T>
void rsDampedSineSum<T>::shiftPhases(T shift)
{
  for(size_t i = 0; i < sines.size(); i++)
    sines[i].p += shift;
}
// ToDo: use functions rsScale, rsAdd for vectors (if they don't exist, implement them). Then move
// the code in the class (it will become a short one-liner then)

template<class T>
void rsDampedSineSum<T>::sortPartials()
{
  std::sort(sines.begin(), sines.end());
  //RAPT::rsHeapSort(&sines[0], (int)sines.size());
}

template<class T>
void rsDampedSineSum<T>::canonicalize()
{
  for(auto & s : sines)
    s.canonicalize();
  sortPartials();
}

template<class T>
rsDampedSineSum<T> rsDampedSineSum<T>::operator-() const 
{ 
  rsDampedSineSum<T> r; 
  r.sines = sines;
  r.scaleAmplitudes(T(-1));
  return r;
}

template<class T>
rsDampedSineSum<T> rsDampedSineSum<T>::operator-(const rsDampedSineSum<T>& q) const
{
  return *this + (-q);
}
// can perhaps be optimized with respect to creating temporaries

template<class T>
rsDampedSineSum<T> rsDampedSineSum<T>::operator+(const rsDampedSineSum<T>& q) const
{
  rsDampedSineSum<T> r; 
  r.sines.reserve(this->sines.size() + q.sines.size());
  rsAppend(r.sines, this->sines);
  rsAppend(r.sines, q.sines);
  return r;
}

template<class T>
rsDampedSineSum<T> rsDampedSineSum<T>::operator*(const rsDampedSineSum<T>& q) const
{
  rsDampedSineSum<T> r; 
  r.sines.resize(2 * this->sines.size() * q.sines.size());
  static const T p2 = 0.5*PI;          // offset for sum/difference phases
  rsDampedSine<T> L, R, S, D;          // temporaries
  size_t k = 0;                        // output array index
  for(size_t i = 0; i < sines.size(); i++) {
    for(size_t j = 0; j < q.sines.size(); j++) {
      L   =   sines[i];                // current left factor
      R   = q.sines[j];                // current right factor
      S.w = L.w + R.w;                 // sum frequency
      S.p = L.p + R.p - p2;            // sum phase
      D.w = L.w - R.w;                 // difference frequency
      D.p = L.p - R.p + p2;            // difference phase
      S.a = D.a = T(0.5) * L.a * R.a;  // amplitudes are both half of the product
      S.d = D.d = L.d + R.d;           // decay rates are both given by the sum
      r.sines[k]   = S;                // write sum sinusoid to output
      r.sines[k+1] = D;                // write difference sinusoid to output
      k += 2;   }}
  return r;
  // Note: Using D.w = R.w - L.w; D.p = R.p - L.p + p2; works also. It negates the frequency and
  // compensates via the phase. We have some freedom of choice here and the choice i have made 
  // seemed to be the most natural one to me. We could probably also negate the amplitude for 
  // compensation (ToDo: try it).
}

template<class T>
rsDampedSineSum<T> rsDampedSineSum<T>::operator*(const T& s) const
{
  rsDampedSineSum<T> r = *this;
  r.scaleAmplitudes(s);
  return r;
}
// implement * for left argument being a scalar



// ToDo: move implementation to cpp file



#endif