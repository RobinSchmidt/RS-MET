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


  //-----------------------------------------------------------------------------------------------
  /** \name Evaluation */

  /** Evaluates the function at the given input x. */
  T evaluate(const T& x) { return a * exp(-d*x) * sin(w*x + p); }

  // ToDo:
  // -implement () operator for evaluation

protected:

  T w = 0, a = 0, d = 0, p = 0; 


  template<class U> friend class rsDampedSineSum;
};

//=================================================================================================

/** Represents a function given by a weighted sum of damped sinusoids. Just like polynomials, the 
set of such functions actually forms the algebraic structure of a ring: we can add, subtract and 
multiply functions of that type and get again a function of the saem type. This is realized here
by implementing the +,-,* arithmetic operators for the function objects. It's not recomended to use
this feature for realtime synthesis, though - its far more efficient to just add, subtract and 
multiply the output signals of the operand functions. But for proof of concept, it has been 
implemented. */

template<class T>
class rsDampedSineSum
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void addSine(T w, T a, T d, T p) { sines.push_back(rsDampedSine<T>(w, a, d, p)); }


  //-----------------------------------------------------------------------------------------------
  /** \name Evaluation */

  /** Evaluates the function at the given input x. */
  T evaluate(const T& x) 
  { 
    T y(0);
    for(size_t i = 0; i < sines.size(); i++)
      y += sines[i].evaluate(x);
    return y;
  }

  /** Adds two damped sine-sums. */
  rsDampedSineSum<T> operator+(const rsDampedSineSum<T>& q) const;

  /** Multiplies two damped sine-sums. */
  rsDampedSineSum<T> operator*(const rsDampedSineSum<T>& q) const;

  // ToDo:
  // -implement () operator for evaluation
  // -implement unary and binary minus operator

protected:

  std::vector<rsDampedSine<T>> sines; /**< The damped sinusoidal components of this function. */

};

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
  static const T p2 = 0.5*PI;
  size_t k = 0;
  rsDampedSine<T> L, R, S, D;
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
  // compensates via the phase.
}





#endif