#ifndef RAPT_DAMPEDSINE_H
#define RAPT_DAMPEDSINE_H

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

protected:

  T w = 0, a = 0, d = 0, p = 0; 


  template<class U> friend class rsDampedSineSum;
};

//=================================================================================================

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
  // implement () operator
  // implement operators +,*. Addition should just concatenate the arrays. Multiplication will 
  // produce for each pair of factors a pair of sines and the sum and difference frequencies with 
  // amplitudes given by the products, decay factors given by the sum, etc....

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
  rsDampedSine<T> L, R, S, D; // left and right factor of current pair, sum and difference tone

  for(size_t i = 0; i < sines.size(); i++)
  {
    for(size_t j = 0; j < q.sines.size(); j++)
    {
      L   =   sines[i];                // current left factor
      R   = q.sines[j];                // current right factor

      S.w = L.w + R.w;                 // sum frequency
      D.w = R.w - L.w;                 // difference frequency
      S.p = L.p + R.p - p2;            // sum phase
      D.p = R.p - L.p + p2;            // difference phase
      S.a = D.a = T(0.5) * L.a * R.a;  // amplitudes are both given by the product
      S.d = D.d = L.d + R.d;           // decay rates are both given by the sum

      r.sines[k]   = S;                // write sum sinusoid to output
      r.sines[k+1] = D;                // write difference sinusoid to output
      k += 2;
    }
  }
  // todo: verify formulas theoretically and experimentally, check especially the signs of the 
  // results of the differences ...swap order of arguments for D.w and D.p and compensate by using 
  // -p2 instead of +p2
  // -we seem to have some freedom with respect to the sign of the difference frequency (which can 
  //  then be compensated by the sign of the difference phase). Maybe implement it in such a way 
  //  that when both input freqs are positive, the output freq is also positive..or maybe such that
  //  the output freq is always positive


  //rsAppend(r.sines, this->sines);
  //rsAppend(r.sines, q.sines);
  return r;
}





#endif