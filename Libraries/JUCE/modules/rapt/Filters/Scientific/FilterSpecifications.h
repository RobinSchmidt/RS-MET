#pragma once

// the classes here are mainly meant for convenient experimentation/prototyping, not for realtime 
// production code

//template <class T>
//struct rsFilterSpecificationBA<T>; // forward dacalation - doesn't work

/** A structure to specify a filter in terms of its zeros, poles, a gain factor and possibly a 
sample rate. If the sample rate is infinity (the default), an analog filter is assumed and the 
zeros and poles are in the s-plane. If it's nonzero, the filter is assumed to be digital and the
zeros and poles are in the z-plane. */

template <class T>
struct rsFilterSpecificationZPK
{
  rsFilterSpecificationZPK() {}
  rsFilterSpecificationZPK(const std::vector<std::complex<T>>& poles,
    const std::vector<std::complex<T>>& zeros, T gain, T sampleRate)
  {
    this->poles = poles;
    this->zeros = zeros;
    this->gain  = gain;
    this->sampleRate = sampleRate;
  }

  //rsFilterSpecificationBA<T> toBA(); // maybe move to FilterCoefficientConverter
  //transferFunctionAt(s_or_z);

  std::vector<std::complex<T>> poles; // rename to p
  std::vector<std::complex<T>> zeros; // rename to z
  std::complex<T> gain = 1;           // rename to k
  T sampleRate = std::numeric_limits<T>::infinity();
};

//=================================================================================================

template <class T>
struct rsFilterSpecificationBA
{
  rsFilterSpecificationBA() {}
  rsFilterSpecificationBA(const std::vector<T>& num, const std::vector<T>& den, T sampleRate)
  {
    this->sampleRate = sampleRate;
    b = num;
    a = den;
  }

  std::vector<std::complex<T>> b; // numerator
  std::vector<std::complex<T>> a; // denominator
  T sampleRate = std::numeric_limits<T>::infinity();
};
// filter specification by polynomial coeffs for numerator (b) and denominator (a)
// H(s) = (b0 + b1*s + ... + bM*s^M) / (a0 + a1*s + ... + aN*s^N)   or
// H(z) = (b0 + b1/z + ... + bM/z^M) / (a0 + a1/z + ... + aN/z^N)