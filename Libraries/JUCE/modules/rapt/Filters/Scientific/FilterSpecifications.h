#pragma once

// the classes here are mainly meant for convenient experimentation/prototyping, not for realtime 
// production code

template <class T>
struct rsFilterSpecificationBA;

/** A structure to specify a filter in terms of its zeros, poles, a gain factor and possibly a 
sample rate. If the sample rate is infinity (the default), an analog filter is assumed and the 
zeros and poles are in the s-plane. If it's nonzero, the filter is assumed to be digital and the
zeros and poles are in the z-plane. */

template <class T>
struct rsFilterSpecificationZPK
{
  rsFilterSpecificationZPK() {}
  rsFilterSpecificationZPK(
    const std::vector<std::complex<T>>& zeros, 
    const std::vector<std::complex<T>>& poles,
    std::complex<T> gain, T sampleRate);

  std::complex<T> transferFunctionAt(std::complex<T> s_or_z);
  rsFilterSpecificationBA<T> toBA(); // maybe move to FilterCoefficientConverter

  std::vector<std::complex<T>> z; // zeros
  std::vector<std::complex<T>> p; // poles
  std::complex<T> k = 1;          // gain
  T sampleRate = std::numeric_limits<T>::infinity();
};

//=================================================================================================

template <class T>
struct rsFilterSpecificationBA
{
  rsFilterSpecificationBA() {}
  rsFilterSpecificationBA(const std::vector<std::complex<T>>& num, 
    const std::vector<std::complex<T>>& den, T sampleRate);

  std::complex<T> transferFunctionAt(std::complex<T> s_or_z);
  rsFilterSpecificationZPK<T> toZPK(); // maybe move to FilterCoefficientConverter

  /** Normalizes the a[0] coefficient to unity by dividing all a- and b-coeffs by a0. That 
  doesn't change the overall transfer function. */
  void normalizeA0();


  std::vector<std::complex<T>> b; // numerator
  std::vector<std::complex<T>> a; // denominator
  T sampleRate = std::numeric_limits<T>::infinity();
};
// filter specification by polynomial coeffs for numerator (b) and denominator (a)
// H(s) = (b0 + b1*s + ... + bM*s^M) / (a0 + a1*s + ... + aN*s^N)   or
// H(z) = (b0 + b1/z + ... + bM/z^M) / (a0 + a1/z + ... + aN/z^N)

// maybe factor out a baseclass rsFilterSpecification that has only a sampleRate member and a 
// function isDigital()...maybe also a virtual member transferFunctionAt that is overriden
// by the two subclasses, maybe protect data members and make the classes mutual friends
// maybe let the FilterPlotter keep an array of baseclass pointers