#pragma once

// the classes here are mainly meant for convenient experimentation/prototyping, not for realtime 
// production code

template <class T>
struct rsFilterSpecificationBA;

/** A structure to specify a filter in terms of its zeros, poles, a gain factor and possibly a 
sample rate. If the sample rate is infinity (the default), an analog filter is assumed and the 
zeros and poles are in the s-plane. If it's nonzero, the filter is assumed to be digital and the
zeros and poles are in the z-plane. In the analog case, the transfer function is represented as:

            (s-q1) * (s-q2) * ... * (s-qM)   
H(s) = k * --------------------------------
            (s-p1) * (s-p2) * ... * (s-pN)

and in the digital case as:

            (1-q1/z) * (1-q2/z) * ... * (1-qM/z)
H(z) = k * --------------------------------------
            (1-p1/z) * (1-p2/z) * ... * (1-pN/z)

where the qi and pi are the (s- or z-plane) zeros and poles respectively and k is an overall gain 
factor. This representation of the transfer function may be called the pole/zero form and is useful
to analyze stability, mininum-phase properties, ringing time, ... */

template <class T>
struct rsFilterSpecificationZPK
{
  rsFilterSpecificationZPK() {}
  rsFilterSpecificationZPK(
    const std::vector<std::complex<T>>& zeros, 
    const std::vector<std::complex<T>>& poles,
    std::complex<T> gain, T sampleRate);


  bool isDigital() { return sampleRate != std::numeric_limits<T>::infinity(); }
  std::complex<T> transferFunctionAt(std::complex<T> s_or_z);
  rsFilterSpecificationBA<T> toBA(); // maybe move to FilterCoefficientConverter

  void sortPolesAndZeros();

  bool equals(const rsFilterSpecificationZPK& other, T tolerance = T(0));

  std::vector<std::complex<T>> z; // zeros
  std::vector<std::complex<T>> p; // poles
  std::complex<T> k = 1;          // gain
  T sampleRate = std::numeric_limits<T>::infinity();
};

//=================================================================================================

/** A structure to specify a filter in terms of the polynomial coefficients of the numerator and 
denominator. In the analog case, the transfer function is:

        B0 + B1*s + B2*s^2 + ... + BM*s^M
H(s) = -----------------------------------
        A0 + A1*s + A2*s^2 + ... + AN*s^N

and in the digital case:

        b0 + b1/z + b2/z^2 + ... + bM/z^M
H(z) = -----------------------------------
        a0 + a1/z + a2/z^2 + ... + aN/z^N

where we typically set AN = 1 in the analog case and a0 = 1 in the digital case (this doesn't 
change the transfer function - it amounts to divide numerator and denominator by AN or a0). 
This representation of the transfer function may be called the coefficient form or direct form and
is useful for obtaining transfer-functions of weighted sums of filters, computing group-delay (i 
guess - because it involves taking a derivative), ... */

template <class T>
struct rsFilterSpecificationBA
{
  rsFilterSpecificationBA() {}
  rsFilterSpecificationBA(const std::vector<std::complex<T>>& num, 
    const std::vector<std::complex<T>>& den, T sampleRate);

  bool isDigital() { return sampleRate != std::numeric_limits<T>::infinity(); }
  std::complex<T> transferFunctionAt(std::complex<T> s_or_z);
  rsFilterSpecificationZPK<T> toZPK(); // maybe move to FilterCoefficientConverter

  /** Normalizes the coefficients such that a[0] in the digital and a[last] in the analog case is
  unity. */
  void normalizeDenominator();


  std::vector<std::complex<T>> b; // numerator
  std::vector<std::complex<T>> a; // denominator
  T sampleRate = std::numeric_limits<T>::infinity();
};

// maybe factor out a baseclass rsFilterSpecification that has only a sampleRate member and a 
// function isDigital()...maybe also a virtual member transferFunctionAt that is overriden
// by the two subclasses, maybe protect data members and make the classes mutual friends
// maybe let the FilterPlotter keep an array of baseclass pointers