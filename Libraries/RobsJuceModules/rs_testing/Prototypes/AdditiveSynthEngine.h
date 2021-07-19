#pragma once

// just a stub at the moment

/** Under construction. Should provide an API and functionaliyt similar to rsSineSweepIterator, but 
use a phasor and direct computaion (either exact or by using approximations). The amplitude 
envelope should probably be cubic in the raw amplitude domain rather than in log-amplitude domain
here...we'll see...  */

template<class T>
class rsSineSweepOsc
{

public:

protected:

  T pos;  // position in 0..1

};

//=================================================================================================

/** Implements a bank of sinusoidal oscillators with cubic envelopes for instantaneous phase 
and for the instantaneous amplitude. It's based on maintaining a phasor and computing the sine 
directly either exactly or by various polynomial approximations delivering different amounts of 
fidelity and eating different amounts of CPU time ...tbc... */

template<class T, int N>
class rsSineSweeperBankDirect
{

public:

protected:

};

//=================================================================================================

/** Implements a bank of recursively implemented sinusoidal oscillators with cubic envelopes for
instantaneous phase and for the log of instantaneous amplitude based on rsSineSweepIterator. It 
uses SIMD processing (based on rsSimdVector) to combine a bunch of such generators together into 
groups that are processed in parallel. The template parameter T is the underlying scalar type and 
N is the size of the groups, i.e. the size of the simd vectors. ...tbc... */

template<class T, int N>
class rsSineSweeperBankIterative
{

public:


  void setMaxNumGroups(int newNumber) { sweeperGroups.resize(newNumber); }


protected:


  using SweeperGroup = rsSineSweepIterator<rsSimdVector<T, N>>;

  std::vector<SweeperGroup> sweeperGroups;

};

//=================================================================================================

/** An additive synthesis engine based on oscillator banks using SIMD processing. It the plural 
"banks" because different implementations are available with different tradeoffs with respect to
accuracy and efficiency. */

template<int N>  // N: size of the SIMD vectors
class rsAdditiveSynthEngine
{

public:


protected:


  // The different implementations/algorithms that user can choose from:
  rsSineSweeperBankDirect<   float,  N> oscsDirectF;
  rsSineSweeperBankIterative<float,  N> oscsIterF;
  rsSineSweeperBankIterative<double, N> oscsIterD;
  // Maybe we should make a baseclass rsSineSweeperBank and only maintain a baseclass pointer here.
  // Th common interface should have functions like processFrame, processBlock for both double and
  // float



};

