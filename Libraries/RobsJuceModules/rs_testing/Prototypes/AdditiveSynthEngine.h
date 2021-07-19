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

/** Abstract baseclass to define the common interface for all our different implementations of a 
bank of sinusoidal oscillators with sweeping frequency and fading amplitude. */

class rsSineSweeperBank
{

public:


  virtual void setMaxNumOscillators(int newLimit) = 0;

  virtual void setNumOscillators(int newNumber) = 0;

  virtual void processFrame(float* left, float* right) = 0; 

  virtual void reset() = 0;
};


//=================================================================================================

/** Implements a bank of sinusoidal oscillators with cubic envelopes for instantaneous phase 
and for the instantaneous amplitude. It's based on maintaining a phasor and computing the sine 
directly either exactly or by various polynomial approximations delivering different amounts of 
fidelity and eating different amounts of CPU time ...tbc... */

template<class T, int N>
class rsSineSweeperBankDirect : public rsSineSweeperBank
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
class rsSineSweeperBankIterative : public rsSineSweeperBank
{

public:

  void setMaxNumOscillators(int newLimit) override;
  void setNumOscillators(int newNumber) override;
  void processFrame(float* left, float* right) override;
  void reset() override;


protected:



  void setMaxNumGroups(int newNumber) { simdGroups.resize(newNumber); }

  using SimdGroup = rsSineSweepIterator<rsSimdVector<T, N>>;

  std::vector<SimdGroup> simdGroups;

};

//=================================================================================================

/** Data structure to represent a patch for a particular key (midi note) for the additive 
synthesis engine. */
struct rsAdditiveKeyPatch
{

  /** Represents the parameters in one of two units, depending on the context. One set of units
  is for presenting them to the user and the other is for internal use in the algorithm. */
  struct SineParams
  {
    float freq  = 0.f;  // in Hz (user), as omega = 2*pi*fs/fs (algo)
    float gain  = 0.f;  // as raw factor (user, direct algos), log of factor (iterative algos)
    float phase = 0.f;  // in degrees (user), radians (algo), some algos may ignore it
  };

  struct Breakpoint
  {
    float time = 0.f;   // in seconds (user), samples (algo)
    std::vector<SineParams> params;
  };


  void addBreakpoint(Breakpoint bp) { breakpoints.push_back(bp); }

  int getNumBreakpoints() const { return (int) breakpoints.size(); }

  const Breakpoint* getBreakpoint(int i) const { return &breakpoints[i]; }

  void clear() { /*startPhases.clear();*/ breakpoints.clear(); }


protected:

  //std::vector<float>      startPhases;
  std::vector<Breakpoint> breakpoints;

  //int key = 0;

};

//=================================================================================================

/** Single voice for the additive synthesis engine based on oscillator banks using SIMD 
processing. It the plural "banks" because different implementations are available with different 
tradeoffs with respect to accuracy and efficiency. */

template<int N>  // N: size of the SIMD vectors
class rsAdditiveSynthVoice
{

public:




  //-----------------------------------------------------------------------------------------------
  // \name Setup


  void setPatch(rsAdditiveKeyPatch* newPatch) { patch = newPatch; }
  // patch must be in the format where the units are for the algorithm

  void setSweeperBank(rsSineSweeperBank* newBank) { sweeperBank = newBank; }


  //-----------------------------------------------------------------------------------------------
  // \name Processing

  void startPlaying();

  void goToBreakpoint(int index, bool reInitAmpAndPhase);

  void processFrame(float* left, float* right);

  void reset();

protected:

  void initSweepers(const rsAdditiveKeyPatch::Breakpoint* bpStart, 
    const rsAdditiveKeyPatch::Breakpoint* bpEnd, bool reInitAmpAndPhase);


  rsSineSweeperBank* sweeperBank = nullptr;

  // The different implementations/algorithms that user can choose from:
  //rsSineSweeperBankDirect<   float,  N> oscsDirectF;
  //rsSineSweeperBankIterative<float,  N> oscsIterF;
  //rsSineSweeperBankIterative<double, N> oscsIterD;
  // Maybe we should make a baseclass rsSineSweeperBank and only maintain a baseclass pointer here.
  // Th common interface should have functions like processFrame, processBlock for both double and
  // float

  rsAdditiveKeyPatch* patch = nullptr;

  int nextBreakpointIndex     = 0;
  int samplesToNextBreakpoint = 0;
  bool playing = false;

  bool alwaysReInit = false; // for test - it's probably not desirable

  //int reInitInterval = 1024;
  /**< Interval for periodic re-initialization of the iterative procedures to avoid indefinite 
  accumulation of roundoff errors. */
  // should be managed by the patch - the patch determines when we re-init

  //double sampleRate = 44100;

};

