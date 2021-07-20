#pragma once

// just a stub at the moment

//=================================================================================================

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

  virtual void init(int index, float t0, float t1, float w0, float w1, float a0, float a1, 
    float p0, float p1, float r0, float r1) = 0;

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
  void init(int index, float t0, float t1, float w0, float w1, float a0, float a1, 
    float p0, float p1, float r0, float r1) override;
  void processFrame(float* left, float* right) override;
  void reset() override;


protected:



  void setMaxNumGroups(int newNumber) { simdGroups.resize(newNumber); }

  using SimdGroup = rsSineSweepIterator<rsSimdVector<T, N>>;

  std::vector<SimdGroup> simdGroups;

};

//=================================================================================================

/** Data structure for 2D arrays. It's a bit like rsMatrix, but much simpler. It  doesn't have all 
the math operations and also doesn't necessarily assume the 2 dimensions to be "rows" and 
"columns".  */

template<class T>
class rsArray2D    // move to RAPT into Data folder
{



  void setShape(int newDim1, int newDim2)
  {
    dim1 = newDim1;
    dim2 = newDim2;
    data.resize(dim1*dim2);
  }
  // todo: optionally retain old data (true by default - rationale: it's safer and optimizations 
  // that sacrifice safety should be explicit)


  /** Read and write access to elements with index pair i,j. */
  T& operator()(const int i, const int j) { return data[flatIndex(i, j)]; }

  /** Converts a pair of indices i,j to a flat array index. */
  int flatIndex(const int i, const int j) const
  {
    rsAssert(i >= 0 && i < dim1, "Invalid index along 1st dimension");
    rsAssert(j >= 0 && j < dim2, "Invalid index along 2nd dimension");
    return dim2*i + j;
  }


protected:

  int dim1 = 0, dim2 = 0;
  std::vector<T> data;

};

//=================================================================================================

/** Data structure to represent a patch for a particular key (midi note) for the additive 
synthesis engine. */
struct rsAdditiveKeyPatch  // turn into class, rename into rsAdditiveKeyPatchEditable
{

  /** Represents the parameters in one of two units, depending on the context. One set of units
  is for presenting them to the user and the other is for internal use in the algorithm. */
  struct SineParams
  {
    SineParams(){}
    SineParams(float freq_, float gain_, float phase_ = 0.f) 
      : freq(freq_), gain(gain_), phase(phase_) { }

    float freq  = 0.f;  // in Hz (user), as omega = 2*pi*fs/fs (algo)
    float gain  = 0.f;  // as raw factor (user, direct algos), log of factor (iterative algos)
    float phase = 0.f;  // in degrees (user), radians (algo), some algos may ignore it
  };

  struct Breakpoint  // maybe mae it a class
  {
    float time = 0.f;   // in seconds (user), samples (algo)

    void addSine(const SineParams& newSine) { params.push_back(newSine); }

    std::vector<SineParams> params;
  };

  void addBreakpoint(Breakpoint bp) { breakpoints.push_back(bp); }

  void clear() { /*startPhases.clear();*/ breakpoints.clear(); }


  void convertUserToAlgoUnits(float sampleRate, bool logAmplitude);



  int getNumBreakpoints() const { return (int) breakpoints.size(); }

  const Breakpoint* getBreakpoint(int i) const { return &breakpoints[i]; }




protected:

  //std::vector<float>      startPhases;
  std::vector<Breakpoint> breakpoints;

  //int key = 0;

};
// make this class internal to rsAdditiveSynthVoice to reduce library surface area

//=================================================================================================

/** This also represents a patch for a particular key like rsAdditiveKeyPatchEditable, but has a
different data layout to facilitate efficient playback. This layout but impedes efficient editing,
so for editing patches, the other class is used and once a patch is finished, it can be converted
into a more "playable" version in one go as preprocessing. */

class rsAdditiveKeyPatchPlayable
{

public:

  /** Sets up this object from the given editable patch. */
  void setupFrom(const rsAdditiveKeyPatch& patch);

protected:

  int simdSize = 0;
  // maybe this needs to be a template parameter because we may need to use 
  // rsSineSweepParameters<rsSimdVector>

  int numBreakpoints = 0;
  int numPartials = 0;
  std::vector<int> timeStamps;


  
  // rsArray2D<rsSineSweepParameters> params;
  // 1st index: breakpoint, 2nd: partial (or maybe simd-group of partials)

};
// make this class internal to rsAdditiveSynthVoice to reduce library surface area

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

  rsAdditiveKeyPatch* patch = nullptr;
  //rsAdditiveKeyPatchPlayable *patch = nullptr; // ..todo

  int nextBreakpointIndex     = 0;
  int samplesToNextBreakpoint = 0;
  bool playing = false;

  bool alwaysReInit = false; 
  // For experimentation. If true, we will re-initialize the instantaneous phase and amplitude at 
  // each breakpoint to the value prescribed in the patch rather than starting from where the osc
  // currently is. It's probably not desirable to ever set this to true as it will just introduce
  // discontinuities in cases of severe numeric drift. If this turns out to be the case, this 
  // field can be removed.

};

