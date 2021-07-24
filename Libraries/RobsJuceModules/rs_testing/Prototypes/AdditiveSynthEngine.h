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

template<class T, int N>
class rsSineSweeperBank
{

public:

  //using CVec = const rsSimdVector<T, N>;

  using Parameters = RAPT::rsSweepParameters<rsSimdVector<T, N>>;

  virtual void setMaxNumOscillators(int newLimit) = 0;

  virtual void setNumActiveGroups(int newNumber) = 0;

  virtual void setup(int simdGroupIndex, const Parameters& params) = 0;

  virtual int getNumActiveGroups() const = 0;

  /** Should return instantaneous phase. */
  virtual rsSimdVector<T, N> getPhase(int groupIndex) const = 0;

  /** Should return instantaneous amplitude. */
  virtual rsSimdVector<T, N> getAmplitude(int groupIndex) const = 0;


  virtual void processFrame(float* left, float* right) = 0;
  // maybe have a double variant, too

  virtual void reset() = 0;

};
// i think this needs a template parameter N for th simd-size, then init should take simd-vector
// parameters and the "index" parameter can go away

//=================================================================================================

/** Implements a bank of sinusoidal oscillators with cubic envelopes for instantaneous phase 
and for the instantaneous amplitude. It's based on maintaining a phasor and computing the sine 
directly either exactly or by various polynomial approximations delivering different amounts of 
fidelity and eating different amounts of CPU time ...tbc... */

template<class T, int N>
class rsSineSweeperBankDirect : public rsSineSweeperBank<T, N>
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
class rsSineSweeperBankIterative : public rsSineSweeperBank<T, N>
{

public:

  using Parameters = RAPT::rsSweepParameters<rsSimdVector<T, N>>;

  void setMaxNumOscillators(int newLimit) override;


  void setNumActiveGroups(int newNumber) override;

  void setup(int i, const Parameters& p) override 
  { 
    simdGroups[i].setup(p); 
  }


  int getNumActiveGroups() const override { return (int) simdGroups.size(); }

  virtual rsSimdVector<T, N> getPhase(int i) const override
  { return simdGroups[i].getPhase(); }

  virtual rsSimdVector<T, N> getAmplitude(int i) const override
  { return simdGroups[i].getAmplitude(); }


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

public:

  void setShape(int newDim1, int newDim2)
  {
    dim1 = newDim1;
    dim2 = newDim2;
    data.resize(dim1*dim2);
  }
  // todo: optionally retain old data (true by default - rationale: it's safer and optimizations 
  // that sacrifice safety should be explicit)

  void fill(const T& value)
  {
    RAPT::rsArrayTools::fillWithValue(&data[0], (int) data.size(), value);
  }

  int getSize() const { return dim1 * dim2; }

  void memsetAllZero()
  {
    memset(&data[0], 0, sizeof(T) * (size_t)getSize()); 
  }


  /** Read and write access to elements with index pair i,j. */
  T& operator()(const int i, const int j) { return data[flatIndex(i, j)]; }

  const T& operator()(const int i, const int j) const { return data[flatIndex(i, j)]; }

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

/** Single voice for the additive synthesis engine based on oscillator banks using SIMD 
processing. It the plural "banks" because different implementations are available with different 
tradeoffs with respect to accuracy and efficiency. */

template<int N>             // N: size of the SIMD vectors use NSimd
class rsAdditiveSynthVoice
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Classes

  /** Data structure to represent a patch for a particular key (midi note). */
  struct EditablePatch  // maybe turn into class
  {
    // ToDo: use double for all parameters

    /** Represents the parameters in one of two units, depending on the context. One set of units
    is for presenting them to the user and the other is for internal use in the algorithm. */
    struct SineParams
    {
      SineParams(){}
      SineParams(double freq_, double gain_, double phase_ = 0) 
        : freq(freq_), gain(gain_), phase(phase_) { }

      double freq  = 0;  // in Hz (user), as omega = 2*pi*fs/fs (algo)
      double gain  = 0;  // as raw factor (user, direct algos), log of factor (iterative algos)
      double phase = 0;  // in degrees (user), radians (algo), some algos may ignore it
      // double fade = 0; // derivative of gain
    };

    struct Breakpoint  // maybe mae it a class
    {
      double time = 0;   // in seconds (user), samples (algo)
      void addSine(const SineParams& newSine) { params.push_back(newSine); }
      int getNumPartials() const { return (int) params.size(); }
      std::vector<SineParams> params;  // rename to partials or sines
    };

    void addBreakpoint(Breakpoint bp) { breakpoints.push_back(bp); }

    void clear() { breakpoints.clear(); }

    /** Fills the phase values of all breakpoints with artificial values obtained from integrating 
    the instantaneous frequency. */
    void createArtificialPhases();

    //void convertUserToAlgoUnits(float sampleRate, bool logAmplitude);
    // seems superfluous now - we should directly convert into a PlayablePatch


    int getNumBreakpoints() const { return (int) breakpoints.size(); }

    int getNumPartials() const;

    bool isWellFormed() const;
    // -time-stamps of breakpoints miust start at zero and be strictly increasing
    // -all breakpoints must have the same number of partials

    const Breakpoint* getBreakpoint(int i) const { return &breakpoints[i]; }

    Breakpoint* getBreakpoint(int i) { return &breakpoints[i]; }

  protected:

    std::vector<Breakpoint> breakpoints;

  };

  /** This also represents a patch for a particular key like EditablePatch, but has a different 
  data layout to facilitate efficient playback. But this layout impedes efficient editing, so for 
  editing patches, EditablePatch is used and once a patch is ready, it can be converted into a
  "playable" version in one go as preprocessing via setupFrom. */
  class PlayablePatch
  {

  public:

    /** Sets up this object from the given editable patch. ...tbc... */
    void setupFrom(const EditablePatch& patch, double sampleRate, bool applyLogToAmplitude);


    int getNumBreakpoints() const { return numBreakpoints; }

    int getNumSimdGroups() const { return numSimdGroups; }

    int getBreakpointTime(int i) const { return timeStamps[i]; }

    const RAPT::rsSweepParameters<rsSimdVector<float, N>>& getParams(int i, int j) const
    {
     
      return params(i, j);
    }


  protected:

    int numBreakpoints = 0;       // number of breakpoints   
    int numSimdGroups  = 0;       // number of SIMD-groups of partials
    std::vector<int> timeStamps;  // time-stamps of the breakpoints in samples
    rsArray2D<RAPT::rsSweepParameters<rsSimdVector<float, N>>> params;
    // 1st index: breakpoint, 2nd: simd-group of partials

  };

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  void setPatch(PlayablePatch* newPatch) { patch = newPatch; }

  void setSweeperBank(rsSineSweeperBank<float, N>* newBank) { sweeperBank = newBank; }

  /** Selects whether the self-correction of instantaneous amplitude and phase to counteract 
  numerical drift should be smooth (default) or not. The non-smooth way is mostly for testing and 
  development. In a real world application, there's probably no good reason to not want it to be 
  smooth. */
  void setSmoothCorrection(bool smooth) { alwaysReInit = !smooth; }

  // void setPhaseless(bool); 
  // todo: 
  // -Implement phaseless synthesis 
  //  -is cheaper (needs only a quadratic polynomial iterator)
  //  -phase can't be controlled directly - it just is what it is
  //  -amplitude correction may be more complicated, we'll see

  //-----------------------------------------------------------------------------------------------
  // \name Processing

  void startPlaying();

  void handleBreakpoint(int index, bool reInitAmpAndPhase); 
  // move to protected

  void processFrame(float* left, float* right);

  void reset();

protected:

  void initSweepers(int startBreakpointIndex, int endBreakpointIndex, bool reInitAmpAndPhase);


  rsSineSweeperBank<float, N>* sweeperBank = nullptr;

  PlayablePatch* patch = nullptr;

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

