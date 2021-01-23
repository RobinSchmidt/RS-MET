#ifndef rosic_SineOscillator_h
#define rosic_SineOscillator_h

namespace rosic
{

  /**

  This class implements a sine oscillator based on a self oscillating (free running) second order 
  recursive filter. 

  */

  class SineOscillator  
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    SineOscillator();   

    /** Destructor */
    ~SineOscillator();  

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the frequency (in Hz) of the sine to be generated. */
    void setFrequency(double newFrequency);

    /** Sets the normalized radian frequency (=2*pi*frequency/sampleRate) of the sine to be 
    generated. */
    void setOmega(double newOmega);

    /** Sets the start phase in radians - this is the phase to wich the oscillator retriggers
    on calls to trigger(). */
    void setStartPhase(double newStartPhase);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the frequency (in Hz) of the sine. */
    double getFrequency() const { return frequency; }

    /** Returns the normalized radian frequency (=2*pi*frequency/sampleRate) of the sine. */
    double getOmega() const { return omega; }

    /** Returns the start phase in radians - this is the phase to wich the oscillator retriggers
    on calls to trigger(). */
    double getStartPhase() const { return startPhase; }

    //---------------------------------------------------------------------------------------------
    // event handling:

    /** Resets the oscillator to its start phase. */
    void trigger(); // rename to reset...maybe

    /** Triggers the oscillator with an arbitrary start phase. */
    void triggerWithPhase(double phase);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output-sample at a time with acquiring the mutex-lock before. */
    INLINE double getSampleThreadSafe();

    /** Calculates one output-sample at a time. */
    INLINE double getSample();

    //=============================================================================================

  protected:

    double a1;         /**< The filter coefficient for the recursion. */
    double s1, s2;     /**< Past sine outputs (delayed by one or two samples respectively). */
    double omega;      /**< Normalized radian frequency of the oscillator. */
    double startPhase; /**< The phase to which the oscillator is reset on trigger(). */

    // these may be factored out into a baseclass:
    double frequency;  /**< The frequency of the oscillator in Hz. */
    double sampleRate; /**< The sample-rate. */

    MutexLock mutex;   
    // to ensure consistency between the variables (otherwise harsh noise may result)
    // todo: get rid of this - maybe the class should be considered obsolete anyway due to having
    // similar functionality in rapt now

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):


  INLINE double SineOscillator::getSampleThreadSafe()
  {
    mutex.lock();
    double result = getSample();
    mutex.unlock();
    return result;
  }

  INLINE double SineOscillator::getSample()
  {
    double tmp = a1*s1 - s2;
    s2         = s1;
    s1         = tmp;
    return tmp;
  }

} // end namespace rosic

#endif // #ifndef rosic_SineOscillator_h
