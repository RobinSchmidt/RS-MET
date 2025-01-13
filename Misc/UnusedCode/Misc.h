/**

Extends BasicDelayLine by keeping information about the sample rate and the delay in
seconds.

\todo: facilitate tempo-sync by maintaining a bpm-value and a sync-flag
...maybe do this in a subclass...

*/

template<class TSig, class TPar>
class rsDelayLine : public rsBasicDelayLine<TSig>
{

public:

  /** \name Construction/Destruction */

  /** Constructor - constructs a delay-line with a given maximum number of samples delay. This
  has to be a power of two minus 1 - otherwise the next power of two minus 1 will be used. */
  rsDelayLine();

  /** Destructor */
  ~rsDelayLine();


  /** \name Setup */

  /** Sets the sample-rate. */
  void setSampleRate(TPar newSampleRate);

  /** Sets the delay-time in samples. */
  void setDelayInSamples(int newDelayInSamples);

  /** Sets the delay-time in seconds. */
  void setDelayInSeconds(TPar newDelayInSeconds);

  /** Sets the delay-time in milliseconds. */
  void setDelayInMilliseconds(TPar newDelayInMilliseconds);


  /** \name Inquiry */

  /** Returns the delay-time in seconds. */
  RS_INLINE TPar getDelayInSeconds() const { return delayInSeconds; }

  /** Returns the delay-time in milliseconds. */
  RS_INLINE TPar getDelayInMilliseconds() const { return 1000.0 * delayInSeconds; }

protected:

  /** \name Data */

  TPar delayInSeconds;
  TPar sampleRate;

};



// Construction/Destruction:

template<class TSig, class TPar>
rsDelayLine<TSig, TPar>::rsDelayLine()
{
  delayInSeconds = TPar(0.001);
  sampleRate     = TPar(44100.0);
  setDelayInSeconds(delayInSeconds);
}

template<class TSig, class TPar>
rsDelayLine<TSig, TPar>::~rsDelayLine()
{

}

// Setup:

template<class TSig, class TPar>
void rsDelayLine<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  sampleRate = newSampleRate;
  setDelayInSeconds(delayInSeconds);
}

template<class TSig, class TPar>
void rsDelayLine<TSig, TPar>::setDelayInSamples(int newDelayInSamples)
{
  rsBasicDelayLine<TSig>::setDelayInSamples(newDelayInSamples);
  delayInSeconds = (TPar) newDelayInSamples / sampleRate;
}

template<class TSig, class TPar>
void rsDelayLine<TSig, TPar>::setDelayInSeconds(TPar newDelayInSeconds)
{
  delayInSeconds = rsMax(0.0, newDelayInSeconds);
  setDelayInSamples( rsRoundToInt(sampleRate*delayInSeconds) );
}

template<class TSig, class TPar>
void rsDelayLine<TSig, TPar>::setDelayInMilliseconds(TPar newDelayInMilliseconds)
{
  setDelayInSeconds(0.001*newDelayInMilliseconds);
}








//==================================================================================================