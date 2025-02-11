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

void chebyPolyConvert
{
  static const int N = 5;

  double a[N+1] = {9, 6, -10, -20, 24, 32}; // polynomial coeffs, gives the Chebychev expansion
                                            // 13*T0 + 11*T1 + 7*T2 + 5*T3 + 3*T4 + 2*T5
  int k;
  int i;
  int s;      // recursion stage

  // For test, we implement the algorithm in a way that stores all the intermediate arrays
  // in a 2D array:

  double B[N][N+1];
  memset(B, 0, N*(N+1)*sizeof(double));
  s = 0;      // recursion stage
  B[s][0] = a[N];
  B[s][1] = a[N-1]; // intitialization
  for(k = N-2; k >= 0; k--)
  {
    s++;
    B[s][0] = a[k]      + 0.5*B[s-1][1];
    B[s][1] = B[s-1][0] + 0.5*B[s-1][2];
    for(i = 2; i <= s+1; i++)
    {
      if( i < s )
        B[s][i] = 0.5*(B[s-1][i-1] + B[s-1][i+1]);
      else
        B[s][i] = 0.5*B[s-1][i-1];
    }


    /*
    for(i = 2; i < s; i++)
      B[s][i] = 0.5*(B[s-1][i-1] + B[s-1][i+1]);
    B[s][i] = 0.5*B[s-1][i-1];  // i == s here
    i++;
    B[s][i] = 0.5*B[s-1][i-1];  // i == s+1 here
    */
  }
  rsArrayTools::copy(B[N-1], b, N+1);
  // looks plausible and seems to work in this case


  //int dummy = 0;

  /*
  // reverse the algorithm: loops run backwards, order of instructions reversed, increments become
  // decremets, left-hand sides and right-hand sides of assignments exchange roles, when the
  // right-hand side contains a combination (i.e. sum) we have to look at the equation, find which
  // values are already known at this point and solve for the unknown which becomes the new
  // left-hand side:
  double C[N][N+1]; // we use C to reconstruct the B matrix
  memset(C, 0, N*(N+1)*sizeof(double));
  rsCopyBuffer(b, C[N-1], N+1); // the last stage is given, we must recostruct previous stages
  s = N-1;
  for(k = 0; k <= N-2; k++)
  {
    C[s-1][s] = 2*C[s][s+1];
    for(i = s; i >= 2; i--)
      C[s-1][i-1] = 2*C[s][i] - C[s-1][i+1];
    C[s-1][0] = C[s][1] - 0.5*C[s-1][2];
    c[k] = C[s][0] - 0.5*C[s-1][1];
    s--;
  }
  c[N-1] = C[s][1];
  c[N]   = C[s][0];
  */

  // \todo: get rid of the 2D-array to store all the intermediate stages - we may re-use a 1D array
  // at each stage when using 1 (or maybe 2) temporary variable
  // todo: move this commented stuff to the "Experiments" project - we may want to have it
  // available for later reference
}