#ifndef rosic_FeedbackDelayNetwork_h
#define rosic_FeedbackDelayNetwork_h

namespace rosic
{

/** Under Construction...maybe move to prototypes

This class ...

\todo try a different ordering of the delaylines using the shuffling function in the cpp file.
this is supposed to be better for the diffusion parameter: short delaylines would mostly
crossfeed long ones and vice versa. if this is used, shuffle the output gains in the same way
(-> introduce output vectors as members and throw the same function at them) ...hmm but maybe
not - it seems that shuffled delayline-length together with non-shuffled output vector work well
- especiall with low diffusion

\todo introduce modulation of the angle "phi" by rectifying the FDNs output signal (maybe use
the mono-sum), pass it into a pair of filters tuned to the modulation-frequency and are 90
degrees out of phase (use the ModalFilter class for this). Two out-of-phase signals are needed in
order to obatin the instantaneous envelope of the modulation signal in order to divide by it to
obtain an amplitude-normalized modulation signal from either of the two filter outputs (perhaps
use the sine-component). */

class FeedbackDelayNetwork
{

public:


  /** \name Construction/Destruction */

  /** Constructor. */
  FeedbackDelayNetwork();

  /** Destructor. */
  ~FeedbackDelayNetwork();


  /** \name Setup */


  void setSampleRate(double newSampleRate);


  /** Sets the number of delaylines to be used for the FDN. It must be a power of two. */
  void setNumDelayLines(int newNumDelayLines);

  /** Sets the amount of diffusion/scattering between the delaylines. The value should be set in
  percent. */
  void setDiffusion(double newDiffusion);




  /** \name Audio Processing */


  /** Calculates a stereo output sample frame. */
  void processFrame(double* inOutL, double* inOutR);
    // \todo: inline this later

  //void processBlock(double *inOutL, double *inOutR, int blockSize);






  /** \name Others */


  /** Resets the internal state of the FDN (empties all delaylines, resets the damping filters,
  etc.). */
  void reset();

  /** Performs a generalized fast Hadamard transform on the input vector x with seed-matrix
  values a, b, c, d. For more details, refer to my paper "A Generalization of the Hadamard
  Transform". When the default values are used, it reduces to the standard Hadamard
  transform without scaling or sequency based ordering. The result will again end up in
  x. N is supposed to be a power of two that gives the vector's dimensionality and log2N should
  be the base-2 logarithm of N. The "work" pointer should point to an array (of the same length
  as x) that can be used as internal workspace.
  \todo move this function out of this class - it might be useful in other contexts as well */
  static void fastGeneralizedHadamardTransform(double* x, int N, int log2N, double* work,
    double a = 1.0, double b = 1.0, double c = 1.0, double d = -1.0);
  // todo: maybe try a generalized Hadamard trafo with complex coeffs
  // oh - and figure out for which choices of a,b,c,d the resulting matrix is unitary - will it 
  // be when using a matrix with c=-b, d=a or c=b, d=-a - oh - i checked with the 2x2 seed matrix
  // it is unitary indeed does this immply the higher order matrices are unitary too?
  // idea - try arbitrary a,b,c,d and rescale the output to have the same length as the input
  // vector - this makes the system nonlinear (really?) - maybe in an interesting way?


  //static void fastInverseGeneralizedHadamardTransform(double* x, int N, int log2N, double* work,
  //  double a = 1.0, double b = 1.0, double c = 1.0, double d = -1.0);

  RS_DEPRECATED_WITH_BODY(
  static void fastInverseGeneralizedHadamardTransform(
    double* x, int N, int log2N, double* work,
    double a = 1.0, double b = 1.0, double c = 1.0, double d = -1.0),
  { RAPT::rsIFGHT(x, N, a,b,c,d); } )   // this is the replacement


  // JUCE_DEPRECATED_WITH_BODY (virtual bool shouldDropFilesWhenDraggedExternally (const String&, Component*, StringArray&, bool&), { return false; })

protected:





  /** Sets up the relative delay-times according to the selected algorithm. */
  void setupRelativeDelayTimes();

  /** Sets up the indices where we read from the delayline */
  void setupReadIndices();

  /** Allocates all the required memory. */
  void allocateMemory();

  /** Allocates the memory fro the delaylines. */
  void allocateDelayLines();

  /** Frees all memory. */
  void freeMemory();

  /** Frees the memory for the delaylines. */
  void freeDelayLines();




  double** delayLines;     // the delaylines themselves
  int* readIndices;        // sample-indices where we read from the delaylines
  int* writeIndices;       // sample-indices where we write into the delaylines
  int* delayLineLengths;   // lengths of the delaylines in samples
  //double *outputGains;
  int    numDelayLines;    // must be a power of 2 and >= 4
  double sampleRate;


  double referenceDelayTime;  // the delaytime of the reference delayline (in seconds)
                              // determines perceived room-size
  double* relativeDelayTimes;

  double diffusion;
  double a, b;                // a, b values for the generalized Hadamard transform


  // todo: use std::vector for all the arrays
};






}

#endif 
