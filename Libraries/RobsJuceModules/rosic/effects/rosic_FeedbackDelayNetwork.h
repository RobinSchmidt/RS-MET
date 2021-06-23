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
  // ToDo: have 3 diffusion coeffs that can be used to compute Euler angles for a 3x3 rotation
  // matrix that will be used in a 3x3 kronecker tarnsform (currently we only provide the 2x2
  // version)
  // maybe rename to density, "diffusion" may refer to some allpass pre- or post-process

  double a, b;
  // a, b values for the generalized Hadamard transform. ToDo: use a,b,c,d as members and assign
  // either c = -b, d = a or c = b,d = -a depending on a user switch. maybe use rsMatrix2x2

  RAPT::rsMatrix3x3<double> seedMat3x3;

  bool use3x3 = false; // switches betwen 2x2 and 3x3 seed matrices

  // todo: use std::vector for all the arrays

  //-----------------------------------------------------------------------------------------------
  // \name Deprecated

public:


   /** The "work" pointer should point to an array (of the same length as x) that can be used as 
   internal workspace.  
   \todo move this function out of this class - it might be useful in other contexts as well */
  RS_DEPRECATED( // replacement: RAPT::rsFGHT
  static void fastGeneralizedHadamardTransform(double* x, int N, int log2N, double* work,
    double a = 1.0, double b = 1.0, double c = 1.0, double d = -1.0));


  RS_DEPRECATED_WITH_BODY(
    static void fastInverseGeneralizedHadamardTransform(
      double* x, int N, int log2N, double* work,
      double a = 1.0, double b = 1.0, double c = 1.0, double d = -1.0),
    { RAPT::rsIFGHT(x, N, a,b,c,d); } )   // this is the replacement

};

}

#endif 
