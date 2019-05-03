#ifndef rosic_RoutingMatrix_h
#define rosic_RoutingMatrix_h

//// rosic-indcludes:
//#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

  /**

  This class implements a mixing matrix that can mix a number of inputs (signal sources) to a 
  number of outputs (signal destinations).

  */

  class RoutingMatrix
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. You must pass the number of inputs and outputs here. */
    RoutingMatrix(int numInputs, int numOutputs); 

    /** Destructor. */
    ~RoutingMatrix();  

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the amount by which one input is mixed into some output (as raw factor). */
    void setMatrixEntry(int inputIndex, int outputIndex, double value);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the amount by which one input is mixed into some output (as raw factor) - if one of
    the indices is out of range, it will return 0.0. */
    double getMatrixEntry(int inputIndex, int outputIndex);

    /** Returns the amount by which one input is mixed into some output (as raw factor) - if one of
    the indices is out of range, it will cause an access violation - no bounds checking is done 
    here, making access faster but less safe. */
    double getMatrixEntryFast(int inputIndex, int outputIndex) 
    { return mixMatrix[inputIndex][outputIndex]; }

    /** Returns the number of input signals that go into the mixer. */
    int getNumInputs() const { return numInputs; }

    /** Returns the number of output signals that go out of the mixer. */
    int getNumOutputs() const { return numOutputs; }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Computes the outputs from the array 'inputs' and returns them in the array 'outputs'.
    The arrays are assumed to be numOutputs and numInputs (as set in the constructor) long. */
    INLINE void computeOutputs(double *inputs, double *outputs);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Initializes the matrix with all zeros. */
    void initializeMatrix();

    //=============================================================================================

  protected:

    int      numInputs, numOutputs;
    double** mixMatrix;
    double*  flatArray;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void RoutingMatrix::computeOutputs(double *inputs, double *outputs)
  {
    for(int o=0; o<numOutputs; o++)
    {
      outputs[o] = 0.0;
      for(int i=0; i<numInputs; i++)
        outputs[o] += mixMatrix[i][o] * inputs[i];
    }
  }

} // end namespace rosic

#endif // rosic_RoutingMatrix_h
