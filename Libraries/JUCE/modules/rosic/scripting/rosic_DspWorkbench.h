#ifndef rosic_DspWorkbench_h
#define rosic_DspWorkbench_h

//// rosic-indcludes:
////#include "../infrastructure/rosic_PresetRememberer.h"
////#include "../infrastructure/rosic_MutexLock.h"
//#include "../basics/rosic_TabulatedFunction.h"
//#include "../modulators/rosic_BreakpointModulator.h"
//#include "../filters/rosic_EllipticSubBandFilter.h"
//#include "../filters/rosic_MultiModeFilter.h"
//#include "rosic_DspScriptInterpreter.h" 
////#include "../math/rosic_ExpressionEvaluator.h"
////#include "../datastructures/rosic_String.h"
//#include <string>

namespace rosic
{

  /**

  This is a ...

  */

  class DspWorkbench //: public PresetRememberer
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    DspWorkbench();  ///< Constructor.

    ~DspWorkbench();  ///< Destructor.

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    void setSampleRate(double newSampleRate);
    /**< Sets the sample-rate for this module. */

    void setOversamplingFactor(int newOversamplingFactor);
    /**< Sets the oversampling factor for the internal signal processing - a factor of one (or 
    lower) will result in no oversampling. */

    bool setAlgorithmString(char* newAlgorithmString);
    /**< Sets a new string for the expression evaluator to be parsed - the string represents the 
    DSP algorithm in a syntax suitable for the ExprEval library. */

    void setParameter(int index, double newValue);
    /**< Assigns the parameter with the given index to a numeric value. */

    void setParameterName(int index, char* newName);
    /**< Assigns a name for the parameter with the given index. */

    //---------------------------------------------------------------------------------------------
    // inquiry:

    int getOversamplingFactor() const;
    /**< Returns the oversampling factor. */

    const char* getAlgorithmString();
    /**< Returns the current algorithm-string as c-string. */

    double getParameter(int index) const;
    /**< Returns the value of the parameter with the given index. */

    const char* getParameterName(int index) const;
    /**< Returns the name of the parameter with the given index as c-string. */

    const char* getInternalParameterName(int index) const;
    /**< Returns the internal name of the parameter (as used in the code) with the given index 
    as c-string. */

    //---------------------------------------------------------------------------------------------
    // audio processing:

    INLINE double getSample(double in);
    ///< Calculates one output sample at a time.

    INLINE void getSampleFrameStereo(double* inL, double* inR, double* outL, double* outR);
    ///< Calculates one streo output frame at a time.

    INLINE void applyAlgorithm(double* inOutL, double* inOutR);
    /**< Applies the actual DSP algorithm by invoking the expression evaluator. */

    //---------------------------------------------------------------------------------------------
    // others:

    //---------------------------------------------------------------------------------------------
    // embedded objects:

    static const int numModulators = 10;
    //static const int numParameters = 128;
    static const int numParameters = 4;

    BreakpointModulator   modulators[numModulators];
    MultiModeFilter       filters[10];
    DspScriptInterpreter  scriptInterpreter;
    EllipticSubBandFilter upsamplerL, upsamplerR, antiAliasFilterL, antiAliasFilterR;

    //=============================================================================================

  protected:

    void initializeArrays();
    /** This function initializes our arrays. */

    void setParameterInternalName(int index, char* newName);
    /**< Assigns a name for the parameter with the given index. */

    // parameter variables:
    double sampleRate;
    double parameters[numParameters];
    char*  parameterClearTextNames[numParameters];
    char*  parameterInternalNames[numParameters];
    int    oversamplingFactor;
    bool   expressionIsValid;

    std::string algorithm;

    // pointers to the input and output variables:
    double *inL, *inR, *outL, *outR;

    // a mutex lock for the exporession evaluator object:
    //MutexLock evaluatorLock;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions
  // which are supposed to be called at audio-rate (they can't be put into
  // the .cpp file):

  INLINE double DspWorkbench::getSample(double in)
  {
    return in;
  }

  INLINE void DspWorkbench::getSampleFrameStereo(double* inL, double* inR, 
    double* outL, double* outR)
  {
    static doubleA tmpL, tmpR;
    static intA i;

    //evaluatorLock.lock();

    tmpL = *inL;
    tmpR = *inR;
    if( oversamplingFactor > 1 )
    {
      // calculate n frames of the oversampled, distorted, and anti-alias filtered signal
      // (we do oversampling by a factor of n):
      tmpL = upsamplerL.getSampleDirect1(oversamplingFactor*tmpL);
      tmpR = upsamplerR.getSampleDirect1(oversamplingFactor*tmpR);
      applyAlgorithm(&tmpL, &tmpR);
      tmpL = antiAliasFilterL.getSampleDirect1(tmpL);
      tmpR = antiAliasFilterR.getSampleDirect1(tmpR);

      for(i=1; i<oversamplingFactor; i++)
      {
        tmpL = upsamplerL.getSampleDirect1(0.0);
        tmpR = upsamplerR.getSampleDirect1(0.0);
        applyAlgorithm(&tmpL, &tmpR);
        tmpL = antiAliasFilterL.getSampleDirect1(tmpL);
        tmpR = antiAliasFilterR.getSampleDirect1(tmpR);
      }
    }
    else
      applyAlgorithm(&tmpL, &tmpR);

    //evaluatorLock.unlock();

    // store the output in its slots:
    *outL = tmpL;
    *outR = tmpR;

    return;
  }
 
  INLINE void DspWorkbench::applyAlgorithm(double* inOutL, double* inOutR)
  {
    // assign variables in the expression evaluator object (this is possible because we have
    // pointers to their memory locations):
    *inL = *inOutL;
    *inR = *inOutR;
    
    // evaluate the expression:
    if( expressionIsValid )
      scriptInterpreter.evaluateExpression();
    else
      return;

    // read out the results by means of our pointers to the memory locations of the outputs:
    *inOutL = *outL;
    *inOutR = *outR;
  }

} // end namespace rosic

#endif // rosic_DspWorkbench_h