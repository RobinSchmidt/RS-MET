#ifndef rosic_DspScriptInterpreter_h
#define rosic_DspScriptInterpreter_h

//#include "rosic_ExpressionEvaluator.h"

namespace rosic
{

  /**

  This is a class for ..

  */

  class DspScriptInterpreter : public ExpressionEvaluator
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    DspScriptInterpreter();    
    //< Constructor.
    
    //---------------------------------------------------------------------------------------------
    // setup:

    void setSampleRate(double newSampleRate);
    /**< Sets up the samplerate. */

    void initFunctionList();
    /** Clears the function list and thereafter adds the standard constants. */


  protected:

    double sampleRate;

  };

} // end namespace rosic

#endif // rosic_DspScriptInterpreter_h