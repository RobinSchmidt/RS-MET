#ifndef rosic_AngelScriptInterpreter_h
#define rosic_AngelScriptInterpreter_h

#include "../_third_party/angelscript/angelscript.h"


namespace rosic
{

  /**

  This is a class for ..

  */

  class AngelScriptInterpreter 
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    AngelScriptInterpreter();    
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

#endif // rosic_AngelScriptInterpreter_h