#ifndef rosic_Snowflake_h
#define rosic_Snowflake_h

namespace rosic
{

/** Snowflake is a turtle graphics based sound generator that uses Lindenmayer systems to generate
the turtle commands. The resulting curves show a self-similar character. The Koch snowflake is one
simple and famous example of such a curve and gave this generator its name.

References:
PGCA: Pattern Generation for Computational Art (Stefan and Richard Hollos)
LSFP: Lindenmayer Systems, Fractals and Plants (Prusinkiewicz, Hanan)
ABoP: The Algorithmic Beauty of Plants (Prusinkiewicz, Lindenmayer)  */

class Snowflake : public TurtleSourceAntiAliased
{

public:

  Snowflake();

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Clears the set of L-system rules. */
  void clearRules();

  /** Adds an L-system rule. */
  void addRule(char input, const std::string& output);

  /** Sets the seed (aka "axiom") for the L-system. */
  void setAxiom(const std::string& newAxiom);

  /** Sets the number of iterations for the L-system. */
  void setNumIterations(int newNumIterations);

  /** Sets the amount of detail to be prduced in the line-drawing from 0 to 100 percent. 0 uses
  and L-system with 0 iterations, 100 an L-system with an iteration number equal to what is set up
  by setNumIterations. Non-integer values will result in a crossfade between 2 numbes of
  iterations. For example, when setNumIterations is set to 10, detail to 35%, you'll get a 50/50
  corssfade between 3 and and 4 iterations. */
  void setDetailAmount(double newAmount);
   // not yet implemented - it's complicated - we may need to always run at least 2 turtles, maybe
   // even numIterations turtles....

  /** Updates all internal variables, so they reflect the user settings. Client code normally
  doesn't have to care, because it does this automatically, when necessary. But sometimes it's
  useful for debugging to call this manually. */
  void updateAllInternals();


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the number of lines that the trutel would produce according to the current command
  string. */
  int getNumTurtleLines() override;

  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Calculates one output-sample frame at a time. */
  INLINE void getSampleFrameStereo(double *outL, double *outR)
  {
    if(!commandsReady)
      updateTurtleCommands();
    TurtleSourceAntiAliased::getSampleFrameStereo(outL, outR);
  }

protected:

  /** Updates our string of turtle-graphics drawing commands. */
  void updateTurtleCommands();

  LindenmayerSystem lindSys;
  std::string axiom;     // maybe rename to seed
  int numIterations = 0; // replace by iteratorString or applicatorString (a string like AAABBAC)
  std::string lindenmayerResult;          // output string of Lindenmayer system
  //std::atomic_bool commandsReady = false; // flag to indicate that "turtleCommands" is up to date
  std::atomic_bool commandsReady; // flag to indicate that "turtleCommands" is up to date

};

}

#endif
