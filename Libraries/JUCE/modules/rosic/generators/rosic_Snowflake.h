#ifndef rosic_Snowflake_h
#define rosic_Snowflake_h

namespace rosic
{

/** Snowflake is a turtle graphics based sound generator that uses Lindenmayer systems to generate
the turtle commands. The resulting curves show a self-similar character. The Koch snowflake is one 
simple example of such a curve. 

References:
PGCA: Pattern Generation for Computational Art (Stefan and Richard Hollos)
LSFP: Lindenmayer Systems, Fractals and Plants (Prusinkiewicz, Hanan)
ABoP: The Algorithmic Beauty of Plants (Prusinkiewicz, Lindenmayer)  */

class Snowflake : public TurtleSource
{

public:

  Snowflake();

  /** Clears the set of L-system rules. */
  void clearRules();

  /** Adds an L-system rule. */
  void addRule(char input, const std::string& output);

  /** Sets the seed (aka "axiom") for the L-system. */
  void setAxiom(const std::string& newAxiom);

  /** Sets the number of iterations for the L-system. */
  void setNumIterations(int newNumIterations);

  /** Updates all internal variables, so they reflect the user settings. Client code normally 
  doesn't have to care, because it does this automatically, when necessary. But sometimes it's
  useful for debugging to call this manually. */
  void updateAllInternals();

  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Calculates one output-sample frame at a time. */
  INLINE void getSampleFrameStereo(double *outL, double *outR)
  {
    if(!commandsReady)  
      updateTurtleCommands();
    TurtleSource::getSampleFrameStereo(outL, outR);
  }

protected:

  /** Updates our string of turtle-graphics drawing commands. */
  void updateTurtleCommands();

  LindenmayerSystem lindSys;
  std::string axiom;
  int numIterations = 0; // replace by iteratorString or applicatorString (a string like AAABBAC)
  std::string lindenmayerResult;          // output string of Lindenmayer system
  std::atomic_bool commandsReady = false; // flag to indicate that "turtleCommands" is up to date

};

}

#endif