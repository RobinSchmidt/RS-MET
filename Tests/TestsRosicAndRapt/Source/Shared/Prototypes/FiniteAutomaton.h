#ifndef FiniteAutomaton_h
#define FiniteAutomaton_h

//#include "rosic/rosic.h"
#include <vector>
#include <string>

/** Implements a finite automaton (aka state machine?). */

class FiniteAutomaton
{

public:


  //void setState

  /** Returns a string representation of this automaton suitable for storing the automaton 
  definition in a file. */
  void getAsString();

  /** Sets this automaton up from a defining string. */
  void setupFromString(const std::string& definition);

  /** Perform the state transition associated with the given symbol and returns the output string 
  that is associated with the new state. */
  std::string readSymbol(char c);

  /** Translates the given string input string to the corresponding output string. This can also be 
  called in the middle of another translation without affecting our state. */
  std::string translate(const std::string& input);

  /** Resets the current state into the start state (index 0). */
  void reset() { state = 0; }

protected:

  std::vector<std::string> stateNames;
  std::vector<std::string> stateOutputs;
  std::vector<char> alphabet; // input alphabet
  //rsMatrix<int> transitions;
  int numStates = 0;
  int numSymbols = 0;

  int state = 0;

};


#endif