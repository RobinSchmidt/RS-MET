#ifndef FiniteAutomaton_h
#define FiniteAutomaton_h

//#include "rosic/rosic.h"
#include <vector>
#include <string>

/** */

class FiniteAutomaton
{

public:


  //void setState


  std::string readSymbol(char c);

  std::string translate(const std::string& input);

  void reset() { state = 0; }

protected:

  std::vector<std::string> stateNames;
  std::vector<std::string> stateOutputs;
  std::vector<char> alphabet;
  //rsMatrix<int> transitions;
  int numStates = 0;
  int numSymbols = 0;

  int state = 0;

};


#endif