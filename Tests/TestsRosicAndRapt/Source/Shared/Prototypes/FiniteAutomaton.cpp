#include "FiniteAutomaton.h"

std::string FiniteAutomaton::readSymbol(char c)
{
  // -find index i belonging to c (maybe have a size 255 array of ints with entries that
  //  are all -1 by default)
  // -if i is in the valid range:
  //  -read out the transitions matrix at [state][i] to update state: 
  //   state = transitions[state][i];
  //  -return stateOutputs[i]

  return std::string(); // return empty string by default
}

std::string FiniteAutomaton::translate(const std::string& input)
{
  int tmpState = state;

  state = 0;
  std::string s;
  for(size_t i = 0; i < input.size(); i++)
    s += readSymbol(input[i]);

  state = tmpState; // we dont want a call to translate affect the state member
  return s;
}