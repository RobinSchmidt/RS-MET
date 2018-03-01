#include "FiniteAutomaton.h"

std::string FiniteAutomaton::readSymbol(char c)
{

  return std::string(); // return empty string by default
}

std::string FiniteAutomaton::translate(const std::string& input)
{
  std::string s;
  for(size_t i = 0; i < input.size(); i++)
    s += readSymbol(input[i]);
  return s;
}