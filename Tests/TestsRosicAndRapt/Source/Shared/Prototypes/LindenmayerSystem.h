#ifndef LindenmayerSystem_h
#define LindenmayerSystem_h

#include <vector>
#include <string>

/** Implements a Lindenmayer system (aka L-system). This is a system that iteratively replaces 
certain characters in a string by replacement strings according to a set of replacement rules. 
When the resulting strings are translated to graphics by interpreting some of the characters as 
drawing commands (such as 'F': go forward, '+': turn right, '-': turn left), curves with 
self-similar features can be generated.

References:
https://de.wikipedia.org/wiki/Lindenmayer-System
https://en.wikipedia.org/wiki/L-system
http://algorithmicbotany.org/papers/abop/abop-ch1.pdf
http://www.abrazol.com/books/patterngen/  

*/

class LindenmayerSystem
{
  
public:

  /** Adds a replacement rule to our set of rules.  */
  void addRule(char input, const std::string& output);

  /** Applies the set of replacement rules to the given character. If a matching rule-input 
  character is found, it returns the corresponding output string, otherwise the character itself
  will be returned as string of length 1. */
  std::string apply(char c);

  /** Applies the set of replacement rules to the given string once and returns the result. */
  std::string apply(const std::string& input);

  /** Iteratively applies the set of replacement rules the given number of times and returns the 
  result. */
  std::string apply(const std::string& input, int numberOfTimes);


protected:

  std::vector<char> ruleInputs;
  std::vector<std::string> ruleOutputs;

};





#endif