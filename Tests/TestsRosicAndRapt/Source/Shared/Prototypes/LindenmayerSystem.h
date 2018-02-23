#ifndef LindenmayerSystem_h
#define LindenmayerSystem_h

#include <vector>
#include <string>

/** Implements a Lindenmayer system (aka L-system). This is a system that iteratively replaces 
certain characters in a string by replacement strings according to a set of replacement rules. 
When the resulting strings are translated to graphics by interpreting some of the characters as 
drawing commands (such as 'F': go forward, '+': turn left, '-': turn right), curves with 
self-similar features can be generated.

References:
https://de.wikipedia.org/wiki/Lindenmayer-System
https://en.wikipedia.org/wiki/L-system
http://algorithmicbotany.org/papers/abop/abop-ch1.pdf
http://www.abrazol.com/books/patterngen/  
http://www.sidefx.com/docs/houdini/nodes/sop/lsystem.html
http://www.sccg.sk/~smolenova/elearning/ks_fmfiuk06.pdf
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
  result. In L-system jargon, the "input" string is also called "axiom" and the numberOfTimes is 
  called order. */
  std::string apply(const std::string& input, int numberOfTimes);


protected:

  std::vector<char> ruleInputs;
  std::vector<std::string> ruleOutputs;

};

//=================================================================================================

/** Implements a "turtle graphics" drawer that can be used to translate strings produced by a
Lindenmayer system into 2D drawings (i.e. a sequence of points in the xy-plane).

References:
https://en.wikipedia.org/wiki/Turtle_graphics  */

class TurtleGraphics
{

public:

  void init(double x, double y, double dx, double dy);

  void goForward();

  void turnLeft();

  void turnRight();

  /** Translates the given string into arrays of x,y coordinates of vertices. */
  void translate(const std::string& str, std::vector<double>& vx, std::vector<double>& vy);

protected:

  double  x = 0,  y = 0;  // current position
  double dx = 1, dy = 0;  // direction vector

};


#endif