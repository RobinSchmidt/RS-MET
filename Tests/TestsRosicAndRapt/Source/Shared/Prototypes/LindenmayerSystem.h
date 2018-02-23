#ifndef LindenmayerSystem_h
#define LindenmayerSystem_h

//#include <vector>
//#include <string>
#include "rosic/rosic.h"

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

  /** Clears the set of replacement rules. */
  void clearRules();

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

todo: maybe make a 3D version

References:
https://en.wikipedia.org/wiki/Turtle_graphics  */

class TurtleGraphics
{

public:

  /** Creates a new TurtleGraphics object using the given angle (in degrees, default is 90°) */
  TurtleGraphics(double angle = 90) { setAngle(angle); }

  /** Sets the turning angle in degrees. */
  void setAngle(double degrees);

  /** Initializes position (x,y) and direction vector (dx,dy). */
  void init(double x, double y, double dx, double dy);

  /** Moves the turtle one step into the current direction. */
  void goForward();

  /** Turns the direction vector to the left. */
  void turnLeft();

  /** Turns the direction vector to the right. */
  void turnRight();

  /** Translates the given string into arrays of x,y coordinates of vertices. */
  void translate(const std::string& str, std::vector<double>& vx, std::vector<double>& vy);

protected:

  double  x = 0,  y = 0;  // current position
  double dx = 1, dy = 0;  // direction vector

  RAPT::rsRotationXY<double> rotLeft, rotRight;

};

//=================================================================================================

/** Encapsulates a LindenmayerSytem (as baseclass) and a TurtleGraphics object (as member) to make
rendering of L-system generatad point-sequences more convenient. */

class LindenmayerRenderer : public LindenmayerSystem
{

public:

  /** Sets the turning angle for the turtle graphics renderer. */
  void setAngle(double degrees) { turtleGraphics.setAngle(degrees); }

  /** Turns normalization on/off. The normalization also includes DC removal. */
  void setNormalization(bool shouldNormalize) { normalize = shouldNormalize; }




  /** Produces the array of verices for a Koch snowflake of given order.
  see: https://en.wikipedia.org/wiki/Koch_snowflake */
  void getKochSnowflake(int order, std::vector<double>& x, std::vector<double>& y);

  /** Produces the array of verices for a Moore curve of given order. 
  see: https://en.wikipedia.org/wiki/Moore_curve */
  void getMooreCurve(int order, std::vector<double>& x, std::vector<double>& y);

  // some more curves taken from: http://mathforum.org/advanced/robertd/lsys2d.html

  void get32SegmentCurve(int order, std::vector<double>& x, std::vector<double>& y);
  void getQuadraticKochIsland(int order, std::vector<double>& x, std::vector<double>& y);

  void getSquareCurve(int order, std::vector<double>& x, std::vector<double>& y);


  // maybe have a "numPoints" parameter that is used for resampling the resulting curve to a given
  // number of points (by linear interpolation)

  void render(const std::string& seed, double angle, int order, 
    std::vector<double>& x, std::vector<double>& y);

  /** Renders the given Lindenmayer string into a sequence of points using the given turning 
  angle. */
  void translate(const std::string& str, double angle, 
    std::vector<double>& x, std::vector<double>& y);

  /** Normalizes the xy coordinates such that both x and y are free of DC (centered around 0) and
  the maximum absolute value is 1 (this second step is done by finding the maximum absolute value
  of both vectors and scaling both by the same value - in order to preserve the aspect ratio). */
  void normalizeXY(std::vector<double>& x, std::vector<double>& y);

protected:

  TurtleGraphics turtleGraphics;

  bool normalize = true;

};


#endif