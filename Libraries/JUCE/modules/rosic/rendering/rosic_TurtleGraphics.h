#ifndef rosic_TurtleGraphics_h
#define rosic_TurtleGraphics_h

namespace rosic
{

/** Implements a "turtle graphics" drawer that can be used to translate strings produced by a
Lindenmayer system into 2D drawings (i.e. a sequence of points in the xy-plane).

References:
https://en.wikipedia.org/wiki/Turtle_graphics  

todo: maybe make a 3D version */

class TurtleGraphics
{

public:

  /** Creates a new TurtleGraphics object using the given angle (in degrees, default is 90°) */
  TurtleGraphics(double angle = 90) { setAngle(angle); }

  /** Sets the turning angle in degrees. */
  void setAngle(double degrees);

  /** Initializes position (x,y) and direction vector (dx,dy). */
  void init(double x = 0, double y = 0, double dx = 1, double dy = 0, bool keepOldCoordsAsStart = false);

  /** Moves the turtle one step into the current direction. */
  void goForward();

  /** Turns the direction vector to the left. */
  void turnLeft();

  /** Turns the direction vector to the right. */
  void turnRight();

  /** Pushes the current state of the turtle onto a stack for branching. */
  void pushState();

  /** Pops the most topmost state form the stack and sets the turtle state accordingly. */
  void popState();

  /** Translates the given string into arrays of x,y coordinates of vertices. */
  void translate(const std::string& str, std::vector<double>& vx, std::vector<double>& vy);

  /** Given a string that may contain drawing commands, this function returns a string that 
  contains only those drawing commands. It may be used to "clean up" a raw output string from a 
  Lindenmayer system, i.e. to strip off the symbols that would be ignored anyway. */
  std::string extractCommands(const std::string& str);

  /** Returns true, if the given character is a drawing command */
  bool isCommand(char c);

  /** Given a string that supposedly contains drawing commands, this function returns the number of
  output lines in the drawing, that this string would produce. It just counts the number of 
  occurrences of characters that draw a line, i.e 'F' and 'G'. */
  int getNumberOfLines(const std::string& commandString);

  /** Returns the current x-coordinate. */
  double getX() { return x; }

  /** Returns the current y-coordinate. */
  double getY() { return y; }

  double getStartX() { return xo; }
  double getStartY() { return yo; }
  double getEndX()   { return x; }
  double getEndY()   { return y; }
  // maybe we need getStartX, getEndX, getStartY, getEndY - returns start and end-points of most
  // recent line

  /** Interprets the command given by character c. This may lead to a change of our internal 
  state, i.e. either (x,y) is changed (in case of 'F' or 'f') or (dx,dy) is changed (in case of
  '+' or '-'). For other characters, nothing happens. The return value informs whether of not an
  edge should be drawn to the new point (i.e. whether the character was an 'F'). If so, you may 
  retrieve the new point via getX, getY and draw a line to it. */
  bool interpretCharacter(char c);


protected:

  double x  = 0, y  = 0;  // current position
  double dx = 1, dy = 0;  // direction vector
  double xo = 0, yo = 0;  // old position (from one step before)

                          /** Structure to store and retrieve the current state of the turtle.  */
  struct TurtleState
  {
    TurtleState(double _x, double _y, double _dx, double _dy, double _xo, double _yo)
      : x(_x), y(_y), dx(_dx), dy(_dy), xo(_xo), yo(_yo) {}
    double x, y, dx, dy, xo, yo;
  };

  std::vector<TurtleState> stateStack;
  RAPT::rsRotationXY<double> rotLeft, rotRight; // use one, have functions apply, applyInverse
                                                // save sin/cos computation when angle changes
};

}

#endif