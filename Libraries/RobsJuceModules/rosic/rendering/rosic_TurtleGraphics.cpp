void TurtleGraphics::setAngle(double degrees)
{
  angle = degrees;
  updateRotations();
}

void TurtleGraphics::setAngleDelta(double degrees)
{
  angleDelta = degrees;
  updateRotations();
}

void TurtleGraphics::init(double _x, double _y, double _dx, double _dy, bool keepOldXY)
{
  if(keepOldXY) { xo =  x; yo = y; } // either keep current position in old position (xo,yo) (to
  else          { xo = _x; _y = 0; } // draw a connection) or reset them to given values as well

  x  = _x;
  y  = _y;
  dx = _dx;
  dy = _dy;

  stateStack.clear();
}

void TurtleGraphics::goForward()
{
  xo = x;
  yo = y;
  if(!reverseMode) { // maybe use a forwardMode flag to avoid negations in "normal" mode (optimization)
    x += dx;
    y += dy; }
  else {
    x -= dx;
    y -= dy; }
}

void TurtleGraphics::goBackward()
{
  xo = x;
  yo = y;
  if(!reverseMode) {
    x -= dx;
    y -= dy; }
  else {
    x += dx;
    y += dy; }
}

void TurtleGraphics::turnLeft()
{
  if(!reverseMode)
    rotLeft.apply(&dx, &dy);
  else
    rotRight.apply(&dx, &dy);
}

void TurtleGraphics::turnRight()
{
  if(!reverseMode)
    rotRight.apply(&dx, &dy); // use rotation.applyInverse
  else
    rotLeft.apply(&dx, &dy);
}

void TurtleGraphics::turnAround()
{
  dx = -dx;
  dy = -dy;
}

void TurtleGraphics::pushState()
{
  stateStack.push_back(TurtleState(x, y, dx, dy, xo, yo));

  // i think, in reverse mode, we should swap push/pop operations, too
}

void TurtleGraphics::popState()
{
  //rsAssert(stateStack.size() > 0); // may happen during editing the rules
  if(stateStack.size() == 0)
    return;

  // factor out into setState(TurtleState& s):
  TurtleState s = RAPT::rsGetAndRemoveLast(stateStack);
  x  = s.x;   y = s.y;
  dx = s.dx; dy = s.dy;
  xo = s.xo; yo = s.yo;
}

void TurtleGraphics::translate(const std::string& str, 
  std::vector<double>& vx, std::vector<double>& vy)
{
  // clear and add initial vertex:
  vx.clear();
  vy.clear();
  vx.push_back(getX()); 
  vy.push_back(getY()); 

  // loop through the string and add vertices as needed:
  for(size_t i = 0; i < str.size(); i++) {
    if(interpretCharacter(str[i])) {
      vx.push_back(getX()); 
      vy.push_back(getY()); 
    }
  }
  
  /*
  // try, if it is faster to pre-compute and pre-allocate the number of lines (given that the 
  // vectors vx, vy already have enough capacity - i.e. they have reserved enough memory via 
  // reserve() but the size() is still lower than the reserved number):
  int numLines = getNumberOfLines(str);
  vx.resize(numLines+1);
  vy.resize(numLines+1);
  vx[0] = getX();
  vy[0] = getY();
  int j = 1;
  for(int i = 0; i < str.size(); i++) {
    if(interpretCharacter(str[i])) {
      vx[j] = getX(); 
      vy[j] = getY();
      j++;
    }
  }
  // hmm...doesn't seem to make much difference
  */
}

bool TurtleGraphics::interpretCharacter(char c)
{
  if(c == 'F') { goForward();  return true;  }
  if(c == '+') { turnLeft();   return false; }
  if(c == '-') { turnRight();  return false; }
  if(c == '|') { turnAround(); return false; }
  if(c == 'f') { goForward();  return false; }
  if(c == 'G') { goForward();  return true;  }
  if(c == '[') { pushState();  return false; }
  if(c == ']') { popState();   return false; }
  return false;

  // write performance test and check, if a switch statement is faster...or maybe make a table
  // that maps characters to a member function pointer (non-command chars map to a dummy function)

  // maybe use commands like +30F-45F to turn 30° left, go forward, turn 45° right, go forward
  // ...so whenever a number comes after a + or -, don't use the default angle that is set up but 
  // the angle given/ by that number - maybe don't call standard turnLeft but a turnLeftBy(angle)
  // function, same for right - but the same effect can be also achieved by using 15 degrees
  // and doing ++F---F, so maybe it doens't really formally extend the possibilities - but some
  // things could be expressed more conveniently...one could use a pentagonal initiator with 
  // triangular generators - expressing this with +++++..., would require to set the angle to 2°
  // and using an excessive number of +ses
}

std::string TurtleGraphics::extractCommands(const std::string& s)
{
  std::string tmp;
  for(size_t i = 0; i < s.size(); i++)
    if(isCommand(s[i])) tmp += s[i];
  return tmp;
}

bool TurtleGraphics::isCommand(char c)
{
  return c=='+' || c=='-' || c=='F' || c=='f' || c=='G' || c=='[' || c==']' || c=='|';
}

bool TurtleGraphics::isLineCommand(char c)
{
  return c == 'F' || c == 'G';
}

int TurtleGraphics::getNumberOfLines(const std::string& s)
{
  int n = 0;
  for(size_t i = 0; i < s.size(); i++)
    if(isLineCommand(s[i])) n++;
    //if(s[i] == 'F' || s[i] == 'G') n++;
  return n;
}

void TurtleGraphics::updateRotations()
{
  double a = RAPT::rsDegreeToRadiant(angle);
  double d = RAPT::rsDegreeToRadiant(angleDelta);
  rotLeft.setAngle(  a + d);
  rotRight.setAngle(-a + d); 
  // entries of the rotation matrices that should be 0 are at 6.123e-17 - probably not an issue in                     
  // practice but not nice - maybe clamp to zero anything below 1.e-16
}