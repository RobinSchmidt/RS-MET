#include "LindenmayerSystem.h"

void LindenmayerSystem::addRule(char input, const std::string& output)
{
  ruleInputs.push_back(input);
  ruleOutputs.push_back(output);
}

std::string LindenmayerSystem::apply(char c)
{
  for(size_t i = 0; i < ruleInputs.size(); i++)
    if( ruleInputs[i] == c )
      return ruleOutputs[i];
  return std::string(1, c);   // 'c' repeated one time
}

std::string LindenmayerSystem::apply(const std::string& s)
{
  std::string r; // result
  for(int i = 0; i < s.length(); i++)
    r += apply(s[i]);
  return r;
}

std::string LindenmayerSystem::apply(const std::string& s, int num)
{
  std::string r = s;
  for(int i = 1; i <= num; i++)
    r = apply(r);
  return r; 
}

//=================================================================================================

void TurtleGraphics::init(double _x, double _y, double _dx, double _dy)
{
  x  = _x;
  y  = _y;
  dx = _dx;
  dy = _dy;
}

void TurtleGraphics::goForward()
{
  x += dx;
  y += dy;
}

void TurtleGraphics::turnLeft()
{
  // todo: use left/right rotation matrices with adjustable angle (currently we turn 90°)
  double tx = dx;
  dx = -dy;
  dy =  tx;
}

void TurtleGraphics::turnRight()
{
  double tx = dx;
  dx =  dy;
  dy = -tx;
}

void TurtleGraphics::translate(const std::string& str, 
  std::vector<double>& vx, std::vector<double>& vy)
{
  // clear and add initial vertex:
  vx.clear();
  vy.clear();
  vx.push_back(x); 
  vy.push_back(y); 

  // loop through the string and add vertices as needed:
  for(int i = 0; i < str.size(); i++) {
    if(str[i] == '+')   
      turnLeft();
    if(str[i] == '-')   
      turnRight();
    if(str[i] == 'F') { 
      goForward(); 
      vx.push_back(x); 
      vy.push_back(y); 
    }
  }
}