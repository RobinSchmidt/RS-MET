#include "LindenmayerSystem.h"

void LindenmayerSystem::addRule(char input, const std::string& output)
{
  ruleInputs.push_back(input);
  ruleOutputs.push_back(output);
}

void LindenmayerSystem::clearRules()
{
  ruleInputs.clear();
  ruleOutputs.clear();
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

void TurtleGraphics::setAngle(double degrees)
{
  double radians = RAPT::rsDegreeToRadiant(degrees);
  //long double factor =  (PI / 180.0);
  //double radians = (double) (factor * (long double) degrees);
  rotLeft.setAngle(  radians);
  rotRight.setAngle(-radians);
  // entries of the rotation matrices that should be 0 are at 6.123e-17 - probably not an issue in 
  // practice but not nice - maybe clamp to zero anything below 1.e-16
}

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
  rotLeft.apply(&dx, &dy);
  //double tx = dx;
  //dx = -dy;
  //dy =  tx;
}

void TurtleGraphics::turnRight()
{
  rotRight.apply(&dx, &dy);
  //double tx = dx;
  //dx =  dy;
  //dy = -tx;
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
    if(str[i] == 'f')    // this is nonstandard, standard would be to go forward but avoid
      goForward();       // drawing a line (that cannot expressed here)
    if(str[i] == 'F') { 
      goForward(); 
      vx.push_back(x); 
      vy.push_back(y); 
    }
  }
}

//=================================================================================================

void LindenmayerRenderer::getKochSnowflake(int N, std::vector<double>& x, std::vector<double>& y)
{
  clearRules();
  addRule('F', "F+F--F+F");
  std::string result = apply("F--F--F", N);

  turtleGraphics.init(0, 0, 1, 0);
  turtleGraphics.setAngle(60);
  turtleGraphics.translate(result, x, y);

  // ...normalize
}

void LindenmayerRenderer::getMooreCurve(int N, std::vector<double>& x, std::vector<double>& y)
{

}

void LindenmayerRenderer::normalizeXY(std::vector<double>& x, std::vector<double>& y)
{

}