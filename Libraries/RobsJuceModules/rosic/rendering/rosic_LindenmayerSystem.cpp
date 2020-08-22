//#include "LindenmayerSystem.h"

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
  for(size_t i = 0; i < s.length(); i++)
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

void LindenmayerRenderer::getKochSnowflake(int N, std::vector<double>& x, std::vector<double>& y)
{
  clearRules();
  addRule('F', "F+F--F+F");
  render("F--F--F", N, 60, x, y);
}

void LindenmayerRenderer::getMooreCurve(int N, std::vector<double>& x, std::vector<double>& y)
{
  clearRules();
  addRule('L', "-RF+LFL+FR-");
  addRule('R', "+LF-RFR-FL+");
  render("LFL+F+LFL", N, 90, x, y);
}

void LindenmayerRenderer::get32SegmentCurve(int N, std::vector<double>& x, std::vector<double>& y)
{
  clearRules();
  addRule('F', "-F+F-F-F+F+FF-F+F+FF+F-F-FF+FF-FF+F+F-FF-F-F+FF-F-F+F+F-F+");
  render("F+F+F+F", N, 90, x, y);
}

void LindenmayerRenderer::getQuadraticKochIsland(int N, 
  std::vector<double>& x, std::vector<double>& y)
{
  clearRules();
  addRule('F', "F-F+F+FFF-F-F+F");
  render("F+F+F+F", N, 90, x, y);
}

void LindenmayerRenderer::getSquareCurve(int N, std::vector<double>& x, std::vector<double>& y)
{
  clearRules();
  addRule('X', "XF-F+F-XF+F+XF-F+F-X");
  render("F+XF+F+XF", N, 90, x, y);
}

void LindenmayerRenderer::getSierpinskiTriangle(int N, 
  std::vector<double>& x, std::vector<double>& y)
{
  clearRules();
  addRule('F', "FF");
  addRule('X', "--FXF++FXF++FXF--");
  render("FXF--FF--FF", N, 60, x, y);
}

void LindenmayerRenderer::getSierpinskiTriangle2(int N, 
  std::vector<double>& x, std::vector<double>& y)
{
  clearRules();
  addRule('F', "F-G+F+G-F");
  addRule('G', "GG");
  render("F-G-G", N, 120, x, y);
}

void LindenmayerRenderer::getPleasantError(int N, std::vector<double>& x, std::vector<double>& y)
{
  clearRules();
  addRule('F', "F-F++F+F-F-F");         // original
  //addRule('F', "F-F++F+F-F-F+F-F");     // variation 1
  //addRule('F', "F-F++F+F-F-F-F+F");
  //addRule('F', "F-F++F+F-F-F+F+F-F-F");
  //addRule('F', "F-F++F+F-F-F-F-F+F+F");
  render("F-F-F-F-F", N, 72, x, y);
  // hmm...it seems, i need to code a module where the user can enter L-system rules, axioms, etc.
  // i cannot hardcode all interesting curves - maybe it should be possible to detune left/right
}

void LindenmayerRenderer::render(const std::string& seed, int order, double angle, 
  std::vector<double>& x, std::vector<double>& y)
{
  turtleGraphics.setAngle(angle);
  render(seed, order, x, y);
}

void LindenmayerRenderer::translate(const std::string& str, double angle,
  std::vector<double>& x, std::vector<double>& y)
{
  turtleGraphics.setAngle(angle);
  translate(str, x, y);
}

void LindenmayerRenderer::render(const std::string& seed, int order, 
  std::vector<double>& x, std::vector<double>& y)
{
  std::string str = apply(seed, order);
  translate(str, x, y);
}

void LindenmayerRenderer::translate(const std::string& str, 
  std::vector<double>& x, std::vector<double>& y)
{
  turtleGraphics.init(0, 0, 1, 0);
  turtleGraphics.translate(str, x, y);
  if(normalize) 
    normalizeXY(x, y);
}

void LindenmayerRenderer::normalizeXY(std::vector<double>& x, std::vector<double>& y)
{
  int N = (int)x.size();

  // remove mean:
  bool loop = true;  // make member
  if(loop == true) { // we need to ignore last element in mean computation
    double m;        // because it just repeats the first and is irrelevant
    m = RAPT::rsArrayTools::mean(&x[0], N-1); RAPT::rsArrayTools::add(&x[0], -m, &x[0], N);
    m = RAPT::rsArrayTools::mean(&y[0], N-1); RAPT::rsArrayTools::add(&y[0], -m, &y[0], N); }
  else {
    RAPT::rsArrayTools::removeMean(&x[0], N);
    RAPT::rsArrayTools::removeMean(&y[0], N); }

  // adjust maximum:
  double maxX = RAPT::rsArrayTools::maxAbs(&x[0], N);
  double maxY = RAPT::rsArrayTools::maxAbs(&y[0], N);
  double scl = 1.0 / max(maxX, maxY);
  RAPT::rsArrayTools::scale(&x[0], &x[0], N, scl);
  RAPT::rsArrayTools::scale(&y[0], &y[0], N, scl);
}

// other closed curves that can be generated:
// http://mathforum.org/advanced/robertd/lsys2d.html (many curves with L-system rules)
// http://www.kevs3d.co.uk/dev/lsystems/ (applet with examples)

// https://www.cut-the-knot.org/do_you_know/hilbert.shtml

// not sure, if doable by L-system:
// https://en.wikipedia.org/wiki/Sierpi%C5%84ski_curve