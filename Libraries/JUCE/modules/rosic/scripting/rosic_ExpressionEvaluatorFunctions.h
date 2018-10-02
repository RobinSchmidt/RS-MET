#ifndef rosic_ExpressionEvaluatorFunctions_h
#define rosic_ExpressionEvaluatorFunctions_h

/**

This file defines custom functions (and mathematical constants) for the ExpressionEvaluator.  

*/

//// standard library includes:
//#include <errno.h>
////extern int errno;  // is set in functions from math.h and can be checked here

// third party includes:
//#include "../_third_party/ExprEval_v3_4/expreval.h"
//using namespace ExprEval;
//using namespace std;

//// rosic-indcludes:
//#include "../math/rosic_ElementaryFunctionsReal.h"
//using namespace rosic;

namespace ExprEval
{

//-------------------------------------------------------------------------------------------------
// Gaussian (aka normal) distribution:

class FunctionNodeGauss : public FunctionNode
{
public:
  FunctionNodeGauss(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(1, 3, 0, 0, 0, 0);
  }

  double DoEvaluate()
  {
    if(m_nodes.size() < 2)
      return rosic::gauss(m_nodes[0]->Evaluate(), 0.0, 1.0);
    else if(m_nodes.size() < 3)
      return rosic::gauss(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate(), 1.0);
    else
      return rosic::gauss(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate(), m_nodes[2]->Evaluate());
  }
};

class FunctionFactoryGauss : public FunctionFactory
{
public:
  std::string GetName() const                    { return "gauss"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeGauss(expression); }
};

//-------------------------------------------------------------------------------------------------
// Chebychev polynomial of arbitrary order:

class FunctionNodeCheby : public FunctionNode
{
public:
  FunctionNodeCheby(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 0, 0, 0, 0);
  }

  double DoEvaluate()
  {
    return rosic::cheby(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
  }
};

class FunctionFactoryCheby : public FunctionFactory
{
public:
  std::string GetName() const                    { return "cheby"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeCheby(expression); }
};

//-------------------------------------------------------------------------------------------------
// power function:

class FunctionNodePow : public FunctionNode
{
public:
  FunctionNodePow(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 0, 0, 0, 0);
  }

  double DoEvaluate()
  {
    return pow(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
  }
};

class FunctionFactoryPow : public FunctionFactory
{
public:
  std::string GetName() const                    { return "pow"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodePow(expression); }
};

//-------------------------------------------------------------------------------------------------
// bipolar power function:

class FunctionNodePowBipolar : public FunctionNode
{
public:
  FunctionNodePowBipolar(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 0, 0, 0, 0);
  }

  double DoEvaluate()
  {
    return rosic::powBipolar(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
  }
};

class FunctionFactoryPowBipolar : public FunctionFactory
{
public:
  std::string GetName() const                    { return "powBipolar"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodePowBipolar(expression); }
};

//-------------------------------------------------------------------------------------------------
// quantization to an arbitrary interval:

class FunctionNodeQuant : public FunctionNode
{
public:
  FunctionNodeQuant(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 0, 0, 0, 0);
  }

  double DoEvaluate()
  {
    return rosic::quant(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
  }
};

class FunctionFactoryQuant : public FunctionFactory
{
public:
  std::string GetName() const                    { return "quant"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeQuant(expression); }
};

//-------------------------------------------------------------------------------------------------
// quantization to some number of bits:

class FunctionNodeQuantToBits : public FunctionNode
{
public:
  FunctionNodeQuantToBits(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 0, 0, 0, 0);
  }

  double DoEvaluate()
  {
    return rosic::quantToBits(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
  }
};

class FunctionFactoryQuantToBits : public FunctionFactory
{
public:
  std::string GetName() const                    { return "quantToBits"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeQuantToBits(expression); }
};

//-------------------------------------------------------------------------------------------------
// sawtooth function:

class FunctionNodeSaw : public FunctionNode
{

public:

  FunctionNodeSaw(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(
      1,               // minimum number of normal arguments
      1,               // maximum number of normal arguments
      0,               // minimum number of reference arguments
      0,               // maximum number of reference arguments
      0,               // minimum number of data arguments
      0                // maximum number of data arguments
      );
  }

  double DoEvaluate()
  {
    return rosic::sawWave(m_nodes[0]->Evaluate());
  }

};

class FunctionFactorySaw : public FunctionFactory
{
public:
  std::string GetName() const                    { return "saw"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeSaw(expression); }
};


//-------------------------------------------------------------------------------------------------
// signum function:

class FunctionNodeSign : public FunctionNode
{
public:
  FunctionNodeSign(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(1, 1, 0, 0, 0, 0);
  }

  double DoEvaluate()
  {
    return RAPT::rsSign(m_nodes[0]->Evaluate());
  }
};

class FunctionFactorySign : public FunctionFactory
{
public:
  std::string GetName() const                    { return "sign"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeSign(expression); }
};


//-------------------------------------------------------------------------------------------------
// soft-clipping:

class FunctionNodeSoftClip : public FunctionNode
{
public:
  FunctionNodeSoftClip(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(1, 5, 0, 0, 0, 0);
  }

  double DoEvaluate()
  {
    if(m_nodes.size() < 2)
      return rosic::softClip(m_nodes[0]->Evaluate());
    else if(m_nodes.size() < 3)
      return rosic::softClip(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
    else if(m_nodes.size() < 4)
      return rosic::softClip(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate(), m_nodes[2]->Evaluate());
    else if(m_nodes.size() < 5)
      return rosic::softClip(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate(), m_nodes[2]->Evaluate(),
        m_nodes[3]->Evaluate());
    else
      return rosic::softClip(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate(), m_nodes[2]->Evaluate(),
        m_nodes[3]->Evaluate(), m_nodes[4]->Evaluate());

  }
};

class FunctionFactorySoftClip : public FunctionFactory
{
public:
  std::string GetName() const                    { return "softClip"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeSoftClip(expression); }
};

//-------------------------------------------------------------------------------------------------
// step function:

class FunctionNodeStep : public FunctionNode
{
public:
  FunctionNodeStep(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(1, 1, 0, 0, 0, 0);
  }

  double DoEvaluate()
  {
    return rosic::step(m_nodes[0]->Evaluate());
  }
};

class FunctionFactoryStep : public FunctionFactory
{
public:
  std::string GetName() const                    { return "step"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeStep(expression); }
};

} // namespace ExprEval



#endif  // rosic_ExpressionEvaluatorFunctions_h