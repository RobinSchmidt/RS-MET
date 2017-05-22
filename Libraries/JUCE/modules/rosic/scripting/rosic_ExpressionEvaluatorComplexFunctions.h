#ifndef rosic_ExpressionEvaluatorComplexFunctions_h
#define rosic_ExpressionEvaluatorComplexFunctions_h

/**

This file defines custom functions (and mathematical constants) for the ExpressionEvaluator for
complex arithmetic. The arguments of the functions represent real and imaginary part of the
operands/input values, the reference arguments are used for the return values.

A unary function such as expC (complex exponential) works like this:
expC(z.re, z.im, &result.re, &result.im)

A binary function, such as addC (complex addition) works like this:
addC(z.re, z.im, w.re, w.im, &result.re, &result.im)


*/

// standard library includes:
//#include <errno.h>
//extern int errno;  // is set in functions from math.h and can be checked here

// third party includes:
#include "../_third_party/ExprEval_v3_4/expreval.h"
//using namespace ExprEval;
//using namespace std;

//// rosic-indcludes:
//#include "../math/rosic_ComplexFunctions.h"
//using namespace rosic;

namespace ExprEval
{

//-------------------------------------------------------------------------------------------------
// absolute value of a complex number:

class FunctionNodeAbsC : public FunctionNode
{
public:
  FunctionNodeAbsC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 0, 0, 0, 0);
  }

  double DoEvaluate()
  {
    double re = m_nodes[0]->Evaluate();
    double im = m_nodes[1]->Evaluate();
    return sqrt(re*re + im*im);
  }
};

class FunctionFactoryAbsC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "absC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeAbsC(expression); }
};

//-------------------------------------------------------------------------------------------------
// complex arcus cosine function:

class FunctionNodeAcosC : public FunctionNode
{
public:
  FunctionNodeAcosC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 2, 2, 0, 0);
  }

  double DoEvaluate()
  {
    rosic::Complex z  = rosic::Complex(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
    rosic::Complex w  = acosC(z);
    *m_refs[0] = w.re;
    *m_refs[1] = w.im;
    return 0.0;
  }
};

class FunctionFactoryAcosC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "acosC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeAcosC(expression); }
};

//-------------------------------------------------------------------------------------------------
// complex arcus hyperbolic cosine function:

class FunctionNodeAcoshC : public FunctionNode
{
public:
  FunctionNodeAcoshC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 2, 2, 0, 0);
  }

  double DoEvaluate()
  {
    rosic::Complex z  = rosic::Complex(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
    rosic::Complex w  = acoshC(z);
    *m_refs[0] = w.re;
    *m_refs[1] = w.im;
    return 0.0;
  }
};

class FunctionFactoryAcoshC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "acoshC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeAcoshC(expression); }
};

//-------------------------------------------------------------------------------------------------
// add two complex numbers:

class FunctionNodeAddC : public FunctionNode
{
public:
  FunctionNodeAddC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(4, 4, 2, 2, 0, 0);
  }

  double DoEvaluate()
  {
    *m_refs[0] = m_nodes[0]->Evaluate() + m_nodes[2]->Evaluate();
    *m_refs[1] = m_nodes[1]->Evaluate() + m_nodes[3]->Evaluate();
    return 0.0;
  }
};

class FunctionFactoryAddC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "addC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeAddC(expression); }
};

//-------------------------------------------------------------------------------------------------
// angle of a complex number:

class FunctionNodeAngleC : public FunctionNode
{
public:
  FunctionNodeAngleC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 0, 0, 0, 0);
  }

  double DoEvaluate()
  {
    double re = m_nodes[0]->Evaluate();
    double im = m_nodes[1]->Evaluate();
    if((re==0.0) && (im==0))
      return 0.0;
    else
      return atan2(im, re);
  }
};

class FunctionFactoryAngleC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "angleC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeAngleC(expression); }
};

//-------------------------------------------------------------------------------------------------
// complex arcus sine function:

class FunctionNodeAsinC : public FunctionNode
{
public:
  FunctionNodeAsinC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 2, 2, 0, 0);
  }

  double DoEvaluate()
  {
    rosic::Complex z  = rosic::Complex(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
    rosic::Complex w  = asinC(z);
    *m_refs[0] = w.re;
    *m_refs[1] = w.im;
    return 0.0;
  }
};

class FunctionFactoryAsinC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "asinC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeAsinC(expression); }
};

//-------------------------------------------------------------------------------------------------
// complex arcus hyperbolic sine function:

class FunctionNodeAsinhC : public FunctionNode
{
public:
  FunctionNodeAsinhC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 2, 2, 0, 0);
  }

  double DoEvaluate()
  {
    rosic::Complex z  = rosic::Complex(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
    rosic::Complex w  = asinhC(z);
    *m_refs[0] = w.re;
    *m_refs[1] = w.im;
    return 0.0;
  }
};

class FunctionFactoryAsinhC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "asinhC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeAsinhC(expression); }
};

//-------------------------------------------------------------------------------------------------
// complex arcus tangent function:

class FunctionNodeAtanC : public FunctionNode
{
public:
  FunctionNodeAtanC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 2, 2, 0, 0);
  }

  double DoEvaluate()
  {
    rosic::Complex z  = rosic::Complex(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
    rosic::Complex w  = atanC(z);
    *m_refs[0] = w.re;
    *m_refs[1] = w.im;
    return 0.0;
  }
};

class FunctionFactoryAtanC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "atanC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeAtanC(expression); }
};

//-------------------------------------------------------------------------------------------------
// complex arcus hyperbolic tangent function:

class FunctionNodeAtanhC : public FunctionNode
{
public:
  FunctionNodeAtanhC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 2, 2, 0, 0);
  }

  double DoEvaluate()
  {
    rosic::Complex z  = rosic::Complex(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
    rosic::Complex w  = atanhC(z);
    *m_refs[0] = w.re;
    *m_refs[1] = w.im;
    return 0.0;
  }
};

class FunctionFactoryAtanhC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "atanhC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeAtanhC(expression); }
};

//-------------------------------------------------------------------------------------------------
// complex cosine function:

class FunctionNodeCosC : public FunctionNode
{
public:
  FunctionNodeCosC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 2, 2, 0, 0);
  }

  double DoEvaluate()
  {
    rosic::Complex z  = rosic::Complex(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
    rosic::Complex w  = cosC(z);
    *m_refs[0] = w.re;
    *m_refs[1] = w.im;
    return 0.0;
  }
};

class FunctionFactoryCosC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "cosC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeCosC(expression); }
};

//-------------------------------------------------------------------------------------------------
// complex hypebolic cosine function:

class FunctionNodeCoshC : public FunctionNode
{
public:
  FunctionNodeCoshC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 2, 2, 0, 0);
  }

  double DoEvaluate()
  {
    rosic::Complex z  = rosic::Complex(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
    rosic::Complex w  = coshC(z);
    *m_refs[0] = w.re;
    *m_refs[1] = w.im;
    return 0.0;
  }
};

class FunctionFactoryCoshC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "coshC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeCoshC(expression); }
};

//-------------------------------------------------------------------------------------------------
// divides one complex number by another complex number:

class FunctionNodeDivC : public FunctionNode
{
public:
  FunctionNodeDivC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(4, 4, 2, 2, 0, 0);
  }

  double DoEvaluate()
  {
    double zRe   = m_nodes[0]->Evaluate();
    double zIm   = m_nodes[1]->Evaluate();
    double wRe   = m_nodes[2]->Evaluate();
    double wIm   = m_nodes[3]->Evaluate();
    double scale = 1.0 / (wRe*wRe + wIm*wIm);
    *m_refs[0]   = scale*(zRe*wRe + zIm*wIm);
    *m_refs[1]   = scale*(zIm*wRe - zRe*wIm);
    return 0.0;
  }
};

class FunctionFactoryDivC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "divC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeDivC(expression); }
};

//-------------------------------------------------------------------------------------------------
// complex exponential funnction (todo: optimize):

class FunctionNodeExpC : public FunctionNode
{
public:
  FunctionNodeExpC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 2, 2, 0, 0);
  }

  double DoEvaluate()
  {
    double tmp = exp(m_nodes[0]->Evaluate());
    rosic::sinCos(m_nodes[1]->Evaluate(), m_refs[1], m_refs[0]);
    *m_refs[0] *= tmp;
    *m_refs[1] *= tmp;
    return 0.0;
  }
};

class FunctionFactoryExpC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "expC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeExpC(expression); }
};

//-------------------------------------------------------------------------------------------------
// complex logarithm function:

class FunctionNodeLogC : public FunctionNode
{
public:
  FunctionNodeLogC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 2, 2, 0, 0);
  }

  double DoEvaluate()
  {
    rosic::Complex z  = rosic::Complex(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
    rosic::Complex w  = logC(z);
    *m_refs[0] = w.re;
    *m_refs[1] = w.im;
    return 0.0;
  }
};

class FunctionFactoryLogC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "logC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeLogC(expression); }
};

//-------------------------------------------------------------------------------------------------
// multiplies two complex numbers:

class FunctionNodeMulC : public FunctionNode
{
public:
  FunctionNodeMulC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(4, 4, 2, 2, 0, 0);
  }

  double DoEvaluate()
  {
    double zRe = m_nodes[0]->Evaluate();
    double zIm = m_nodes[1]->Evaluate();
    double wRe = m_nodes[2]->Evaluate();
    double wIm = m_nodes[3]->Evaluate();
    *m_refs[0] = zRe*wRe-zIm*wIm;
    *m_refs[1] = zRe*wIm+zIm*wRe;
    return 0.0;
  }
};

class FunctionFactoryMulC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "mulC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeMulC(expression); }
};


//-------------------------------------------------------------------------------------------------
// complex power function:

class FunctionNodePowC : public FunctionNode
{
public:
  FunctionNodePowC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(4, 4, 2, 2, 0, 0);
  }

  double DoEvaluate()
  {
    rosic::Complex basis    = rosic::Complex(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
    rosic::Complex exponent = rosic::Complex(m_nodes[2]->Evaluate(), m_nodes[3]->Evaluate());
    rosic::Complex result   = powC(basis, exponent);
    *m_refs[0]       = result.re;
    *m_refs[1]       = result.im;
    return 0.0;
  }
};

class FunctionFactoryPowC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "powC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodePowC(expression); }
};

//-------------------------------------------------------------------------------------------------
// complex sine function:

class FunctionNodeSinC : public FunctionNode
{
public:
  FunctionNodeSinC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 2, 2, 0, 0);
  }

  double DoEvaluate()
  {
    rosic::Complex z  = rosic::Complex(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
    rosic::Complex w  = sinC(z);
    *m_refs[0] = w.re;
    *m_refs[1] = w.im;
    return 0.0;
  }
};

class FunctionFactorySinC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "sinC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeSinC(expression); }
};

//-------------------------------------------------------------------------------------------------
// complex hyperbolic sine function:

class FunctionNodeSinhC : public FunctionNode
{
public:
  FunctionNodeSinhC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 2, 2, 0, 0);
  }

  double DoEvaluate()
  {
    rosic::Complex z  = rosic::Complex(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
    rosic::Complex w  = sinhC(z);
    *m_refs[0] = w.re;
    *m_refs[1] = w.im;
    return 0.0;
  }
};

class FunctionFactorySinhC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "sinhC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeSinhC(expression); }
};

//-------------------------------------------------------------------------------------------------
// complex square root function:

class FunctionNodeSqrtC : public FunctionNode
{
public:
  FunctionNodeSqrtC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 2, 2, 0, 0);
  }

  double DoEvaluate()
  {
    rosic::Complex z  = rosic::Complex(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
    rosic::Complex w  = sqrtC(z);
    *m_refs[0] = w.re;
    *m_refs[1] = w.im;
    return 0.0;
  }
};

class FunctionFactorySqrtC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "sqrtC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeSqrtC(expression); }
};

//-------------------------------------------------------------------------------------------------
// subtract two complex numbers:

class FunctionNodeSubC : public FunctionNode
{
public:
  FunctionNodeSubC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(4, 4, 2, 2, 0, 0);
  }

  double DoEvaluate()
  {
    *m_refs[0] = m_nodes[0]->Evaluate() - m_nodes[2]->Evaluate();
    *m_refs[1] = m_nodes[1]->Evaluate() - m_nodes[3]->Evaluate();
    return 0.0;
  }
};

class FunctionFactorySubC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "subC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeSubC(expression); }
};

//-------------------------------------------------------------------------------------------------
// complex tangent function:

class FunctionNodeTanC : public FunctionNode
{
public:
  FunctionNodeTanC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 2, 2, 0, 0);
  }

  double DoEvaluate()
  {
    rosic::Complex z  = rosic::Complex(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
    rosic::Complex w  = tanC(z);
    *m_refs[0] = w.re;
    *m_refs[1] = w.im;
    return 0.0;
  }
};

class FunctionFactoryTanC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "tanC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeTanC(expression); }
};

//-------------------------------------------------------------------------------------------------
// complex hyperbolic tangent function:

class FunctionNodeTanhC : public FunctionNode
{
public:
  FunctionNodeTanhC(ExprEval::Expression *expression) : FunctionNode(expression)
  {
    SetArgumentCount(2, 2, 2, 2, 0, 0);
  }

  double DoEvaluate()
  {
    rosic::Complex z  = rosic::Complex(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
    rosic::Complex w  = tanhC(z);
    *m_refs[0] = w.re;
    *m_refs[1] = w.im;
    return 0.0;
  }
};

class FunctionFactoryTanhC : public FunctionFactory
{
public:
  std::string GetName() const                    { return "tanhC"; }
  FunctionNode* DoCreate(ExprEval::Expression *expression) { return new FunctionNodeTanhC(expression); }
};

} // namespace ExprEval

#endif  // rosic_ExpressionEvaluatorComplexFunctions_h
