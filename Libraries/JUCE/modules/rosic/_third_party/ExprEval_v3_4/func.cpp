// File:    func.cpp
// Author:  Brian Allen Vanderburg II
// Purpose: ExprEval internal functions
//------------------------------------------------------------------------------

//#define M_PI 3.1415926535897932384626433832795 // added by Robin Schmidt

// Includes
#include <new>
#include <memory>

// define macro added by Robin Schmidt:
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <cerrno>
#include <ctime>

#include "defs.h"
#include "funclist.h"
#include "node.h"
#include "except.h"

using namespace std;
using namespace ExprEval;

// Anonymous namespace for items
namespace
    {
    // Absolute value
    //--------------------------------------------------------------------------
    class abs_FunctionNode : public FunctionNode
        {
        public:
            abs_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(1, 1, 0, 0);
                }

            double DoEvaluate()
                {
                return fabs(m_nodes[0]->Evaluate());
                }
        };

    class abs_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "abs";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new abs_FunctionNode(expr);
                }
        };

    // Modulus
    //--------------------------------------------------------------------------
    class mod_FunctionNode : public FunctionNode
        {
        public:
            mod_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(2, 2, 0, 0);
                }

            double DoEvaluate()
                {
                double result;

                // Check for math errors
                errno = 0;

                result = fmod(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());
                if(errno)
                    throw(MathException(GetName()));

                return result;
                }
        };

    class mod_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "mod";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new mod_FunctionNode(expr);
                }
        };

    // Integer part
    //--------------------------------------------------------------------------
    class ipart_FunctionNode : public FunctionNode
        {
        public:
            ipart_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(1, 1, 0, 0);
                }

            double DoEvaluate()
                {
                double result;

                modf(m_nodes[0]->Evaluate(), &result);

                return result;
                }
        };

    class ipart_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "ipart";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new ipart_FunctionNode(expr);
                }
        };

    // Fraction part
    //--------------------------------------------------------------------------
    class fpart_FunctionNode : public FunctionNode
        {
        public:
            fpart_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(1, 1, 0, 0);
                }

            double DoEvaluate()
                {
                double dummy;

                return modf(m_nodes[0]->Evaluate(), &dummy);
                }
        };

    class fpart_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "fpart";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new fpart_FunctionNode(expr);
                }
        };

    // Minimum
    //--------------------------------------------------------------------------
    class min_FunctionNode : public FunctionNode
        {
        public:
            min_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(2, -1, 0, 0);
                }

            double DoEvaluate()
                {
                vector<Node*>::size_type pos;

                double result = m_nodes[0]->Evaluate();

                for(pos = 1; pos < m_nodes.size(); pos++)
                    {
                    double tmp = m_nodes[pos]->Evaluate();
                    if(tmp < result)
                        result = tmp;
                    }

                return result;
                }
        };

    class min_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "min";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new min_FunctionNode(expr);
                }
        };

    // Maximum
    //--------------------------------------------------------------------------
    class max_FunctionNode : public FunctionNode
        {
        public:
            max_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(2, -1, 0, 0);
                }

            double DoEvaluate()
                {
                vector<Node*>::size_type pos;

                double result = m_nodes[0]->Evaluate();

                for(pos = 1; pos < m_nodes.size(); pos++)
                    {
                    double tmp = m_nodes[pos]->Evaluate();
                    if(tmp > result)
                        result = tmp;
                    }

                return result;
                }
        };

    class max_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "max";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new max_FunctionNode(expr);
                }
        };

    // Square root
    //--------------------------------------------------------------------------
    class sqrt_FunctionNode : public FunctionNode
        {
        public:
            sqrt_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(1, 1, 0, 0);
                }

            double DoEvaluate()
                {
                errno = 0;

                double result = sqrt(m_nodes[0]->Evaluate());

                if(errno)
                    throw(MathException(GetName()));

                return result;
                }
        };

    class sqrt_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "sqrt";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new sqrt_FunctionNode(expr);
                }
        };

    // Sine
    //--------------------------------------------------------------------------
    class sin_FunctionNode : public FunctionNode
        {
        public:
            sin_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(1, 1, 0, 0);
                }

            double DoEvaluate()
                {
                errno = 0;

                double result = sin(m_nodes[0]->Evaluate());

                if(errno)
                    throw(MathException(GetName()));

                return result;
                }
        };

    class sin_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "sin";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new sin_FunctionNode(expr);
                }
        };

    // Cosine
    //--------------------------------------------------------------------------
    class cos_FunctionNode : public FunctionNode
        {
        public:
            cos_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(1, 1, 0, 0);
                }

            double DoEvaluate()
                {
                errno = 0;

                double result = cos(m_nodes[0]->Evaluate());

                if(errno)
                    throw(MathException(GetName()));

                return result;
                }
        };

    class cos_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "cos";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new cos_FunctionNode(expr);
                }
        };

    // Tangent
    //--------------------------------------------------------------------------
    class tan_FunctionNode : public FunctionNode
        {
        public:
            tan_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(1, 1, 0, 0);
                }

            double DoEvaluate()
                {
                errno = 0;

                double result = tan(m_nodes[0]->Evaluate());

                if(errno)
                    throw(MathException(GetName()));

                return result;
                }
        };

    class tan_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "tan";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new tan_FunctionNode(expr);
                }
        };

    // Hyperbolic Sine
    //--------------------------------------------------------------------------
    class sinh_FunctionNode : public FunctionNode
        {
        public:
            sinh_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(1, 1, 0, 0);
                }

            double DoEvaluate()
                {
                errno = 0;

                double result = sinh(m_nodes[0]->Evaluate());

                if(errno)
                    throw(MathException(GetName()));

                return result;
                }
        };

    class sinh_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "sinh";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new sinh_FunctionNode(expr);
                }
        };

    // Hyperbolic Cosine
    //--------------------------------------------------------------------------
    class cosh_FunctionNode : public FunctionNode
        {
        public:
            cosh_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(1, 1, 0, 0);
                }

            double DoEvaluate()
                {
                errno = 0;

                double result = cosh(m_nodes[0]->Evaluate());

                if(errno)
                    throw(MathException(GetName()));

                return result;
                }
        };

    class cosh_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "cosh";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new cosh_FunctionNode(expr);
                }
        };

    // Hyperbolic Tangent
    //--------------------------------------------------------------------------
    class tanh_FunctionNode : public FunctionNode
        {
        public:
            tanh_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(1, 1, 0, 0);
                }

            double DoEvaluate()
                {
                errno = 0;

                double result = tanh(m_nodes[0]->Evaluate());

                if(errno)
                    throw(MathException(GetName()));

                return result;
                }
        };

    class tanh_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "tanh";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new tanh_FunctionNode(expr);
                }
        };

    // Arc Sine
    //--------------------------------------------------------------------------
    class asin_FunctionNode : public FunctionNode
        {
        public:
            asin_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(1, 1, 0, 0);
                }

            double DoEvaluate()
                {
                errno = 0;

                double result = asin(m_nodes[0]->Evaluate());

                if(errno)
                    throw(MathException(GetName()));

                return result;
                }
        };

    class asin_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "asin";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new asin_FunctionNode(expr);
                }
        };

    // Arc Cosine
    //--------------------------------------------------------------------------
    class acos_FunctionNode : public FunctionNode
        {
        public:
            acos_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(1, 1, 0, 0);
                }

            double DoEvaluate()
                {
                errno = 0;

                double result = acos(m_nodes[0]->Evaluate());

                if(errno)
                    throw(MathException(GetName()));

                return result;
                }
        };

    class acos_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "acos";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new acos_FunctionNode(expr);
                }
        };

    // Arc Tangent
    //--------------------------------------------------------------------------
    class atan_FunctionNode : public FunctionNode
        {
        public:
            atan_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(1, 1, 0, 0);
                }

            double DoEvaluate()
                {
                errno = 0;

                double result = atan(m_nodes[0]->Evaluate());

                if(errno)
                    throw(MathException(GetName()));

                return result;
                }
        };

    class atan_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "atan";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new atan_FunctionNode(expr);
                }
        };

    // Arc Tangent 2: atan2(y, x)
    //--------------------------------------------------------------------------
    class atan2_FunctionNode : public FunctionNode
        {
        public:
            atan2_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(2, 2, 0, 0);
                }

            double DoEvaluate()
                {
                errno = 0;

                double result = atan2(m_nodes[0]->Evaluate(), m_nodes[1]->Evaluate());

                if(errno)
                    throw(MathException(GetName()));

                return result;
                }
        };

    class atan2_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "atan2";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new atan2_FunctionNode(expr);
                }
        };

    // Log
    //--------------------------------------------------------------------------
    class log_FunctionNode : public FunctionNode
        {
        public:
            log_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(1, 1, 0, 0);
                }

            double DoEvaluate()
                {
                errno = 0;

                double result = log10(m_nodes[0]->Evaluate());

                if(errno)
                    throw(MathException(GetName()));

                return result;
                }
        };

    class log_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "log";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new log_FunctionNode(expr);
                }
        };

    // Ln
    //--------------------------------------------------------------------------
    class ln_FunctionNode : public FunctionNode
        {
        public:
            ln_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(1, 1, 0, 0);
                }

            double DoEvaluate()
                {
                errno = 0;

                double result = log(m_nodes[0]->Evaluate());

                if(errno)
                    throw(MathException(GetName()));

                return result;
                }
        };

    class ln_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "ln";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new ln_FunctionNode(expr);
                }
        };

    // Exp
    //--------------------------------------------------------------------------
    class exp_FunctionNode : public FunctionNode
        {
        public:
            exp_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(1, 1, 0, 0);
                }

            double DoEvaluate()
                {
                errno = 0;

                double result = exp(m_nodes[0]->Evaluate());

                if(errno)
                    throw(MathException(GetName()));

                return result;
                }
        };

    class exp_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "exp";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new exp_FunctionNode(expr);
                }
        };

    // Logn
    //--------------------------------------------------------------------------
    class logn_FunctionNode : public FunctionNode
        {
        public:
            logn_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(2, 2, 0, 0);
                }

            double DoEvaluate()
                {
                errno = 0;

                // Check for division by zero
                double tmp = log(m_nodes[1]->Evaluate());

                if(tmp == 0.0)
                    throw(MathException(GetName()));

                // Calculate result
                double result = log(m_nodes[0]->Evaluate()) / tmp;

                if(errno)
                    throw(MathException(GetName()));

                return result;
                }
        };

    class logn_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "logn";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new logn_FunctionNode(expr);
                }
        };

    // Ceil
    //--------------------------------------------------------------------------
    class ceil_FunctionNode : public FunctionNode
        {
        public:
            ceil_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(1, 1, 0, 0);
                }

            double DoEvaluate()
                {
                return ceil(m_nodes[0]->Evaluate());
                }
        };

    class ceil_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "ceil";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new ceil_FunctionNode(expr);
                }
        };

    // Floor
    //--------------------------------------------------------------------------
    class floor_FunctionNode : public FunctionNode
        {
        public:
            floor_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(1, 1, 0, 0);
                }

            double DoEvaluate()
                {
                return floor(m_nodes[0]->Evaluate());
                }
        };

    class floor_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "floor";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new floor_FunctionNode(expr);
                }
        };

    // Rand
    //--------------------------------------------------------------------------

    // Get next random (0,1)
    inline double NextRandom(double *seed)
        {
        long a = (long)(*seed) * 214013L + 2531011L;
        *seed = (double)a;

        return (double)((a >> 16) & 0x7FFF) / 32767.0;
        }

    class rand_FunctionNode : public FunctionNode
        {
        public:
            rand_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(0, 0, 1, 1);
                }

            double DoEvaluate()
                {
                return NextRandom(m_refs[0]);
                }
        };

    class rand_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "rand";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new rand_FunctionNode(expr);
                }
        };

    // Random
    //--------------------------------------------------------------------------
    class random_FunctionNode : public FunctionNode
        {
        public:
            random_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(2, 2, 1, 1);
                }

            double DoEvaluate()
                {
                double a = m_nodes[0]->Evaluate();
                double b = m_nodes[1]->Evaluate();

                return NextRandom(m_refs[0]) * (b - a) + a;
                }
        };

    class random_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "random";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new random_FunctionNode(expr);
                }
        };

    // Randomize
    //--------------------------------------------------------------------------
    class randomize_FunctionNode : public FunctionNode
        {
        public:
            randomize_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(0, 0, 1, 1);
                }

            double DoEvaluate()
                {
                static long curcall = 1;

                *m_refs[0] = (double)(((clock() + 1024 + curcall) * time(NULL)) % 2176971487);
                curcall++;

                return 0.0;
                }
        };

    class randomize_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "randomize";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new randomize_FunctionNode(expr);
                }
        };

    // Radians to degrees
    //--------------------------------------------------------------------------
    class deg_FunctionNode : public FunctionNode
        {
        public:
            deg_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(1, 1, 0, 0);
                }

            double DoEvaluate()
                {
                return (m_nodes[0]->Evaluate() * 180.0) / 3.1415926535897932384626433832795;
                }
        };

    class deg_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "deg";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new deg_FunctionNode(expr);
                }
        };

    // Degrees to radians
    //--------------------------------------------------------------------------
    class rad_FunctionNode : public FunctionNode
        {
        public:
            rad_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(1, 1, 0, 0);
                }

            double DoEvaluate()
                {
                return (m_nodes[0]->Evaluate() * 3.1415926535897932384626433832795) / 180.0;
                }
        };

    class rad_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "rad";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new rad_FunctionNode(expr);
                }
        };

    // Rectangular to polar: rect2pol(x, y, &distance, &angle)
    //--------------------------------------------------------------------------
    class rect2pol_FunctionNode : public FunctionNode
        {
        public:
            rect2pol_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(2, 2, 2, 2);
                }

            double DoEvaluate()
                {
                errno = 0;

                double x = m_nodes[0]->Evaluate();
                double y = m_nodes[1]->Evaluate();

                double d = sqrt(x * x + y * y);
                double a = atan2(y, x);

                if(errno)
                    throw(MathException(GetName()));

                *m_refs[0] = d;
                if(a < 0.0)
                    *m_refs[1] = a + (2.0 * 3.1415926535897932384626433832795);
                else
                    *m_refs[1] = a;

                return d;
                }
        };

    class rect2pol_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "rect2pol";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new rect2pol_FunctionNode(expr);
                }
        };

    // Polar to rectangular: pol2rect(distance, angle, &x, &y)
    //--------------------------------------------------------------------------
    class pol2rect_FunctionNode : public FunctionNode
        {
        public:
            pol2rect_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(2, 2, 2, 2);
                }

            double DoEvaluate()
                {
                errno = 0;

                double d = m_nodes[0]->Evaluate();
                double a = m_nodes[1]->Evaluate();

                double x = d * cos(a);
                double y = d * sin(a);

                if(errno)
                    throw(MathException(GetName()));

                *m_refs[0] = x;
                *m_refs[1] = y;

                return x;
                }
        };

    class pol2rect_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "pol2rect";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new pol2rect_FunctionNode(expr);
                }
        };

    // If
    //--------------------------------------------------------------------------
    class if_FunctionNode : public FunctionNode
        {
        public:
            if_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(3, 3, 0, 0);
                }

            double DoEvaluate()
                {
                double c = m_nodes[0]->Evaluate();

                if(c == 0.0)
                    return m_nodes[2]->Evaluate();
                else
                    return m_nodes[1]->Evaluate();
                }
        };

    class if_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "if";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new if_FunctionNode(expr);
                }
        };

    // Select
    //--------------------------------------------------------------------------
    class select_FunctionNode : public FunctionNode
        {
        public:
            select_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(3, 4, 0, 0);
                }

            double DoEvaluate()
                {
                double c = m_nodes[0]->Evaluate();

                if(c < 0.0)
                    return m_nodes[1]->Evaluate();
                else if(c == 0.0)
                    return m_nodes[2]->Evaluate();
                else
                    {
                    if(m_nodes.size() == 3)
                        return m_nodes[2]->Evaluate();
                    else
                        return m_nodes[3]->Evaluate();
                    }
                }
        };

    class select_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "select";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new select_FunctionNode(expr);
                }
        };

    // Equal
    //--------------------------------------------------------------------------
    class equal_FunctionNode : public FunctionNode
        {
        public:
            equal_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(2, 2, 0, 0);
                }

            double DoEvaluate()
                {
                if(m_nodes[0]->Evaluate() == m_nodes[1]->Evaluate())
                    return 1.0;
                else
                    return 0.0;
                }
        };

    class equal_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "equal";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new equal_FunctionNode(expr);
                }
        };

    // Above
    //--------------------------------------------------------------------------
    class above_FunctionNode : public FunctionNode
        {
        public:
            above_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(2, 2, 0, 0);
                }

            double DoEvaluate()
                {
                if(m_nodes[0]->Evaluate() > m_nodes[1]->Evaluate())
                    return 1.0;
                else
                    return 0.0;
                }
        };

    class above_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "above";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new above_FunctionNode(expr);
                }
        };

    // Below
    //--------------------------------------------------------------------------
    class below_FunctionNode : public FunctionNode
        {
        public:
            below_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(2, 2, 0, 0);
                }

            double DoEvaluate()
                {
                if(m_nodes[0]->Evaluate() < m_nodes[1]->Evaluate())
                    return 1.0;
                else
                    return 0.0;
                }
        };

    class below_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "below";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new below_FunctionNode(expr);
                }
        };

    // Clip
    //--------------------------------------------------------------------------
    class clip_FunctionNode : public FunctionNode
        {
        public:
            clip_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(3, 3, 0, 0);
                }

            double DoEvaluate()
                {
                double v = m_nodes[0]->Evaluate();
                double a = m_nodes[1]->Evaluate();
                double b = m_nodes[2]->Evaluate();

                if(v < a)
                    return a;
                else if(v > b)
                    return b;
                else
                    return v;
                }
        };

    class clip_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "clip";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new clip_FunctionNode(expr);
                }
        };

    // Clamp
    //--------------------------------------------------------------------------
    class clamp_FunctionNode : public FunctionNode
        {
        public:
            clamp_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(3, 3, 0, 0);
                }

            double DoEvaluate()
                {
                double v = m_nodes[0]->Evaluate();
                double a = m_nodes[1]->Evaluate();
                double b = m_nodes[2]->Evaluate();

                if(a == b)
                    return a;
                else
                    {
                    double tmp = fmod(v - a, b - a);

                    if(tmp < 0)
                        return tmp + b;
                    else
                        return tmp + a;
                    }
                }
        };

    class clamp_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "clamp";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new clamp_FunctionNode(expr);
                }
        };

    // Rescale
    //--------------------------------------------------------------------------
    class rescale_FunctionNode : public FunctionNode
        {
        public:
            rescale_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(5, 5, 0, 0);
                }

            double DoEvaluate()
                {
                double pnt = m_nodes[0]->Evaluate();
                double o1 = m_nodes[1]->Evaluate();
                double o2 = m_nodes[2]->Evaluate();
                double n1 = m_nodes[3]->Evaluate();
                double n2 = m_nodes[4]->Evaluate();

                double odiff = o2 - o1;
                if(odiff == 0.0)
                    return n1;
                else
                    {
                    return (pnt - o1) * (n2 - n1) / odiff + n1;
                    }
                }
        };

    class rescale_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "rescale";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new rescale_FunctionNode(expr);
                }
        };

    // Poly: poly(x, c3, c2, c1, c0): c3*x^3 + c2*x^2 + c1*x + c0
    //--------------------------------------------------------------------------
    class poly_FunctionNode : public FunctionNode
        {
        public:
            poly_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(2, -1, 0, 0);
                }

            double DoEvaluate()
                {
                double total = 0.0;
                double curpow;

                vector<Node*>::size_type pos, count;
                count = m_nodes.size();

                curpow = (double)(count - 2);

                // Value of x
                double x = m_nodes[0]->Evaluate();

                errno = 0;

                for(pos = 1; pos < count; pos++)
                    {
                    total += (m_nodes[pos]->Evaluate() * pow(x, curpow));
                    curpow -= 1.0;
                    }

                if(errno)
                    throw(MathException(GetName()));

                return total;
                }
        };

    class poly_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "poly";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new poly_FunctionNode(expr);
                }
        };

    // And
    //--------------------------------------------------------------------------
    class and_FunctionNode : public FunctionNode
        {
        public:
            and_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(2, 2, 0, 0);
                }

            double DoEvaluate()
                {
                if(m_nodes[0]->Evaluate() == 0.0 || m_nodes[1]->Evaluate() == 0.0)
                    return 0.0;
                else
                    return 1.0;
                }
        };

    class and_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "and";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new and_FunctionNode(expr);
                }
        };

    // Or
    //--------------------------------------------------------------------------
    class or_FunctionNode : public FunctionNode
        {
        public:
            or_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(2, 2, 0, 0);
                }

            double DoEvaluate()
                {
                if(m_nodes[0]->Evaluate() == 0.0 && m_nodes[1]->Evaluate() == 0.0)
                    return 0.0;
                else
                    return 1.0;
                }
        };

    class or_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "or";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new or_FunctionNode(expr);
                }
        };

    // Not
    //--------------------------------------------------------------------------
    class not_FunctionNode : public FunctionNode
        {
        public:
            not_FunctionNode(Expression *expr) : FunctionNode(expr)
                {
                SetArgumentCount(1, 1, 0, 0);
                }

            double DoEvaluate()
                {
                if(m_nodes[0]->Evaluate() == 0.0)
                    return 1.0;
                else
                    return 0.0;
                }
        };

    class not_FunctionFactory : public FunctionFactory
        {
        public:
            string GetName() const
                {
                return "not";
                }

            FunctionNode *DoCreate(Expression *expr)
                {
                return new not_FunctionNode(expr);
                }
        };


    } // namespace


// Initialize default functions
void FunctionList::AddDefaultFunctions()
    {
      /*
    #define ADDFUNCTION(name) \
        auto_ptr<FunctionFactory> name ## _func(new name ## _FunctionFactory()); \
        m_functions.push_back(name ## _func.get()); \
        name ## _func.release();
      */

    // according to http://stackoverflow.com/questions/2404115/is-auto-ptr-deprecated
    // unique_ptr should be used instead of auto_ptr
    #define ADDFUNCTION(name) \
        unique_ptr<FunctionFactory> name ## _func(new name ## _FunctionFactory()); \
        m_functions.push_back(name ## _func.get()); \
        name ## _func.release();

    ADDFUNCTION(abs);
    ADDFUNCTION(mod);

    ADDFUNCTION(ipart);
    ADDFUNCTION(fpart);

    ADDFUNCTION(min);
    ADDFUNCTION(max);
    ADDFUNCTION(sqrt);

    ADDFUNCTION(sin);
    ADDFUNCTION(cos);
    ADDFUNCTION(tan);

    ADDFUNCTION(sinh);
    ADDFUNCTION(cosh);
    ADDFUNCTION(tanh);

    ADDFUNCTION(asin);
    ADDFUNCTION(acos);
    ADDFUNCTION(atan);
    ADDFUNCTION(atan2);

    ADDFUNCTION(log);
    ADDFUNCTION(ln);
    ADDFUNCTION(exp);
    ADDFUNCTION(logn);

    ADDFUNCTION(ceil);
    ADDFUNCTION(floor);

    ADDFUNCTION(rand);
    ADDFUNCTION(random);
    ADDFUNCTION(randomize);

    ADDFUNCTION(deg);
    ADDFUNCTION(rad);

    ADDFUNCTION(rect2pol);
    ADDFUNCTION(pol2rect);

    ADDFUNCTION(if); // Preprocess will take care of this beforehand
    ADDFUNCTION(select);

    ADDFUNCTION(equal)
    ADDFUNCTION(above);
    ADDFUNCTION(below);

    ADDFUNCTION(clip);
    ADDFUNCTION(clamp);
    ADDFUNCTION(rescale);

    ADDFUNCTION(poly);

    ADDFUNCTION(and);
    ADDFUNCTION(or);
    ADDFUNCTION(not);
    }
