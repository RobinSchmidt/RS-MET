#ifndef RAPT_DIFFERENTIALEQUATIONSYSTEM_H
#define RAPT_DIFFERENTIALEQUATIONSYSTEM_H

/** This is a class for representing a system of ordinary differential equations of the form:

\f[ \frac{d\mathbf{y}}{dx} = \mathbf{f}(x, \mathbf{y}) \f].

It is intended to be used as baseclass. A specific system of differential equations is
implemented by creating a subclass which must implement the vector-valued function
\f$\mathbf{f}\f$. The baseclass itself calls that overriden function inside the
baseclass-provided numerical integration routines in a kind of "template-method" pattern.

The class is templatized on the type of the independent scalar variable \f$x\f$ and the
element-type of the dependent vector-variable \f$\mathbf{y}\f$. Typically, you'll want to use
floats/doubles for both, so an explicit instantiation exists for double. But if needed, you may
also instantiate it for complex numbers or whatever other datatype you may need.

The class maintains the current values of the independent variable \f$x\f$ and the dependent
variable \f$\mathbf{y}\f$. These can be set and retrieved from client-code as well as advanced by
some stepsize \f$h\f$ using different numerical integration methods.

\section Background Modeling of Dynamical Systems
Often, the independent variable is time in which case the equation can be interpreted as a
prescription to calculate the velocity of a point that moves around in some multidimensional
space. If the velocity depends only on the position itself (the \f$\mathbf{y}\f$-vector) but
not explicitly on the time instant (the scalar \f$x\f$), the system is said to be autonomous.
You would implement an autonomous system by simply ignoring the value of \f$x\f$ inside your
implementation of \f$\mathbf{f}(x, mathbf{y})\f$.

\section References
(1): Numerical Recipes in C, 2nd Ed, Ch. 16
(2): Höhere Mathematik für Naturwissenschaftler und Ingenieure, 2nd Ed, Ch. 6

The implementation of this class stays close to the notation used in (1), so in order to
understand how it works internally, it is recommended to consult that book.

...actually, this implementation should be used only as prototype - it does a lot of unnecessary
memory allocations - this makes the code neat and convenient but is really inefficient
-> move to prototypes and use the implementation from the GNUPlotCPP codebase for production code

\todo implement more iteration formulas, maybe even implicit methods (those require solution of
a linear system per step), adpative stepsize control, predictor/corrector methods, whatever.....
but maybe defer this to subclasses */

template<class TypeX, class TypeY> // rename to TX, TY
class rsDifferentialEquationSystem
{

public:

  /** \name Construction/Destruction */

  /** Constructor. Initializes the independent variable to zero. The dependent vector-variable y
  is still zero-dimensional in this baseclass, so initialization (including setting its
  dimensionality) is deferred to subclasses. */
  rsDifferentialEquationSystem()
  {
    x = TypeX(0);
  }


  /** \name Setup */

  /** Sets the current value of the independent variable. */
  inline void setX(const TypeX &newX)
  {
    x = newX;
  }

  /** Sets the current position in phase-space. \todo maybe rename to setState */
  inline void setY(const rsVector<TypeY> &newY)
  {
    y = newY;
  }

  /** Sets the i-th element of the current coordinate y-vector to a new value. */
  inline void setElementOfY(int i, TypeY newValue)
  {
    y[i] = newValue;
  }


  /** \name Inquiry */

  /** Returns the current value of the independent variable. */
  inline TypeX getX()
  {
    return x;
  }

  /** Returns the current position in phase-space. */
  inline rsVector<TypeY> getY()
  {
    return y;
  }

  /** Returns the i-th element of the current y-vector. */
  inline TypeY getElementOfY(int i)
  {
    return y[i];
  }

  /** Returns the number of dimensions of the phase-space. */
  inline int getNumDimensions()
  {
    return y.dim;
  }


  /** \name Derivative Calculation */

  /** This is supposed to be overriden in subclasses to provide the derivative vector
   \f[\mathbf{f}(x, mathbf{y}) \f]. */
  virtual rsVector<TypeY> f(const TypeX &x, const rsVector<TypeY> &y) = 0;


  /** \name Iteration Functions */

  /** Performs a forward Euler step using stepsize h.  See (1), p. 710, Eq. 16.1.1 */
  void stepEuler(TypeX h)
  {
    y += h * f(x, y);
    x += h;
  }

  /** Performs a midpoint-method step using stepsize h. See (1), p. 710, Eq. 16.1.2 */
  void stepMidpoint(TypeX h)
  {
    rsVector<TypeY> k1, k2;

    k1 = h * f(x, y);
    k2 = h * f(x+h/2, y+k1/2);

    y += k2;
    x += h;
  }

  /** Performs a 2nd order Heun step using stesize h. See (2), p. 493, Eq. 6.102
  ...hmmm is this the same as Midpoint? experiments suggest so */
  void stepHeun2(TypeX h)
  {
    rsVector<TypeY> k1, k2;

    k1 = f(x, y);
    k2 = f(x+h, y+h*k1);

    y += (h/2) * (k1 + k2);
    x += h;
  }

  /*
  // commented out because it produces a very large error. it's probably still buggy.
  void stepHeun3(TypeX h)
  {
    rsVector<TypeY> k1, k2, k3;

    k1 = f(x              , y                 );
    k2 = f(x + (1.0/3.0)*h, y + h*(1.0/3.0)*k1);
    k3 = f(x + (2.0/3.0)*h, y + h*(2.0/3.0)*k2);

    y += h * ((1.0/4.0)*k1 + (3.0/4.0)*k3);
    x += h;
  }
  */

  /** Performs a 4-th order Runge/Kutta step using stepsize h. See (1), p. 711, Eq. 16.1.3 */
  void stepRungeKutta4(TypeX h)
  {
    rsVector<TypeY> k1, k2, k3, k4;

    k1 = h * f(x, y);
    k2 = h * f(x+h/2, y+k1/2);
    k3 = h * f(x+h/2, y+k2/2);
    k4 = h * f(x+h, y+k3);

    y += k1/6 + k2/3 + k3/3 + k4/6;
    x += h;
  }

  /** Performs a 5-th order Runge-Kutta step and returns a local error estimate. The estimated
  error is the difference between the actual step done (the error of which is 5th order) and an
  embedded Runge-Kutta step of 4th order. So, the returned error estimate actually estimates the
  error of the 4th order step, nonetheless we take the 5th order step. That amounts to taking a
  step that is actually supposed to have higher precision than the error estimate indicates but
  we can't really count on that because we can't estimate the error of the 5th order step itself.
  So the error-estimate is pessimistic, which is typically what we want. See (1), p. 714 ff. */
  rsVector<TypeY> stepCashKarpWithErrorEstimate(TypeX h)
  {
    rsVector<TypeY> k1, k2, k3, k4, k5, k6;

    // the constants from (1), p. 717:
    static const double a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875;
    static const double c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0; //c2=c5=0
    static const double
      d1=c1-2825.0/27648.0, d3=c3-18575.0/48384.0, d4=c4-13525.0/55296.0, d5=-277.0/14336.0,
      d6=c6-0.25;  // di = (ci - ci*), d2 = 0
    static const double
      b21=0.2,
      b31=3.0/40.0, b32=9.0/40.0,
      b41=0.3, b42 = -0.9, b43=1.2,
      b51 = -11.0/54.0, b52=2.5, b53 = -70.0/27.0, b54=35.0/27.0,
      b61=1631.0/55296.0, b62=175.0/512.0, b63=575.0/13824.0, b64=44275.0/110592.0,
      b65=253.0/4096.0;

    // (1), p. 716, Eq. 16.2.4:
    k1 = h * f(x, y);
    k2 = h * f(x+a2*h, y+b21*k1);
    k3 = h * f(x+a3*h, y+b31*k1+b32*k2);
    k4 = h * f(x+a4*h, y+b41*k1+b42*k2+b43*k3);
    k5 = h * f(x+a5*h, y+b51*k1+b52*k2+b53*k3+b54*k4);
    k6 = h * f(x+a6*h, y+b61*k1+b62*k2+b63*k3+b64*k4+b65*k5);

    y += c1*k1+c3*k3+c4*k4+c6*k6;   // c2*k2, c5*k5 missing because c2=c5=0
    x += h;

    // return estimated error via (1), p. 716, Eq. 16.2.6:
    return d1*k1+d3*k3+d4*k4+d5*k5+d6*k6;  // d2*k2 missing because d2=0
  }

protected:


  /** \name Data Members */

  TypeX           x;  ///< current value of the independent variable
  rsVector<TypeY> y;  ///< current position in phase-space

};

// explicit instantiation for doubles:
typedef rsDifferentialEquationSystem<double, double> rsDifferentialEquationSystemDbl;

/*

maybe implement variations of the steppers based on non-newtonian derivatives:
https://hal.archives-ouvertes.fr/file/index/docid/945788/filename/nncam.pdf
https://books.google.de/books?id=RLuJmE5y8pYC&printsec=frontcover&dq=%22Non-Newtonian+Calculus%22&hl=en&sa=X&ved=0ahUKEwjD7b3x6P3iAhUDZ1AKHXmHBFMQ6AEIKjAA#v=onepage&q=%22Non-Newtonian%20Calculus%22&f=false
i think, the Euler step should be modified as follows:
original: y[n+1] = y[n] + h * f(y[n])
new:      y[n+1] = y[n] * exp(h * f(y[n]) / y[n])
...and for the higher order steps (Runge-Kutta, etc.), similar formulas can be derived, see here:
https://www.hindawi.com/journals/aaa/2015/594685/
https://arxiv.org/pdf/1402.2877.pdf

maybe, we should have a fully general Euler step, based on Eq.21 here:
"Generalized Runge-Kutta Method with respect to the Non-Newtonian Calculus":
https://www.hindawi.com/journals/aaa/2015/594685/
...and similarly for other steps such as Runge Kutta
it needs 4 functions: alpha, alphaInverse, beta, betaInverse - in case of the bigeometric calculus
they woud be exp,log,exp,log, in case of geometric calculus: id,id,exp,log (or exp,log,id,id?)



*/


#endif
