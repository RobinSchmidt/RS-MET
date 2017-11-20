#ifndef RS_DOUBLEPENDULUM_H
#define RS_DOUBLEPENDULUM_H

namespace RSLib
{

  /**

  This class implements the differential equation system of a double pendulum. This system exhibits
  chaotic motion.

  References: 
  (1) http://scienceworld.wolfram.com/physics/DoublePendulum.html

  \todo
  When using different masses, like l1 = 1, l2 = 0.8, m1 = 1, m2 = 1.2, it seems like the algorithm 
  becomes unstable. Maybe this is a "stiff" problem and we should try implicit methods? That would 
  require to extend the differential equation solver - see Numerical Recipies...

  \todo:
  -include a damping:
   -subtract a term proportional to the respective angular velocity from the derivative with repect
    to the momenta?
  -include a driving force (input signal)

  */

  class RSLib_API rsDoublePendulum : public rsDifferentialEquationSystemDbl
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. */
    rsDoublePendulum();   


    /** \name Setup */

    /** Sets the length of the first edge. */
    void setLength1(double newLength);

    /** Sets the length of the second edge. */
    void setLength2(double newLength);

    /** Sets the mass of the first bob. */
    void setMass1(double newMass);

    /** Sets the mass of the second bob. */
    void setMass2(double newMass);

    /** Sets the integration step size */
    void setStepSize(double newStepSize); // maybe factor out

    /** Sets up the state of the pendulum. It consists of the the two angles theta1/2 and momenta
    momentum1/2. This function can be used to set up the initial values before the system is
    iterated. */
    void setState(double theta1, double theta2, double momentum1, double momentum2);


    /** \name Inquiry */

    /** Returns the values of the state variables. */
    void getState(double *theta1, double *theta2, double *momentum1, double *momentum2);

    /** Returns the two angles in a1, a2 */
    void getAngles(double *a1, double *a2);

    /** Returns the two momenta in m1, m2 */
    void getMomenta(double *m1, double *m2);

    /** Returns the cartesian coordinates of both bobs of the pendulum. */
    void getPendulumCoordinates(double *x1, double *y1, double *x2, double *y2);


    /** \name Processing */

    /** Updates the state of the pendulum, i.e. advances the time by the stepsize. */
    void updateState();



    /** \name Misc */


    /** Overriden function, inherited from rsDifferentialEquationSystemDbl - computes the value
    of the partial derivatives at the current state. */
    virtual rsVector<double> f(const double &x, const rsVector<double> &y);




  protected:

    /** \name Internal Functions */



    /** \name Data */

    // user parameters:
    double l1, l2;  // the two lengths
    double m1, m2;  // the two masses
      // todo: include r1, r2 - two coefficients of friction

    double h;       // step size
      // maybe factor out


  };

  //-----------------------------------------------------------------------------------------------
  // inlined functions:



  // todo: make a second implementation based on teh equations here:
  // http://www.myphysicslab.com/dbl_pendulum.html
  // they operate directly on w1, w1, t1, t2 instead (w1, w2: derivatives of t1, t2)
  // this site has more physics simulations:
  // http://www.myphysicslab.com/index.html


}

#endif
