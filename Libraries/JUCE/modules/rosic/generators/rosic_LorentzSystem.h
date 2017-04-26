#ifndef rosic_LorentzSystem_h
#define rosic_LorentzSystem_h

// rosic-indcludes:
#include "rosic_SineOscillatorStereo.h"

namespace rosic
{

  /** Implements a discrete time version of the Lorentz system defined by the system of 
  differential equations: dx/dt = sigma*(y-x), dy/dt = x*(rho-z) - y, dz/dt = x*y - beta*z where
  sigma, rho, beta are the parameters of the system


  
  \todo use LaTeX markup in the comment
  \todo write more sophisticated state iteration functions based on runge/kutta or something

  */

  class LorentzSystem
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. */
    LorentzSystem();

    /** \name Setup */

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets up a kind of pseudo-frequency that determines how fast the system will run. */
    void setPseudoFrequency(double newPseudoFrequency);

    /** Sets the internal state of the system consisting of the 3 state variables x, y, z. */
    INLINE void setState(double x, double y, double z)
    {
      this->x = x;
      this->y = y;
      this->z = z;
    }

    INLINE void setSigma(double sigma)
    {
      this->sigma = sigma;
    }

    INLINE void setRho(double rho)
    {
      this->rho = rho;
    }

    INLINE void setBeta(double beta)
    {
      this->beta = beta;
    }

    /** \name Inquiry */

    INLINE void getState(double *x, double *y, double *z)
    {
      *x = this->x;
      *y = this->y;
      *z = this->z;
    }

    /** \name Processing */

    /** Iterates the internal state one step forward. */
    INLINE void iterateState()
    {
      double dx = sigma*(y-x);
      double dy = x*(rho-z) - y;
      double dz = x*y - beta*z;
      x += h*dx;
      y += h*dy;
      z += h*dz;
    }


  protected:

    INLINE void updateStepSize()
    {
      const double c = 1.0;    
        // \todo tweak this such that pseudoFrequency indeed coincides with perceived frequency

      h = c * pseudoFrequency / sampleRate; 
    }

    double x, y, z;                        // state variables
    double sigma, rho, beta, h;            // internal parameters
    double sampleRate, pseudoFrequency;    // user parameters

  };

}

#endif
