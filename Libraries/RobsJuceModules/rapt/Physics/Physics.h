#ifndef RAPT_PHYSICS_H_INCLUDED
#define RAPT_PHYSICS_H_INCLUDED

namespace RAPT
{

/** A list of physical constants. Physical simulations may keep a pointer to an object of that type
and outlying code may modify this object. ...hmm...not sure if that's a good idea - maybe each 
simulation should keep the constants it actually needs as members.

References:
-https://physics.nist.gov/cuu/Constants/Table/allascii.txt

todo:
-add more constants
*/

template<class T>
class rsPhysicalConstants
{

public:

  T boltzmann   = 1.38064852e-23;   // J K^-1, Boltzmann constant
  T gravitation = 6.67408e-11;      // m^3 kg^-1 s^-2, Newtonian constant of gravitation 
  T lightspeed  = 299792458;        // m s^-1, speed of light in vacuum
  T planck      = 6.626070040e-34;  // J s, Planck constant
                                    // ...more to come
};


#include "ParticleSystem.h"

}

#endif