#ifndef RAPT_CONIC_SECTION_H_INCLUDED
#define RAPT_CONIC_SECTION_H_INCLUDED

/** This is a class for dealing with conic sections, represented by the implicit equation: 
A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0. */

template<class T>
class ConicSection
{

public:

  //static void 

protected:

  T A = 1, B = 0, C = 1, D = 0, E = 0, F = -1; // init as unit circle

};

#endif
