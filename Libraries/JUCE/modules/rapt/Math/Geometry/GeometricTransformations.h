#ifndef RAPT_GEOMETRICTRANSFORMATIONS_H_INCLUDED
#define RAPT_GEOMETRICTRANSFORMATIONS_H_INCLUDED

/** Rotates a 3D vector given by x,y,z around the three axes by the given angles rx, ry, rz (in
that order - which is important because rotations are not commutative). */
template<class T>
void rsRotateXYZ(T* x, T* y, T* z, T rx, T ry, T rz);
 // not yet tested

#endif