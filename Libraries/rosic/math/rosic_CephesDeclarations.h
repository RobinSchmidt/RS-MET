/**

This file includes all relevant declarations for special functions from the Cephes math library.

*/

#ifndef rosic_CephesDeclarations_h
#define rosic_CephesDeclarations_h

extern "C" int ellpj(double u, double m, double *sn, double *cn, double *dn, double *ph);

extern "C" double jv(double v, double x);

//extern "C" int sgngam = 0;
//extern int sgngam;

#ifndef LINUX // interferes with /usr/include/bits/mathcalls.h on linux -> define linux somewhere in the compiler-settings on linux
extern "C" int sincos(double x, double *s, double *c, int flg = 0);
#endif

#endif // end of #ifndef rosic_Math_h
