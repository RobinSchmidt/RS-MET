#pragma once

namespace LaPackCPP {

typedef long int integer;
//typedef unsigned long int uinteger; // added by Robin Schmidt (guesswork!!)
//typedef char *address;
//typedef short int shortint;
typedef float f2c_real;
typedef double doublereal;
//typedef struct { f2c_real r, i; } f2c_complex;
//typedef struct { doublereal r, i; } doublecomplex;
typedef long int logical;
//typedef short int shortlogical;
typedef long ftnlen;

#define TRUE_ (1)
#define FALSE_ (0)
#define VOID void

//// try to get rid of them:
//#define min(a,b) ((a) <= (b) ? (a) : (b))
//#define max(a,b) ((a) >= (b) ? (a) : (b))


template<class T> T d_lg10(T *);
template<class T> T d_sign(T *, T *);
template<class T> integer i_dnnt(T *);
integer i_len(char *, ftnlen);
template<class T> integer i_nint(T *);
template<class T> T pow_di(T *, integer *);
integer s_cmp(char *, char *, ftnlen, ftnlen);
void s_copy(char *, char *, ftnlen, ftnlen);


}