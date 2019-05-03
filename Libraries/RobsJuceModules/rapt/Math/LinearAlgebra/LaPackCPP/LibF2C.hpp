#pragma once

namespace LaPackCPP {

typedef long int integer;    // get rid - replace with int
//typedef unsigned long int uinteger; // added by Robin Schmidt (guesswork!!)
//typedef char *address;
//typedef short int shortint;
typedef float f2c_real;     // replace with float
typedef double doublereal;  // replace with double
//typedef struct { f2c_real r, i; } f2c_complex;
//typedef struct { doublereal r, i; } doublecomplex;
typedef long int logical;   // replace with int or bool
//typedef short int shortlogical;
typedef long ftnlen;        // replace with int

//#define TRUE_ (1)
//#define FALSE_ (0)
//#define VOID void

static const logical TRUE_ = 1;
static const logical FALSE_ = 0;


//// try to get rid of them:
//#define min(a,b) ((a) <= (b) ? (a) : (b))
//#define max(a,b) ((a) >= (b) ? (a) : (b))


template<class T> T d_lg10(T *);
template<class T> T d_sign(T *, T *);
template<class T> integer i_dnnt(T *);

//template<class T> integer i_nint(T *)
template<class T>
inline integer i_nint(T *x)
{
  return (integer)(*x >= 0 ? floor(*x + .5) : -floor(.5 - *x));
}

template<class T> T pow_di(T *, integer *);

integer i_len(char *, ftnlen);
integer s_cmp(char *, char *, ftnlen, ftnlen);
void s_copy(char *, char *, ftnlen, ftnlen);


}