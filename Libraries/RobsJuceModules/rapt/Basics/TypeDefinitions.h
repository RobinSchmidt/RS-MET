#ifndef RAPT_TYPEDEFINITIONS_H_INCLUDED
#define RAPT_TYPEDEFINITIONS_H_INCLUDED

// todo: move them to rosic - rapt is not supposed to use typed code anyway - everything should be
// generic/templated

/** 8-bit signed integer */
typedef signed char rsInt8;

/** 8-bit unsigned integer */
typedef unsigned char rsUint8;

/** 16-bit signed integer */
typedef signed short rsInt16;

/** 16-bit unsigned integer */
typedef unsigned short rsUint16;

/** 32-bit signed integer */
typedef signed int rsInt32;

/** 32-bit unsigned integer */
typedef unsigned int rsUint32;

/** 64-bit signed integer  */
typedef signed long long rsInt64;  // not supported by ISO C++

/** 64-bit unsigned integer */
typedef unsigned long long rsUint64;

/** 32-bit floating point number */
typedef float rsFloat32;

/** 64-bit floating point number */
typedef double rsFloat64;

/** 80-bit floating point number */
typedef long double rsFloat80;

/** Pointer to a function that takes a double parameter and returns a double */
typedef double(*FunctionPointerDoubleToDouble) (double);

/** Pointer to a function that takes a 2 double parameters and returns a double */
typedef double(*FunctionPointer2DoublesToDouble) (double, double);

/** Pointer to a function that takes a 3 double parameters and returns a double */
typedef double(*FunctionPointer3DoublesToDouble) (double, double, double);

///** doubles, aligned at 64-bit (8 byte) boundaries - conflict with rosic - not needed in RAPT */
//#ifdef _MSC_VER
//typedef __declspec(align(8)) double doubleA;
//#else
//typedef double doubleA; // something to do here...
//#endif

// todo: define SIMD types rsFloat32x4, rsFloat64x2, etc. (vector of 4 floats, 2 doubles, ...)
// these should be wrapped into classes which should be subclasses of the respective compiler
// intrinsic vector types. (arithmetic) operators can be defined (inline) which map directly to
// the corresponding intrinsic vector functions

/** An empty struct meant to be used in places where the compiler demands a template parameter, but 
we actually do not want any data to be stored in a field of that type. An example could be an 
instantiation of rsGraph, where we don't need any data to be associated with edges and/or vertices. 
We could pass the rsEmptyType for the template parameter for edge- and/or vertex data type. */
struct rsEmptyType 
{ 
  rsEmptyType()    {};
  rsEmptyType(int) {};
};
// The constructors are meant to be able to define default arguments for function calls and make
// the work also with the empty struct (todo: check, if that is actually necessarry)

#endif
