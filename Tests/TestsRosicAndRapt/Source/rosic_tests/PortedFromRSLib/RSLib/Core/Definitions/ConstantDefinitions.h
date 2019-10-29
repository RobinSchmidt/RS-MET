#ifndef RS_CONSTANTDEFINITIONS_H
#define RS_CONSTANTDEFINITIONS_H

namespace RSLib
{

// mathematical constants
// \todo perhaps we should use static const doubles instead of macros?

#define PI 3.1415926535897932384626433832795
#define EULER 2.7182818284590452353602874713527  // maybe rename to EULER_NUMBER
#define EULER_CONSTANT 0.5772156649015328606065120900824024310421
#define SQRT2 1.4142135623730950488016887242097

#define SQRT2_INV 0.70710678118654752440084436210485  // rename: INV_SQRT2 - get rid
#define RS_SQRT2_INV 0.70710678118654752440084436210485

#define LN10 2.3025850929940456840179914546844
#define LN10_INV 0.43429448190325182765112891891661   // rename: INV_LN10
#define LN2 0.69314718055994530941723212145818
#define LN2_INV 1.4426950408889634073599246810019     // rename: INV_LN2
#define SEMITONE_FACTOR 1.0594630943592952645618252949463  // move to the Audio module
#define TINY FLT_MIN
#define EPS DBL_EPSILON  // deprecate - use RS_EPS(double) instead

static const double rsInfDouble          = std::numeric_limits<double>::infinity();
static const double rsQuietNaNDouble     = std::numeric_limits<double>::quiet_NaN();
static const double rsSignalingNaNDouble = std::numeric_limits<double>::signaling_NaN();

#define RS_INF(T) (std::numeric_limits<T>::infinity())
#define RS_MIN(T) (std::numeric_limits<T>::min())
#define RS_MAX(T) (std::numeric_limits<T>::max())
#define RS_NAN(T) (std::numeric_limits<T>::quiet_NaN())
#define RS_EPS(T) (std::numeric_limits<T>::epsilon())

// powers of two
#define RS_POW2_0 1
#define RS_POW2_1 2
#define RS_POW2_2 4
#define RS_POW2_3 8
#define RS_POW2_4 16
#define RS_POW2_5 32
#define RS_POW2_6 64
#define RS_POW2_7 128
#define RS_POW2_8 256
#define RS_POW2_9 512
#define RS_POW2_10 1024
#define RS_POW2_11 2048
#define RS_POW2_12 4096
#define RS_POW2_13 8192
#define RS_POW2_14 16384
#define RS_POW2_15 32768


inline double rsSqrt(double x) { return ::sqrt(x); }

}

#endif // #ifndef RSLib_ConstantDefinitions_h
