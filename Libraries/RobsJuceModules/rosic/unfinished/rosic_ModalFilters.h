#pragma once
namespace rosic
{

/** A filter representing a mode that uses four parallel decaying sine filters with different decay
rates but otherwise the same parameters. The different decay rates mix to create interesting 
envelope shapes. It uses a single rsFloat32x4 variable for the signal, i.e. it computes four 
parallel filters at once using SSE2 vector instructions. */

class rsModalFilterFloatSSE2
{

public:

  /** Sets up the mode parameters. Omega is the radian frequency (2*pi*f/fs), the phase is in 
  radians and the decay time constants are in samples. See also rsDampedSineFilter. */
  void setParametersTwoEnvs(
    double omega,   double amplitude, double phase, 
    double attack1, double attack2,   double attackBlend,
    double decay1,  double decay2,    double decayBlend);

  // maybe try another parametrization:
  // freq, amp, phase, att1, dec1, freqSpread, phaseSpread, scaleAtt2, scaleDec2, blend
  // Freq, Amp, Phase, Attack, Decay
  // FreqSpread, PhaseSpread, AttackSpread, DecaySpread, blend
  void setParameters(double omega, double amplitude, double phase, double attack, double decay,
    double deltaOmega = 0, double phaseDelta = 0, double blend = 0.5, 
    double attackScale = 1.0, double decayScale = 1.0);



  /** \name Processing */

  /** Produces the vector of the 4 outputs of the 4 individual decaying sine filters. The actual 
  scalar output sample would be the sum of these 4. */
  //inline rsFloat32x4 getSample(rsFloat32x4 in)
  inline rsFloat32x4 getSample(const rsFloat32x4& in)
  {                                  // CPU-cycles per call (on Athlon 5050e, todo: Core i3-7100U)
    //return getSampleDF1(in);       //   21.2
    return getSampleTDF1(in);        //   12.1      -> most efficient version
    //return getSampleDF2(in);       //   21.3
    //return getSampleTDF2(in);      //   16.1
  }

  /** Produces a scalar output sample that adds up all the 4 decaying sines. When a single mode is
  synthesiszed, you can use this function. When many modes are added, it makes more sense to just
  call the vector function and accumulate the vectors and do just a single sum after the 
  accumulation. */
  inline float getSample(float in) { return getSample(rsFloat32x4(in)).getSum(); }

  /** Convenience function, if you need the output in double-precision format (of, course, the 
  actual precision is single anyway - it's just converted to double). */
  inline double getSample(double in) { return double(getSample(float(in))); }




  /** Resets the state variables to all zeros. */
  void reset() { x1 = v1 = v2 = 0; }

  /** Computes a 4-vector of output samples using a direct form 1 implementation. 
  see: https://ccrma.stanford.edu/~jos/fp/Direct_Form_I.html  */
  //inline rsFloat32x4 getSampleDF1(rsFloat32x4 in)
  inline rsFloat32x4 getSampleDF1(const rsFloat32x4& in)
  {
    rsFloat32x4 y = b0*in + b1*x1 - a1*v1 - a2*v2;
    x1 = in;
    v2 = v1;
    v1 = y;
    return y;
  }

  /** Computes a 4-vector of output samples using a transposed direct form 1 implementation. 
  see: https://ccrma.stanford.edu/~jos/fp/Transposed_Direct_Forms.html  */
  //inline rsFloat32x4 getSampleTDF1(rsFloat32x4 t)  // old
  inline rsFloat32x4 getSampleTDF1(const rsFloat32x4& in)  // new
  {
    //t += v1;               // old
    rsFloat32x4 t = in + v1; // new

    v1 = -a1*t + v2;
    v2 = -a2*t;
    rsFloat32x4 y = b0*t + x1; // x1 == s3
    //y  = b0*t + x1; // use member instead of stack variable - no good idea, increases cpu load
    x1 = b1*t;
    return y;
  }

  /** Computes a 4-vector of output samples using a direct form 2 implementation. 
  see: https://ccrma.stanford.edu/~jos/fp/Direct_Form_II.html  */
  //inline rsFloat32x4 getSampleDF2(rsFloat32x4 tmp) // old
  inline rsFloat32x4 getSampleDF2(const rsFloat32x4& in) // new
  {
    //tmp -= (a1*v1 + a2*v2);  // old
    rsFloat32x4 tmp = in - (a1*v1 + a2*v2); // new

    rsFloat32x4 out = b0*tmp + b1*v1;
    v2 = v1;
    v1 = tmp;
    return out;
  }

  /** Computes a 4-vector of output samples using a transposed direct form 2 implementation. 
  see: https://ccrma.stanford.edu/~jos/fp/Transposed_Direct_Forms.html */
  //inline rsFloat32x4 getSampleTDF2(rsFloat32x4 in)
  inline rsFloat32x4 getSampleTDF2(const rsFloat32x4& in)
  {
    rsFloat32x4 y = b0*in + v1;
    v1 = b1*in - a1*y + v2;
    v2 = -a2*y;
    return y;
  }

protected:

  //rsFloat32x4 y; // trying to reduce cpu load in TDF1 by using a member instead of stack 
  // allocated variable - but that actually increases cpu load

  rsFloat32x4 x1 = 0, v1 = 0, v2 = 0, b0 = 0, b1 = 0, a1 = 0, a2 = 0;
  // Maybe make a subclass for DF1 and TDF1 that need to store x1 - save one memory cell for the
  // other forms. Actually, it turned out that TDF1 is the most cpu efficient version and this one
  // needs the additional storage cell. but maybe when we later use thousands of these filters, the
  // memory footprint becomes more relevant. Do performance tests under such realistic conditions
  // later - it may turn out to be advantegeous to switch to TDF2 with its lower memory footprint.
};

// -maybe templatize and make two explicit instantiations for rsFloat32x4 (SSE2) and 
//  rsFloat32x16 (AVX2), maybe have two template arguments - for the scalar and vector type
// -maybe also make a rsFloat64x2 version - of course, that can't use 4 filters, but maybe 2 is 
//  enough in some situations. An/or make a rsFloat64x4 version
// -maybe rename to rsModalFilterFloat32x4

}