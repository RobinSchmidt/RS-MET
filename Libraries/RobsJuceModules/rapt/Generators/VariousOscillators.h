#ifndef RAPT_VARIOUSOSCILLATORS_H_INCLUDED
#define RAPT_VARIOUSOSCILLATORS_H_INCLUDED



template<class T>
class rsSineOscillatorNaive
{

public:


  void setAmplitude(T newAmplitude) { amp = newAmplitude; }

  void setOmega(T newOmega) { omega = newOmega; }

  void setPhase(T newPhase) { phase = newPhase;   }

  void modulatePhase(T amount) { phase += T(2*PI)*amount; } // verify factor 2*PI

  inline T getSample()
  {
    T y = amp * sin(phase);
    phase += omega;
    return y;
  }

  void reset(T startPhase = T(0))
  {
    phase = startPhase;
  }

protected:

  T amp   = T(1);
  T phase = T(0);
  T omega = T(0);

};


// maybe rename to rsSineOscillatorRecursive - we may also have a naive one - they should both 
// behave the same way, but be different with regard to which operations are efficient and which 
// are expensive - both should allow FM and PM - maybe even have a class that somehow automatically
// dispatches between both implementations depending on the conditions - PM/FM is cheaper with a 
// naive implementation, just producing a steady sine is cheaper with the recursive implementation
// so we need to detect, if PM/FM is going on - maybe that should be done in the jura processor.
// there we can inquire, if a modulator is connected to freq and/or phase

template<class T>
class rsSineOscillatorRecursive : public rsSineIterator<T>
{

public:

  //using rsSineIterator::rsSineIterator;  // inherit constructors


  void setAmplitude(T newAmplitude) { amp = newAmplitude; }

  void setOmega(T w, bool fixPhase = true) 
  { 
    this->a1 = 2.0*cos(w);
    if(fixPhase)
    {
      T p = this->getPhase();
      this->s1 = this->a1*sin(p-    w);
      this->s2 = this->a1*sin(p-2.0*w);
      // but what if w < 0? i think, we should then swap s1,s2 - this needs tests - we want to be
      // able to do through-zero FM
    }


    //setup(newOmega, getPhase(), T(1)); 
  }
  // getPhase should reconstruct the phase from the state - needs asin and then figure out if we 
  // are in the ascending or descending part and possibly add an offset of pi, maybe:
  //   y = getValue();
  //   p = asin(y);
  //   if(y < s1)
  //     p += PI;
  //

  inline T getSample()
  {
    return amp * rsSineIterator<T>::getValue();
  }


protected:

  T amp = T(1);

};


//=================================================================================================

/** An oscillator based on morphing between saw-up/triangle/saw-down waveforms. 

todo: produce info for blep/blamp anti-aliasing, make a DualTriSawOsc - drive that controls two
TriSaws and allows them to interact - in particular hardsync...or maybe even a TripleTriSawOsc - it
may be cool to have two independent sync-masters ...or maybe that sync-stuff could even be 
implemented in liberty
*/

template<class T>
class rsTriSawOscillator
{

public:

  rsTriSawOscillator() { updateTriSawCoeffs(); }

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets the phase-increment per sample. Should be set to frequency/sampleRate. (times 2?) */
  inline void setPhaseIncrement(T newIncrement) { inc = newIncrement; }

  /** Sets the asymmetry parameter of the waveform which should be a number in -1..+1. For -1, you
  will get a downward sawtooth, for 0 a triangle and for +1 an upward sawtooth wave (assuming the
  bending and sigmoidity parameters are at their neutral settings of 0). */
  inline void setAsymmetry(T newAsymmetry) 
  { 
    h = T(0.5)*(newAsymmetry+1);
    updateTriSawCoeffs();
  }

  /** Sets the bending parameter for the upward half-wave. The value is in -1..+1 where negative 
  values let the curve bend inward (toward the time axis) and positive values outward (toward a 
  square wave). */
  inline void setAttackBending(T newParam) { t1 =  newParam; }

  /** Bending parameter for downward half-wave. */
  inline void setDecayBending(T newParam) { t2 = -newParam; } 

  /** Sets the amount of sigmoidity or s-shapedness of the upward half-wave. A value of +1 leads to
  an s-shape with zero derivative at the corner point (i.e. perfectly rounded corners) and negative 
  values make the curve more spikey. */
  inline void setAttackSigmoid(T newParam) { s1 = T(-0.5)*newParam; }

  /** Sigmoidity of downward half-wave. */
  inline void setDecaySigmoid(T newParam) { s2 = T(-0.5)*newParam; }


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Computes an asymmetry parameter which results in a (downward) transition with the given 
  number of samples according to the current increment setting. If you want to have the upward 
  transition last this number of samples, use minus this value. It's useful for ensuring a minimum
  absolute transition time. */
  inline T asymForTransitionSamples(T numTrasitionSamples)
  {
    return rsMax(T(0), T(1) - T(2)*numTrasitionSamples*inc);
  }

  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Updates the phase-variable (called in getSample) */
  inline void updatePhase()
  {
    p += inc;
    while(p >= 1) p -= 1;
    while(p <  0) p += 1;
  }
  // maybe factor out into a Phasor class - that could also have a phase-input (for 
  // phase-modulation)...but this would require the wrap-around code in the update *and* in the
  // read-out code

  /** Given a value between -1..1 (supposed to be the raw TriSaw output), this applies the shape
  formula with parameters t for the tension and s for the sigmoidity. */
  static inline T shape(T x, T t, T s)
  {
    x = ((1-s) + s*x*x)*x;  // apply sigmoidity
    T d = t*x+1;            // denominator
    if(d != T(0))
      return (x+t) / d;     // apply bending
    else
      return T(1);          // is this limit value always the correct one? ...verify...
  }
  // used in romos::TriSawModule::process
  // without the zero-check, there was a crash when AttackBend == 1, strangely DecayBend = +-1 gave 
  // no problems and AttackBend = -1 also worked

  inline T shape1(T x) { return shape(x, t1, s1); }
  inline T shape2(T x) { return shape(x, t2, s2); }
  //inline T shape1(T x) { return x; }
  //inline T shape2(T x) { return x; }

  //inline T getFromPhase(T 

  inline T getSample()
  {
    T y;
    if(p < h)
      y = shape1(a0 + a1*p);  // upward section
    else 
      y = shape2(b0 + b1*p);  // downward section

    // y = clip(y/(1-s), T(-1), T(+1)); // just an idea - s: squarishness parameter

    rsAssert(rsIsFiniteNumber(y));
    updatePhase();
    return y;
  }
  // maybe factor out a raw TriSawOscillator (without the shape stuff) and realize the shape stuff
  // in a subclass...or have a member function getSampleRaw that doesn't apply shaping
  // 

  inline void reset() { p = startPhase; }

  /*
  // used in romos:
  static void getTriSawCoeffs(T h, double *a0, double *a1, double *b0, double *b1)
  {
    *a0 = -1;
    *a1 = 2 / h;
    *b0 = (1+h)/(1-h);
    *b1 = -1 - *b0;

    // todo: we need to catch cases when h = 0 or 1-h = 0
    // if(h < eps)
    //   a1 = 0;
    // if(1-h < eps)
    //   b1 = 0;
    // or something
  }
  static T getFromTime(T in, T Asym, T AttBend, T AttSigm, T DecBend, T DecSigm)
  {
    double p = fmod(in, 1);

    return p; // preliminary
  }
  */

protected:

  /** Updates the coefficients for the creation of basic triangle/sawtooth waveform before the 
  waveshaper is applied. */
  void updateTriSawCoeffs();



  T startPhase = 0;  
  T p = 0;     // current phase in 0..1
  T h = 0.5;   // time instant of the "half-cycle" point, i.e. the switch between up/down
  T inc = 0;   // phase increment

  // variables for the trisaw waveform before applying the waveshaper
  T a0 =  0; 
  T a1 =  2;
  T b0 =  2;
  T b1 = -2;
  // maybe rename to u0,u1, d0,d1 (for up/down)


  T t1 = 0, s1 = 0;  // tension and sigmoidity for first (upward) half-wave
  T t2 = 0, s2 = 0;  // same for downward half-wave
};

//=================================================================================================

/** Oscillator based on an ellipse in the xy-plane.

Credits:
 -Xoxos

more info:
 https://www.kvraudio.com/forum/viewtopic.php?p=6656238#p6656238
 https://gitlab.com/Hickler/Soundemote/issues/67
 https://github.com/RobinSchmidt/RS-MET/issues/72
play with parameters:
  https://www.desmos.com/calculator/7h9mknbv3q
i think, it works as follows:
-create x,y values on a circle (standard rotating phasor in the plane)
-convert x,y to values on an arbitrary ellipse
-(-project onto the x- and y-axis (i.e. take the x- and y-value))...really? 

the name is actually not really suitable - maybe rename to XoxosOscillator

...and make an EllipseOscillator as well

*/

template<class T>
class rsEllipseOscillator
{

public:

  inline void setA(T newA) { A = newA; }  // rename to setOffset

  inline void setB(T newB) { B = newB; }  // rename to setRotation

  inline void setC(T newC) { C = newC; }  // rename to setScale

  inline void setOmega(T newOmega) { w = newOmega; }

  inline void setRotationDegrees(T newRotation) { B = rsDegreeToRadiant(newRotation); }
   // no - this is not a rotation...maybe phase-shift?


  inline void updatePhase()
  {
    p += w;
    while(p > 2*PI)
      p -= 2*PI;

    // maybe do also a wraparound at 0 -> allow negative frequencies
  }

  // maybe have a function that returns x and y separately - maybe that's useful as stereo signal?
  // if not, it is certainly helpful to figure out what the osc is doing - and then we can do 
  // mid/side (re)mixing

  inline void getSamplePair(T* x, T* y)
  {
    T c  = cos(p);  // x
    T s  = sin(p);  // y

    // precompute these:
    T sB = sin(B);
    T cB = cos(B);

    T Ac = A + c;    // x += a -> shift x
    T Cs = C * s;    // y *= c -> scale y


    T scl = 1 / sqrt(Ac*Ac + Cs*Cs);  // normalizer
    if(renorm != 1)
      scl = pow(scl, renorm);

    *x = scl*Ac*cB;  // phase(?) and
    *y = scl*Cs*sB;  // normalize

    updatePhase();
  }

  inline T getSample()
  {
    T x, y;
    getSamplePair(&x, &y);
    return 2*((1-mix)*x + mix*y); // precompute coeffs 2*(1-mix), 2*mix
    //return x + y;   // use mix*x + (1-mix)*y
  }

  inline void reset() { p = startPhase; }

protected:

  T startPhase = 0;  
  T p = 0;   // current phase
  T w = 0;   // normalized radian frequency
  // have independent phases and frequencies for x and y

  // waveshape parameters:
  T A = 0; 
  T B = 0;
  T C = 1;

  T renorm = 1;

  T mix = 0.5;
};

#endif