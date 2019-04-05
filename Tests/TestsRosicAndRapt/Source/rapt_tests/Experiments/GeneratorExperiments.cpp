#include "GeneratorExperiments.h"
using namespace RAPT;
//using namespace rosic;

void blep()  // rename to blit
{
  //double f  = 1000;    // signal frequency
  //double fs = 44100;   // sample rate
  //double length = 0.02; // length in seconds
  //double period = fs/f;
  //int N = (int) (fs*length);

  int N = 10;
  double period = 10.25;

  typedef rsStepBandLimiter<double, double> SBL;
  SBL sbl;
  sbl.setLength(3);

  std::vector<double> x(N), y(N); // naive and anti-aliased signal


  // preliminary test:
  sbl.setLength(3);
  rsSetZero(y);
  sbl.addImpulse(0.75, 1);       // delay is 0.75, fractional position of spike is 0.25
  y[0] = sbl.getSample(1);
  for(int n = 1; n < N; n++)
    y[n] = sbl.getSample(0);
  rsPlotVector(y);
  // y looks wrong - one of the two positive samples from the mainlobe is missing - which one?
  // the one from the delayline or the one from the corrector? the height of the existing spike is 
  // 0.9 - this is the value from the delayline

  // with sincLength = 3:
  // we want to see an impulse at 2.25 represented by a large positive value at sample 2, a smaller
  // positive sample at sample 3, 1 and 4 should be negative and 0 and 5 positive again

  // maybe try with sincLength = 1 - we should see only two nonzero samples
  // i think, the delayline is overwritten before its sample is consumed






  int numSpikes = 0;
  int nextSpike = 0;
  double ts = 0.0;  // exact time of spike
  double tf = 0.0;  // fractional part of ts
  rsSetZero(y);
  sbl.reset();
  for(int n = 0; n < N; n++) {
    if(n == nextSpike) {
      x[n] = 1;

      sbl.addImpulse(1-tf, 1);
      y[n] = sbl.getSample(1);
      numSpikes++;

      ts += period;
      nextSpike = (int) ceil(numSpikes*period);
      tf = ts - (nextSpike-1);
    }
    else {
      x[n] = 0;
      y[n] = sbl.getSample(0);
    }
  }


  // compensate delay for easier comparison:
  rsArray::shift(&y[0], N, -sbl.getDelay());

  //rsPlotVector(x);
  rsPlotVectors(x, y);


  //GNUPlotter plt;
}

void particleForceDistanceLaw()
{
  // Plots the force vs the distance of the rsPartcielSystem class for various choices of the
  // parameters.

  rsParticleSystemF ps(1);
  ps.setForceLawExponent(2);
  //ps.setForceLawOffset(0.1);

  static const int N = 1000;
  static const int numExponents = 7;
  float exponents[numExponents] = { -3, -2, -1, 0, 1, 2, 3 };
  float size1 = 0.01f;
  float size2 = 0.01f;

  float test = ps.getForceByDistance(1, size1, size2);

  float dMin = 0;
  float dMax = 2;
  float d[N];
  float f[numExponents][N];
  rsArray::fillWithRangeLinear(d, N, dMin, dMax);
  for(int p = 0; p < numExponents; p++)
  {
    ps.setForceLawExponent(exponents[p]);
    for(int n = 0; n < N; n++)
      f[p][n] = ps.getForceByDistance(d[n], size1, size2);
  }

  GNUPlotter plt;
  plt.addDataArrays(N, d, f[0], f[1], f[2], f[3], f[4], f[5], f[6]);
  plt.plot();

  // hmm...somehow the force-distance laws currently implemented are not really good. Ideally, we 
  // would like the law to satisfy:
  // f(d=0) = y0 (y0: user parameter)
  // f(d=1) = y1 (a1: another user parameter)
  // f(d->inf) -> f^p (p: power, also user parameter, asymptote should work also for p < 0)
  // for y0 = inf, y1 = 1, p = -2, the physically correct (for gravitaion, etc) inverse-square law
  // should come out

  // try: f(d) = (a+d)^p, 1/(a + d^-p), a / (b + d^-p)

  // or maybe the function f(d) = d / (c + d^-p) that we currently effectively use is not so bad
  // it may correspond to cloud-like particles that exert no force on each other when they are at
  // the same location

  // for wolfram: derivative of x / (c + x^p) with respect to x
  // finds also the location where the derivative is 0, evaluating f(x) at that location gives the
  // maximum, solving for c gives c in terms of a desired force-limit

  // maybe give each particle a size and then use f = d / (size1 + size2 + d^p) - that corresponds 
  // to a mental model where each particle is a spherically symmetric cloud of matter, the clouds
  // can peneterate and fall through each other
  // for the force of a uniform cloud, see section 5.7 here:
  // http://www.feynmanlectures.caltech.edu/II_05.html
  // the function above is qualitatively similar bur more smooth, i think (verify)

}

void getTwoParticleTrajectories(rsParticleSystemF& ps, int N, float* x1, float* y1, float* z1,
  float* x2, float* y2, float* z2, float* Ek, float* Ep, float* Et, float* Eg, float* Ee, 
  float* Em)
{
  ps.reset();
  for(int n = 0; n < N; n++)
  {
    x1[n] = ps.particles[0].pos.x;
    y1[n] = ps.particles[0].pos.y;
    z1[n] = ps.particles[0].pos.z;

    x2[n] = ps.particles[1].pos.x;
    y2[n] = ps.particles[1].pos.y;
    z2[n] = ps.particles[1].pos.z;

    Ek[n] = ps.getKineticEnergy();
    Ep[n] = ps.getPotentialEnergy();
    Et[n] = ps.getTotalEnergy();

    Eg[n] = ps.getGravitationalPotentialEnergy();
    Ee[n] = ps.getElectricPotentialEnergy();
    Em[n] = ps.getMagneticPotentialEnergy();

    ps.updateState();
  }
}
void particleSystem()
{
  // We simulate a simple system of two particles with unit mass and unit charge to see, if they
  // behave as physically expected (i.e. to see, if the force euqations ae plausible)

  static const int N = 2000; // number of steps in the simulations
  float stepSize = 0.0002f;

  // create and set up the particle system:
  rsParticleSystemF ps(2);

  // both particles have unit mass, charge and size:
  ps.particles[0].mass   = 1.f;
  ps.particles[0].charge = 1.f;
  ps.particles[0].size   = 1.f;
  ps.particles[1].mass   = 1.f;
  ps.particles[1].charge = 1.f;
  ps.particles[1].size   = 1.f;


  // place them at (-1,0,0) and (+1,0,0) with zero velocity initially:
  ps.initialPositions[0]  = rsVector3DF(-0.5f, +0.0f,  +0.0f);
  ps.initialPositions[1]  = rsVector3DF(+0.5f, -0.0f,  +0.0f);
  ps.initialVelocities[0] = rsVector3DF( 0.0f, -0.01f, -0.0f);
  ps.initialVelocities[1] = rsVector3DF( 0.0f, +0.01f, +0.0f);

  // in a first run, we only let them interact via gravitation - they should attract each other:
  ps.setGravitationalConstant(1.0f);
  ps.setElectricConstant(0.0f);
  ps.setMagneticConstant(0.2f);
  ps.setStepSize(stepSize);
  ps.setForceLawExponent(2.0f); // 2: inverse-square law (asymptotic), 0: force distance-independent
  //ps.setForceLawOffset(1.0);   // 0: non-asymptotic inverse power law

  // record trajectories and energies:
  float x1[N], y1[N], z1[N], x2[N], y2[N], z2[N];  // coordinates
  float Ek[N], Ep[N], Et[N];     // kinetic, potential, total energy
  float Eg[N], Ee[N], Em[N];     // gravitational, electric and magnetic potential energy


  getTwoParticleTrajectories(ps, N, x1, y1, z1, x2, y2, z2, Ek, Ep, Et, Eg, Ee, Em);

  // they initially approach each other, fly through each other and then drift apart forever
  // i suppose, this is due to the singularity, when they are very close (division by zero
  // -> infinite force)


  // test cross-product formula - move to unit-tests:
  rsVector3DF test = cross(rsVector3DF(2,3,5), rsVector3DF(7,11,13)); 
  // bool result = test == rsVector3DF(-16,9,1); // ok - cross-product is correct


  GNUPlotter plt;
  float t[N];
  createTimeAxis(N, t, 1/stepSize);
  plt.addDataArrays(N, t, x1, y1, z1, x2, y2, z2);
  //plt.addDataArrays(N, t, x1, x2);
  //plt.addDataArrays(N, t, y1, y2);
  //plt.addDataArrays(N, t, z1, z2);
  //plt.addDataArrays(N, t, Et, Ek, Ep);
  //plt.addDataArrays(N, t, Et, Ek, Eg, Ee, Em);
  //plt.addDataArrays(N, t, Et, Ek, Eg);
  //plt.addDataArrays(N, t, Ek);
  plt.plot();

  // Observations:
  // -when starting at p1=(-1,0.0),p2=(+1,0,0),v1=v2=(0,0.2,0) with only magnetic forces active,
  //  i would expect the x-coordinates converge (initially) and the y-coordinates just linearly 
  //  increase for both partcile (and the z-coordinates stay zero) - but something weird happens
  //  -either the magnetic force computaion is wrong, or the cross-product formula or my 
  //   expectation is wrong
  //  -it seems to depend on the stepsize
  //  -maybe we should apply the stepSize to the velocity update?
}

void bouncillator()
{
  rsBouncillatorF rb;
  rb.setFloor(    -0.0f);
  rb.setCeil(     +1.0f);
  rb.setIncrement( 0.05f);
  rb.setDecrement( 0.02f);
  rb.setShape(    +0.016f);

  // create output sequence:
  static const int N = 500;   // number of output samples
  float x[N];
  rb.resetToMin();
  for(int n = 0; n < N; n++)
    x[n] = rb.getSample();

  GNUPlotter plt;
  plt.addDataArrays(N, x);
  plt.plot();
}

void bouncillatorFormula()
{
  // test, whether the explicit formula for the particle poisition x[n] gives correct results

  float min   = -0.2f;
  float max   = +1.5f;
  float inc   =  0.04f;
  float dec   =  0.04f;
  float shape =  0.01f;
  //float start = -0.2f;

  rsBouncillatorF bnc;
  bnc.setFloor(    min);
  bnc.setCeil(     max);
  bnc.setIncrement(inc);
  bnc.setDecrement(dec);
  bnc.setShape(shape);


  // create output sequence using the bouncillator object:
  static const int N = 60;   // number of output samples
  float x[N], xp[N];          // computed and predicted output
  bnc.resetToMin();
  //bnc.resetToMax();
  for(int n = 0; n < N; n++)
  {
    x[n]  = bnc.getSample();
    xp[n] = bnc.predictOutput(float(n), min, inc, 1+shape);
    //xp[n] = bnc.predictOutput(float(n), max, -dec, 1+shape);
  }

  // predict instant of hitting the wall:
  float nw = bnc.getInstantForHitting(max, min, inc, 1+shape);
  //float nw = bnc.getInstantForHitting(min, max, -dec, 1+shape);

  GNUPlotter plt;
  plt.addDataArrays(N, x);
  plt.addDataArrays(N, xp);
  plt.plot();
}

// test inputs for FM/PM to replace sin
double saw2(double phi)
{
  return saw(phi, 2);
}
double sqr3(double phi)
{
  return sqr(phi, 3);
}
void constPowerFade(double x, double* s1, double* s2)
{
  double p = x*0.5*PI;
  *s1 = cos(p);
  *s2 = sin(p);
}
void freqVsPhaseMod()
{
  // Compares frequency modulation with phase modulation, also creates a mixed-modulation signal
  // with adjustable morphing factor (0..1) and tries to transform FM into PM and vice versa by 
  // integration/differentiation of the modulator signal
  // use as carrier:    c(t) = sin(wc*t) + sin(2*wc*t)/2
  // and as modulator:  m(t) = sin(wm*t) + sin(3*wm*t)/3
  // or vice versa - each has the first two sine components of a saw or square wave - they 
  // waveforms are still simple but maybe complex enough to expose the difference between FM and PM


  // experiment parameters:
  static const int N = 882;     // number of samples to produce
  //static const int N = 2000;      // number of samples to produce
  double fs = 44100;              // sample rate
  double fc = 300;                // carrier freq
  double fm = 200;                // modulator freq
  double depth = 1.0;             // modulation depth
  double fmVsPm = 0.5;            // 0: pure FM, 1: pure PM
  double lpCutoff = 10;           // cutoff freq for integrator filter
  double hpCutoff = 10;           // cutoff freq for differentiator filter
  //double (*modWave) (double x) = sqr3;  // modulator waveform function (like sin, cos, ...)
  //double (*carWave) (double x) = saw2;  // carrier waveform function
  double (*modWave) (double x) = sin;
  double (*carWave) (double x) = sin;


  typedef RAPT::rsOnePoleFilter<double, double> Flt;
  Flt lpf, hpf;
  lpf.setMode(Flt::LOWPASS_IIT);  // later use BLT
  lpf.setSampleRate(fs);
  lpf.setCutoff(lpCutoff);
  hpf.setMode(Flt::HIGHPASS_MZT); // later use BLT
  hpf.setSampleRate(fs);
  lpf.setCutoff(hpCutoff);
  // 1st order lowpass and highpass filters have a phase of 45° at the cutoff frequency but we
  // actually need 90° to turn a sine into a (positive or negative) cosine
  // an allpass filter, on the other hand, has 90° phase shift at the cutoff, so just putting an
  // allpass after a lowpass give too much phase shift
  // maybe a lowpass and an allpass with a different cutoff freq could make sense
  // orr...no...the integrator should actually have a cutoff of zero or almost zero


  // create signals:
  double t[N], yC[N], yM[N];      // time axis, unmodulated carrier, modulator
  double yFM[N], yPM[N], yMM[N];  // FM output, PM output, mixed-mod output
  double yMd[N], yMi[N];          // differentiated and integrated modulator signal
  //double yPFM, yPPM;              // pseudo-FM and pseudo-PM outputs
  double wc  = 2*PI*fc/fs;
  double wm  = 2*PI*fm/fs;
  //double pm  = 0;   // instantaneous phase of modulator
  double pFM = 0;   // instantaneous phase of FM signal
  double pPM = 0;   // instantaneous phase of PM signal (before phase-modulation)
  double pMM = 0;   // instantaneous phase of mixed modulation signal
  double scaleFM, scalePM;
  constPowerFade(fmVsPm, &scaleFM, &scalePM);
  for(int n = 0; n < N; n++)
  {
    // time axis, unmodulated carrier and modulator:
    double tn = n/fs;    // current time instant (in seconds?)
    t[n]   = tn;
    yM[n] = modWave(wm*n);  // modulator output (needed in equations below)
    yC[n] = carWave(wc*n);  // unmodulated carrier output (needed just for the plot)

    // frequency modulation:
    double wInst = wc + depth * yM[n] * wc; // instantaneous omega
    yFM[n] = carWave(pFM);
    pFM += wInst;

    // phase modulation:
    double pInst = pPM + depth * yM[n];     // instantaneous phase
    yPM[n] = carWave(pInst);
    pPM += wc;
    // do we need a constant scale-factor, like 2*PI, to scale the phase-offset? maybe compare
    // spectra of FM and PM signals with the same index - they should have the same bandwidth..or
    // maybe even equal in terms of magnitudes?

    // mixed modulation:
    wInst = wc  + scaleFM * depth * yM[n] * wc;   // instantaneous omega
    pInst = pMM + scalePM * depth * yM[n];        // instantaneous phase
    yMM[n] = carWave(pInst);
    pMM += wInst;



    // try to obtain pseudo-FM via PM and pseudo PM via FM by using an integrated or 
    // differentiated modulator signal
    yMd[n] = hpf.getSample(yM[n]);
    yMi[n] = lpf.getSample(yM[n]);
    // we need gain factors equal to the reciprocal of the magnitude response of the filters
    // at the modulator frequency


  }

  // todo: maybe to analyze the spectra, choose frequencies such that the period of the signals is 
  // a power of two, suitable for FFT (or use arbitrary length FFT, tuned to actual cycle length)
  // i think, the period is given by the lowest common multiple of the carrier and modulator freq

  // plot:
  GNUPlotter plt;
  //plt.addDataArrays(N, t, yC,  yM);
  plt.addDataArrays(N, t, yM, yMd, yMi); // modulator, differentiated mod, integrated mod
  //plt.addDataArrays(N, t, yFM, yPM, yMM);
  plt.plot();
}

void rayBouncer()
{
  rsRayBouncerF rb;
  rb.setEllipseParameters(1, 2, float(PI/4), 0.1f, 0.2f);
  rb.setInitialPosition(0.2f, -0.4f);
  rb.setSpeed(0.02f);
  rb.setLaunchAngle(float(PI) * 30.f / 180.f);

  // create output sequence:
  static const int N = 8000;   // number of output samples
  float x[N], y[N];
  rb.reset();
  for(int n = 0; n < N; n++)
    rb.getSampleFrame(x[n], y[n]);

  GNUPlotter plt;
  plt.setRange(-1.5, +1.5, -1.5, +1.5);
  plt.setPixelSize(600, 600);
  plt.addCommand("set size square");           // set aspect ratio to 1:1 ..encapsulate in GNUPlotter
  plt.addDataArrays(N, x, y);
  //plt.addDataArrays(Ne, xe, ye);
  plt.plot();
}

void circleFractals()
{
  // http://benice-equation.blogspot.de/2012/01/fractal-spirograph.html
  // from the comments:
  // Speed(n) = k^(n-1), k=±2, ±3, ±4, ....
  // The most general parametric equation is:
  // x(t) = ?(k=1 to n) R(k)*cos(a(k)*t)
  // y(t) = ?(k=1 to n) R(k)*sin(a(k)*t)

  // R(k): radius of k-th circle, a(k): angular velocity?

  // hmm - does not yet really work well (pictures not as interesting as on benice-equation)
  // ...more research needed

  int N = 501; // number of points
  std::vector<double> x(N), y(N);


  std::vector<double> r, w; 
  //r = { 1, .5, .25, .125 };  // radii
  //w = { 1,  2,   4, 8    };  // relative frequencies

  //r = { 1, .5};  // cardioid
  //w = { 1,  2};

  r = { 1, .7, .5};
  w = { 1, 2,  4};

  //r = { 1,  1/3.}; // nephroid
  //w = { 1,  3};

  //r = { 1,  1/4.}; // 
  //w = { 1,  4};

  //r = { 1, .7, .5, .4, .3, .2};
  //w = { 1, 5,  10, 15, 20, 25};

  r = { 1, .5, .25,};
  w = { 1,  5,  25};

  for(size_t n = 0; n < N; n++)
  {
    x[n] = y[n] = 0;
    for(size_t k = 0; k < r.size(); k++)
    {
      double wk = 2*PI*w[k] / (N-1); // N for open loop, N-1 for closed loop
      x[n] += r[k] * cos(wk*n);
      y[n] += r[k] * sin(wk*n);
    }
  }

  GNUPlotter plt;
  plt.addDataArrays(N, &x[0], &y[0]);
  //plt.setRange(-2, +2, -2, +2);
  plt.setPixelSize(400, 400);
  plt.addCommand("set size square");  // set aspect ratio to 1:1
  plt.plot();
}

// from https://en.wikipedia.org/wiki/Hilbert_curve
void rot(int n, int *x, int *y, int rx, int ry) //rotate/flip a quadrant appropriately
{
  if (ry == 0) {
    if (rx == 1) {
      *x = n-1 - *x;
      *y = n-1 - *y;
    }
    //Swap x and y
    int t  = *x;
    *x = *y;
    *y = t;
  }
}
int xy2d (int n, int x, int y)           // convert (x,y) to d
{
  int rx, ry, s, d=0;
  for (s=n/2; s>0; s/=2) {
    rx = (x & s) > 0;
    ry = (y & s) > 0;
    d += s * s * ((3 * rx) ^ ry);
    rot(s, &x, &y, rx, ry);
  }
  return d;
}
void d2xy(int n, int d, int *x, int *y)  // convert d to (x,y)
{
  int rx, ry, s, t=d;
  *x = *y = 0;
  for (s=1; s<n; s*=2) {
    rx = 1 & (t/2);
    ry = 1 & (t ^ rx);
    rot(s, x, y, rx, ry);    
    *x += s * rx;
    *y += s * ry;
    t /= 4;
  }
}
void hilbertCurve()
{
  int order = 3;
  int n = (int) pow(4, order); // wiki says, it could be a power of 2 (but then it'a not a square)
  int max = (int) pow(2, order);

  std::vector<int> x(n), y(n);
  int i;
  for(i = 0; i < n; i++)
    d2xy(n, i, &x[i], &y[i]);

  GNUPlotter plt;
  plt.addDataArrays(n, &x[0], &y[0]);
  plt.setRange(-1, max+1, -1, max+1);
  plt.setPixelSize(400, 400);
  plt.addCommand("set size square");  // set aspect ratio to 1:1
  plt.plot();

  // see:
  // https://www.reddit.com/r/visualizedmath/comments/7xtxgb/hilbert_curve/
  // http://www4.ncsu.edu/~njrose/pdfFiles/HilbertCurve.pdf
  // http://www.fundza.com/algorithmic/space_filling/hilbert/basics/index.html
  // https://en.wikipedia.org/wiki/Hilbert_curve
  // https://en.wikipedia.org/wiki/Moore_curve
  // https://github.com/adishavit/hilbert
  // https://marcin-chwedczuk.github.io/iterative-algorithm-for-drawing-hilbert-curve
  // http://wwwmayr.informatik.tu-muenchen.de/konferenzen/Jass05/courses/2/Valgaerts/Valgaerts_paper.pdf
  // http://people.csail.mit.edu/jaffer/Geometry/HSFC
  // https://arxiv.org/pdf/1211.0175.pdf
  // http://www.dcs.bbk.ac.uk/~jkl/thesis.pdf
  // http://www.fractalcurves.com/Taxonomy.html
  // http://archive.org/stream/BrainfillingCurves-AFractalBestiary/BrainFilling#page/n27/mode/2up
  // http://www.diss.fu-berlin.de/diss/servlets/MCRFileNodeServlet/FUDISS_derivate_000000003494/2_kap2.pdf?hosts=

  // maybe i should implement a Lindenmayer system ("L-system")
}

void lindenmayer()
{
  // Uses a Lindenmayer system to produce various 2D curves.
  
  std::vector<double> x, y;
  rosic::LindenmayerRenderer lr;

  // uncomment the curve, you want to render:
  //lr.getKochSnowflake(4, x, y);
  //lr.getMooreCurve(4, x, y);
  //lr.get32SegmentCurve(2, x, y);
  //lr.getQuadraticKochIsland(3, x, y);
  //lr.getSquareCurve(4, x, y);
  //lr.getSierpinskiTriangle(6, x, y);

  //lr.getSierpinskiTriangle2(2, x, y); // doesn't work
  //lr.getPleasantError(3, x, y);

  // some of my own experiments - get rid of passing the angle to each call of render, call
  // setAngle once instead (after setting the seed)

  // shapes based on a triangle seed:
  std::string seed = "F+F+F";
  //lr.clearRules(); lr.addRule('F', "F+F-F"); lr.render(seed, 7, 120, x, y); // sort of triangular grid
  //lr.clearRules(); lr.addRule('F', "F+F-FF"); lr.render(seed, 6, 120, x, y); // thin spiral arms
  //lr.clearRules(); lr.addRule('F', "F+F-F+"); lr.render(seed, 8, 120, x, y); 
  //lr.clearRules(); lr.addRule('F', "F+F-F-F"); lr.render(seed, 5, 120, x, y); 

  // shapes based on a square seed:
  seed = "F+F+F+F";
  //lr.clearRules(); lr.addRule('F', "F+F-F-FF+F+F-F"); lr.render(seed, 3, 90, x, y);
  //lr.clearRules(); lr.addRule('F', "F+FF-FF-FF+F+F-F"); lr.render(seed, 3, 90, x, y);
  //lr.clearRules(); lr.addRule('F', "F+FF-F-FF+F+F-F"); lr.render(seed, 4, 90, x, y);
  //lr.clearRules(); lr.addRule('F', "F+FF-FF-FF+FF+FF-F"); lr.render(seed, 3, 90, x, y);
  //lr.clearRules(); lr.addRule('F', "F+F-F-FF+FF+FF-F"); lr.render(seed, 3, 90, x, y);
  //lr.clearRules(); lr.addRule('F', "F+FF-FF-FFFF+F+FF-F"); lr.render(seed, 3, 90, x, y);
  //lr.clearRules(); lr.addRule('F', "F+FF-FF-FFF+F+F-F"); lr.render(seed, 4, 90, x, y);   // nice
  //lr.clearRules(); lr.addRule('F', "F+FF-FF-FfF+F+F-F"); lr.render(seed, 4, 90, x, y);
  lr.clearRules(); lr.addRule('F', "f+FF-FF-FFF+F+F-F"); lr.render(seed, 4, 90, x, y); 
  //lr.clearRules(); lr.addRule('F', "F+FF-FF-FFF+f+F-F"); lr.render(seed, 4, 90, x, y);
  //lr.clearRules(); lr.addRule('F', "F+FF-FF-FFF+F+f-F"); lr.render(seed, 4, 90, x, y); // swastika?
  //lr.clearRules(); lr.addRule('F', "F+FF-FF-FFF+F+F-f"); lr.render(seed, 4, 90, x, y);

  // shapes based on sort of s-shape seed:
  seed = "F+F-FF+FF+F+F-FF+FF";
  //lr.clearRules(); lr.addRule('F', "F+F-F-FF+F+F-F"); lr.render(seed, 3, 90, x, y); 
  //lr.clearRules(); lr.addRule('F', "F+f-F-FF+f+F-F"); lr.render(seed, 4, 90, x, y); 
  //lr.clearRules(); lr.addRule('F', "F+F-F"); lr.render(seed, 5, 90, x, y); 
  //lr.clearRules(); lr.addRule('F', "+F-f-F+"); lr.render(seed, 8, 90, x, y); // digital noise

  // shapes based on an octagon seed:
  seed = "F+F+F+F+F+F+F+F";
  //lr.clearRules(); lr.addRule('F', "F+F---FF+++F-"); lr.render(seed, 5, 45, x, y);
  //lr.clearRules(); lr.addRule('F', "+F---FF+++F-"); lr.render(seed, 5, 45, x, y);
  //lr.clearRules(); lr.addRule('F', "F+F---FF+++F-"); lr.render(seed, 5, 45, x, y);
  //lr.clearRules(); lr.addRule('F', "+F--F"); lr.render(seed, 8, 45, x, y);
  //lr.clearRules(); lr.addRule('F', "+F--F+"); lr.render(seed, 7, 45, x, y); // flowerish
  //lr.clearRules(); lr.addRule('F', "+FF--FF+"); lr.render(seed, 5, 45, x, y); 

  // attempt to do a sierpinski triangle (using 60°) - not very successfuly yet:
  seed = "F++F++F";
  //lr.clearRules(); lr.addRule('F', "+F-F+F-F"); lr.render(seed, 5, 60, x, y); // nope, but cool hexgonal shape
  //lr.clearRules(); lr.addRule('F', "+F-F+F-F+"); lr.render(seed, 6, 60, x, y);
  //lr.clearRules(); lr.addRule('F', "+F--F+F--F"); lr.render(seed, 5, 60, x, y); 
  //lr.clearRules(); lr.addRule('F', "+F--F++F--F"); lr.render(seed, 2, 60, x, y); 
  //lr.clearRules(); lr.addRule('F', "+F--F++F--F--FF"); lr.render(seed, 4, 60, x, y); 

  seed = "FX++FX++FX";
  //lr.clearRules(); lr.addRule('X', "+FX++FX++FX+"); lr.render(seed, 4, 60, x, y); 
  //lr.clearRules(); lr.addRule('F', "X"); lr.addRule('X', "FX++FX++FX"); lr.render(seed, 4, 60, x, y); 

  // strategies for variation: replace a single F with f in a good looking pic




  // plot 2D:
  GNUPlotter plt2;
  plt2.addDataArrays((int)x.size(), &x[0], &y[0]);
  plt2.setRange(-1.1, +1.1, -1.1, 1.1);
  plt2.setPixelSize(600, 600);
  plt2.addCommand("set size square");  // set aspect ratio to 1:1
  plt2.plot();

  // plot 1D:
  GNUPlotter plt1;
  plt1.addDataArrays((int)x.size(), &x[0]);
  plt1.addDataArrays((int)y.size(), &y[0]);
  plt1.plot();
}

void triSawOsc()
{
  static const int N = 1000;   // number of output samples
  float T = 250;               // period in samples

  // set up osc:
  RAPT::rsTriSawOscillator<float> osc;
  osc.setPhaseIncrement(1.f/T);
  osc.setAsymmetry(-0.5);

  // test:
  float a = osc.asymForTransitionSamples(20);
  osc.setAsymmetry(a);

  // generate signal:
  float y[N];
  for(int n = 0; n < N; n++)
    y[n] = osc.getSample();

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, y);
  plt.plot();
}

void triSawOscAntiAlias()
{
  // We create an exponential sweep with the TriSaw oscillator and write it into a wave file for 
  // investigating the aliasing behavior. As a crude and simple form of anti-aliasing, we limit the
  // asymmetry parameter in a way to achieve some minimum absolute number of samples for the 
  // transition.

  double fs = 44100;   // sample rate
  double fL = 500;     // lower frequency
  double fU = fs/2;    // upper frequency
  double length = 4;   // length in seconds
  double asym   =  1;  // asymmetry
  double sigmo  =  0;  // sigmoidity
  double amp    = 0.8; // amplitude

  int N = (int) ceil(length*fs);

  std::vector<double> f(N), x(N), xa(N);
  RAPT::rsArray::fillWithRangeExponential(&f[0], N, fL, fU);

  RAPT::rsTriSawOscillator<double> osc;
  osc.setAsymmetry(asym);
  osc.setAttackSigmoid(sigmo);
  int n;

  for(n = 0; n < N; n++) {
    osc.setPhaseIncrement(f[n]/fs);
    x[n] = amp * osc.getSample();
  }
  rosic::writeToMonoWaveFile("TriSawSweep-SawNoAA.wav", &x[0], N, 44100, 16);

  double trans  =  2;  // minimum number of samples for transition
  double tmp;
  osc.reset();
  for(n = 0; n < N; n++) {
    osc.setPhaseIncrement(f[n]/fs);
    tmp = osc.asymForTransitionSamples(trans);
    tmp = RAPT::rsClip(asym, -tmp, tmp);
    osc.setAsymmetry(tmp);
    x[n] = amp * osc.getSample();
  }
  rosic::writeToMonoWaveFile("TriSawSweep-SawLimitedAsym.wav", &x[0], N, 44100, 16);

  // maybe the minimum number of transition samples should itself be a function of frequency (i.e.
  // in create with frequency). For low frequencies, there should be no limit in order to not have 
  // too much lowpass character in the signal...
  double k1 = 8, k2 = 0;
  double inc;
  osc.reset();
  for(n = 0; n < N; n++) {
    inc   = f[n]/fs;
    trans = k1*inc + k2*inc*inc;
    osc.setPhaseIncrement(inc);
    tmp = osc.asymForTransitionSamples(trans);
    tmp = RAPT::rsClip(asym, -tmp, tmp);
    osc.setAsymmetry(tmp);
    x[n] = amp * osc.getSample();
  }
  rosic::writeToMonoWaveFile("TriSawSweep-SawDynamicallyLimitedAsym.wav", &x[0], N, 44100, 16);
  // hmm...no the function trans(inc) = k1*inc + k2*inc^2 seems not suitable. in the low frequency 
  // range there's not enough limiting - we need a function that rises faster intitially - maybe 
  // once again, the rational map could be used.
  // ...maybe experiment with that in ToolChain
  // inc = k1 * ratMap(2*inc, k2); 
  // use 2*inc because at fs/2 the increment is 0.5 - at that freq, we want the limit to be k1


  // to do something similar with the sigmoidity, we would need to set a lower limit for the 
  // sigmoidity
  // sigmo = rsMax(sigmoSetting, minSigmo)
  // where minSigmo should be some function of the increment, we should probably want:
  // minSigmo(inc=0.25) = 1 such that a triangle turns into a sine when inc=0.25
  // minSigmo = rsMin(1, k0 + k1*inc)
  // k0 should probably be = -2, so we may reach sigmo = -2 at inc = 0 (which is the lowest 
  // meaningful setting)
  // OK - let's create a sweep with 0 asymmetry and full negative sigmoidity for reference:
  sigmo = -2;
  osc.setAsymmetry(0);
  osc.setAttackSigmoid(-2);
  osc.setDecaySigmoid( -2);
  osc.reset();
  for(n = 0; n < N; n++) {
    osc.setPhaseIncrement(f[n]/fs);
    x[n] = amp * osc.getSample();
  }
  rosic::writeToMonoWaveFile("TriSawSweep-AntiSigmoidNoAA.wav", &x[0], N, 44100, 16);

  // now with sigmoidity limiting:
  k1 = 20; // with k1=12 at inc=0.25, minSigmo becomes 1 (12*0.25=3), 20 sounds good
  double minSigmo;
  osc.reset();
  for(n = 0; n < N; n++) {
    inc   = f[n]/fs;
    minSigmo = RAPT::rsMin(1.0, -2 + k1*inc);
    tmp = rsMax(sigmo, minSigmo);
    osc.setAttackSigmoid(tmp);
    osc.setDecaySigmoid( tmp);
    osc.setPhaseIncrement(inc);
    x[n] = amp * osc.getSample();
  }
  rosic::writeToMonoWaveFile("TriSawSweep-AntiSigmoidLimited.wav", &x[0], N, 44100, 16);
  // 20 sounds good, but it probably also be a function of asymmetry setting -> more asymmetry 
  // should give more limiting...maybe use minSigmo = -2 + (k1+k2*|asym|)*inc where asym is the
  // already limited asymmetry

}

void xoxosOsc()
{
  // Oscillator based on an ellipse in the xy-plane
  // more info:
  // https://www.kvraudio.com/forum/viewtopic.php?p=6656238#p6656238
  // https://gitlab.com/Hickler/Soundemote/issues/67
  // https://github.com/RobinSchmidt/RS-MET/issues/72
  // play with parameters:
  // https://www.desmos.com/calculator/7h9mknbv3q
  // i think, it works as follows:
  // -create x,y values on a circle (standard rotating phasor in the plane)
  // -convert x,y to values on an arbitrary ellipse
  //  (-project onto the x- and y-axis (i.e. take the x- and y-value))...really?


  static const int N = 1000;   // number of output samples
  float T = 250;               // period in samples
  float A = 0.8f;              // -1..+1,      Upper/Lower Bias (?)
  float B = 0.5f;              // -PI..+PI     Left/Right Bias (?)
  float C = 0.5f;              // -inf..+inf   Upper/Lower Pinch (?)
  float w = float(2*PI/T);

  rosic::rsEllipseOscillator osc;
  osc.setA(A);
  osc.setB(B);
  osc.setC(C);
  osc.setOmega(w);

  // generate signals:
  float x[N], y[N], sum[N];
  float sB = sin(B);
  float cB = cos(B);
  //float s, c, Ac, Cs, a;
  for(int n = 0; n < N; n++)
  {
    double xt, yt;
    osc.getSamplePair(&xt, &yt);
    x[n] = (float) xt;
    y[n] = (float) yt;
    sum[n] = x[n] + y[n];

    /*
    s  = sin(w*n);
    c  = cos(w*n);
    Ac = A + c;
    Cs = C * s;
    a  = 1 / sqrt(Ac*Ac + Cs*Cs);  // normalizer
    x[n]   = a*Ac*cB;
    y[n]   = a*Cs*sB;
    sum[n] = x[n] + y[n]; // maybe it could be further flexibilized by taking a weighted sum?
    */
  }

  // plot outputs:
  GNUPlotter plt;
  //plt.addDataArrays(N, x, y); // plot the ellipse
  plt.addDataArrays(N, sum);
  plt.addDataArrays(N, x);
  plt.addDataArrays(N, y);
  plt.plot();

  // Observations:
  // B=C=0: pulse wave, A controls pulse-width
  // A=1, B=PI/2, C = +-0.8: saw wave, higher C makes long transition more sigmoid

  // todo: check the parameter meanings, i think, they are y-offset, x-scale and rotation-angle

  // make a 3D version of it:
  // -obtain new point x,y,z via 3D rotation matrix (this has 3 frequencies)
  // -transform point (which is on sphere) to point on ellipsoid (which contains origin)
  // -project onto unit sphere by dividing x,y,z by  sqrt(x^2 + y^2 + z^2)
  // -take linear combination of x,y,z as output
}
