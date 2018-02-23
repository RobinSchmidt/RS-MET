#include "GeneratorExperiments.h"
using namespace RAPT;
//using namespace rosic;

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

void lindenmayerKoch()
{
  // Uses a Lindenmayer system to produce a Koch snowflake.
  // https://en.wikipedia.org/wiki/Koch_snowflake#Representation_as_Lindenmayer_system

  int order = 3;

  // set up and run the L-system:
  LindenmayerSystem ls;
  ls.addRule('F', "F+F--F+F");
  std::string result = ls.apply("F--F--F", order);

  // translate into curve:
  TurtleGraphics tg;
  tg.setAngle(60);
  tg.init(0, 0, 1, 0);
  std::vector<double> x, y;
  tg.translate(result, x, y);

  // plot:
  GNUPlotter plt;

  //// 1D:
  //plt.addDataArrays(x.size(), &x[0]);
  //plt.addDataArrays(y.size(), &y[0]);

  // 2D:
  plt.addDataArrays((int)x.size(), &x[0], &y[0]);
  plt.setPixelSize(400, 400);
  plt.addCommand("set size square");  // set aspect ratio to 1:1

  plt.plot();

  // todo: factor out a LindenmayerRenderer class that encapsulates creation of the string, 
  // translation into a sequence of points and have functions for generating some standard curves, 
  // like: ls.getKochSnowFlake(order, x, y), getMooreCurve(order, x, y), etc.
}

void lindenmayerMoore()
{
  // Uses a Lindenmayer system to produce a Moore curve.
  // see: https://en.wikipedia.org/wiki/Moore_curve

  int order = 4; // order of the Moore curve (for audio, up to 5 makes sense)
  int max = (int) pow(2, order+1);

  // set up the L-system:
  LindenmayerSystem ls;
  ls.addRule('L', "-RF+LFL+FR-");
  ls.addRule('R', "+LF-RFR-FL+");
  std::string seed = "LFL+F+LFL"; // original Moore curve

  // iterate the L-system:
  std::string result = ls.apply(seed, order);

  // translate the result string into a Moore curve:
  TurtleGraphics tg;
  tg.init(0, max/2-1, 1, 0);
  std::vector<double> x, y;
  tg.translate(result, x, y);

  // plot:
  GNUPlotter plt;

  //// 1D:
  //plt.addDataArrays(x.size(), &x[0]);
  //plt.addDataArrays(y.size(), &y[0]);

  // 2D:
  plt.addDataArrays((int)x.size(), &x[0], &y[0]);
  plt.setRange(-1, max, -1, max);
  plt.setPixelSize(400, 400);
  plt.addCommand("set size square");  // set aspect ratio to 1:1

  plt.plot();

  // other closed curves that can be generated:
  // http://mathforum.org/advanced/robertd/lsys2d.html (many curves with L-system rules)
  // http://www.kevs3d.co.uk/dev/lsystems/ (applet with examples)

  // https://www.cut-the-knot.org/do_you_know/hilbert.shtml

  // not sure, if doable by L-system:
  // https://en.wikipedia.org/wiki/Sierpi%C5%84ski_curve
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
  float s, c, Ac, Cs, a;
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
