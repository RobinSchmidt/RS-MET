#include "Experiments.h"

#include "PrototypesForAPE.cpp"

#include "../rapt_tests/Experiments/ScratchPad.cpp"
#include "../rapt_tests/Experiments/MathExperiments.cpp"
#include "../rapt_tests/Experiments/FilterExperiments.cpp"
#include "../rapt_tests/Experiments/GeneratorExperiments.cpp"
#include "../rapt_tests/Experiments/ModulatorExperiments.cpp"
#include "../rapt_tests/Experiments/GraphicsExperiments.cpp"

// get rid:
using namespace RAPT;
using namespace rosic;

#include "../rosic_tests/PortedFromRSLib/Experiments/MathExperiments.cpp"
#include "../rosic_tests/PortedFromRSLib/Experiments/AnalysisExperiments.cpp"
#include "../rosic_tests/PortedFromRSLib/Experiments/CellularExperiments.cpp"
#include "../rosic_tests/PortedFromRSLib/Experiments/DelayExperiments.cpp"
#include "../rosic_tests/PortedFromRSLib/Experiments/FilterExperiments.cpp"
#include "../rosic_tests/PortedFromRSLib/Experiments/MiscAudioExperiments.cpp"
#include "../rosic_tests/PortedFromRSLib/Experiments/ModalExperiments.cpp"
#include "../rosic_tests/PortedFromRSLib/Experiments/ModulatorExperiments.cpp"
#include "../rosic_tests/PortedFromRSLib/Experiments/OscillatorExperiments.cpp"
#include "../rosic_tests/PortedFromRSLib/Experiments/PartialExtractionExperiments.cpp"
#include "../rosic_tests/PortedFromRSLib/Experiments/PhaseVocoderExperiments.cpp"
#include "../rosic_tests/PortedFromRSLib/Experiments/PhysicsExperiments.cpp"
#include "../rosic_tests/PortedFromRSLib/Experiments/ResamplingExperiments.cpp"
#include "../rosic_tests/PortedFromRSLib/Experiments/SaturationExperiments.cpp"



void particleBouncerExperiment()
{
  static const int N = 8000;   // number of output samples

  // create and set up particle bouncer:
  ParticleBouncer bouncer;
  bouncer.setEnclosureEllipseAspectRatio(1.15);
  //bouncer.setEnclosureEllipseAspectRatio(1.0);
  //bouncer.setEnclosureEllipseAspectRatio(0.8);
  //bouncer.setInitialPosition(0.2, -0.4);
  bouncer.setInitialPosition(0.1, -0.2);
  //bouncer.setInitialPosition(0.0, 0.0);
  bouncer.setSpeed(0.02);
  bouncer.setLaunchAngle(30.0);
  //bouncer.setLaunchAngle(0.0);

  //// boring, for test:
  //bouncer.setSpeed(0.03248536741);
  //bouncer.setInitialPosition(0.0, 0.0);
  //bouncer.setLaunchAngle(10.0);
  //bouncer.setSpeed(0.02121315726);
  ////bouncer.setLaunchAngle(30.0);


  // create output sequence:
  double x[N], y[N];
  bouncer.reset();
  for(int n = 0; n < N; n++)
    bouncer.getSampleFrame(x[n], y[n]);

  // create enclosing ellipse for reference inside the plot:
  static const int Ne = 100;
  double a = bouncer.getEllipseA();
  double b = bouncer.getEllipseB();
  double xe[Ne], ye[Ne];
  for(int n = 0; n < Ne; n++)
  {
    double t = 2 * PI * double(n) / (Ne-1);
    xe[n] = a * cos(t);
    ye[n] = b * sin(t);
  }
  double areaOverPi = a*b; // actual area is PI*a*b

  // plot sequence:
  GNUPlotter plt;
  plt.setRange(-1.2, +1.2, -1.2, +1.2);
  plt.setPixelSize(600, 600);
  plt.addCommand("set size square");           // set aspect ratio to 1:1 ..encapsulate in GNUPlotter
  plt.addDataArrays(N, x, y);
  plt.addDataArrays(Ne, xe, ye);
  plt.plot();
  //plt.plotFunctionTables(N, x, y);

  // Observations:
  // -when we use a circle of radius 1 and the speed is some value such that the particle hits the
  //  boundary at some sample instant (say speed = 0.02), we seem to run into numerical errors - 
  //  the reflection angles detoriate: try starting at 0.0 with speed 0.02 - it should give a line
  //  but it doesn't. using some more random value value like 0.02121315726, we see the line as
  //  expected
  // -sometimes, the emerging shape is an ellipse, sometimes 2 parabolas
  //  -parabola: ratio=1.15, pos = (0.1, -0.2), angle = 30.0, speed = 0.02
  //  -ellipse:  ratio=1.15, pos = (0.2, -0.4), angle = 30.0, speed = 0.02
  //  -maybe there's a pair of hyperbolas too, for some other values?
}