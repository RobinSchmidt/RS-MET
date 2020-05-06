#include "PhysicsExperiments.h"

void doublePendulum()
{
  int freqDivider = 4;

  int N  = freqDivider*150000;    // number of samples
  int fs = 44100;

  // create and set up the double pendulum object:
  rsDoublePendulumDD dp;
  dp.setLength1(1.0);
  dp.setLength2(0.5);
  dp.setMass1(1.0);
  dp.setMass2(0.5);
  dp.setStepSize(0.01/freqDivider);
  //dp.setStepSize(0.0025);


  // let the pendulum swing and record the locations of both bobs over time:
  vector<double> a1(N), a2(N);                 // angles
  vector<double> m1(N), m2(N);                 // momenta
  vector<double> x1(N), y1(N), x2(N), y2(N);   // x- and y-coordinates
  for(int n = 0; n < N; n++)
  {
    dp.getAngles(&a1[n], &a2[n]);
    dp.getMomenta(&m1[n], &m2[n]);
    dp.getPendulumCoordinates(&x1[n], &y1[n], &x2[n], &y2[n]);
    dp.updateState();
  }

  GNUPlotter plt;
  //plt.addDataArrays(N, &a1[0]);
  //plt.addDataArrays(N, &a2[0]);

  plt.addDataArrays(N, &m1[0]);  // maybe momenta should be scaled by the mass
  plt.addDataArrays(N, &m2[0]);

  //plt.addDataArrays(N, &x1[0]);
  //plt.addDataArrays(N, &x2[0]);

  //plt.addDataArrays(N, &y1[0]);
  //plt.addDataArrays(N, &y2[0]);
  //plt.plot();


  // normalize observables and write into wavefiles:
  RAPT::rsArrayTools::normalize(&m1[0], N, 1.0, true);
  writeToMonoWaveFile("DoublePendulumM1.wav", &m1[0], N, fs, 16);
  RAPT::rsArrayTools::normalize(&m2[0], N, 1.0, true);
  writeToMonoWaveFile("DoublePendulumM2.wav", &m2[0], N, fs, 16);
  RAPT::rsArrayTools::normalize(&x1[0], N, 1.0, true);
  writeToMonoWaveFile("DoublePendulumX1.wav", &x1[0], N, fs, 16);
  RAPT::rsArrayTools::normalize(&x2[0], N, 1.0, true);
  writeToMonoWaveFile("DoublePendulumX2.wav", &x2[0], N, fs, 16);

  // Observations:
  // Notation: m1, m2: masses of the bobs; l1, l2: lengths of the arms
  // We start with the 1st arm horizontal and the 2nd arm upward
  // when  l2 << l1 and m2 << m1, like l1=1, l2=0.01, m1=1, m2=0.01, x2 looks like a distorted
  // sinewave and x2 almost the same.
}
// I think, in order to ensure reproducible results when using chaotic systems, it is necessarry to
// ensure that the rounding behavior is always the same (for example on different platforms, 
// compilers, etc.). There's a compiler switch for the floating point model which can be set to 
// "precise", but i'm not sure, if that's enough. The standard library functions may have different
// implementations. In any case, I think there shouldn't be any algebraically equivalent different
// formulas between versions of a software that uses such chaotic systems.
//
// Idea: 
// -start with a simple linear system with well understood and musically useful behavior, like the 
//  (driven) damped harmonic oscillator: 
//  m*a = -k*x - r*v + F 
//  where m: mass, a: acceleration, k: spring hardness, x: excursion, r: damping/resistence, 
//        v: velocity, F: driving force
// -add nonlinear terms one by one and study their influence on the system, for example a 
//  progressive spring hardening term: -k3*x^3 - this should distort the waveform and increase the 
//  frequency at higher amplitudes
// -generally, instead of using just the linear contributions (m*a, k*x, r*v), we could use 
//  polynomials in x, v, a
// -maybe we could also consider more derivatives like the increase/decrease of acceleration 
//  ("jerk"?)


void heatEquationNyquistBug()
{
  // We try to figure out what causes the parasitic osciallation at the Nyquist freq by giving
  // it an initial distribution that alternates - could it be that the Nyquist freq dies out
  // more slowly than other high frequencies? is this filter maybe not really a lowpass?
  // ...if this turns out to be the case, the easy fix is to filter it out as post-process - but
  // first, let's try to get rid of it directly in the algo

  rsHeatEquation1D<double> hteq;
  hteq.setMaxCycleLength(8);

  double initDist[6] = {1,-1,1,-1,1,-1};  // initial distribution
  hteq.setHeatDistribution(initDist, 6);

  hteq.setDiffusionCoefficient(0.1);
  // 1: nothig ever happens, 0.5: oscillation completely dies out in first iteration
  // 0.45: it dies by a factor of 10 in each iteration, 0.1: dies out with reasonable speed


  static const int N = 100;
  double x;
  for(int n = 0; n < N; n++)
    x = hteq.getSample2();



  int dummy = 0;
}

void heatEquation1D()
{
  //heatEquationNyquistBug();


  int fs           = 44100;
  int numCycles    = 128;    // plot: ~50, audio: ~1000
  int cycleLength  = 1024;   // plot: ~40, audio: ~200
  int N            = numCycles * cycleLength;    // number of samples
  double amplitude = 0.5;


  rsHeatEquation1D<double> hteq;
  hteq.setMaxCycleLength(2048);
  //hteq.setDiffusionCoefficient(0.90); // it seems to decay more slowly, the higher the coeff - ?
  hteq.setDiffusionCoefficient(0.5);
  hteq.setRandomHeatDistribution(5, cycleLength); // nice: 5,9
  //hteq.setTwoValueDistribution(0.25, cycleLength); 
  hteq.normalizeHeatDistribution();

  std::vector<double> y(N);
  for(int n = 0; n < N; n++)
    y[n] = amplitude * hteq.getSample();

  //rsPlotVector(y);

  RAPT::rsArrayTools::normalize(&y[0], N, 1.0, true);
  rosic::writeToMonoWaveFile("HeatEquation1D.wav", &y[0], N, fs);

  // todo: plot it as a 3D plot - each cycle is shown in its own plane in the (inward) 
  // z-direction

  // it's buzzy and there's a parasitic oscillation at the Nyquist freq.
  // -buzz is probably because of end-handling (try cyclic end-handling to get rid of the buzz)
  // -or maybe the buzz is due to the slight discontinuity when going from one cycle to the next
  //  (because we update the rod only once per cycle - not per sample - which wouldn't make sense
  //  in this setup)
  // -i think, the parasitic oscillation was due to choosing 1.0 as diffusion coeff - maybe it must
  //  be strictly less than 1
  //  hmmm...with 0.95, there's still s little bit of tha oscillation in the transient
  //  ...i think, with 1, the rod element would be set to exactly the average of its neighbours 
  //  and with > 1, it would overshoot the average? -> try with very short strings and with
  //  and initially alternating configuration - maybe: 1,-1,1,-1,1,-1
  // -i think, 0.5 is the maximum to avoid the nyquist oscillations
  // -different seeds give wildly different sounds - maybe try a random phase-spectrum with defined
  //  magnitude sprectrum
  // -with cyclic end-handling, it converges to some non-zero DC (fixed?)
  // -it generally sounds karplus-strong-ish - which is no surprise - the update of the rod acts as
  //  a sort of lowpass 

  // todo: produce sets of samples:
  // -to go one octave up, we should use the same (random) distribution as in the lower octave but
  //  downsampled, such that the samples are consistent
  
  //GNUPlotter plt;
}

void waveEquation1D()
{
  int numGridPoints = 65;    // 2^k + 1 are "nice" because factors are exactly representable
  int numTimeSteps  = 200;
  int width         = 10;    // width of initial impulse/displacement

  double timeStep   = 1.0 / (numGridPoints-1);  
  // optimum value - makes numerical solution exact, when waveSpeed is set to unity, if it's less,
  // we see numerical dispersion, if it's higher, the scheme becomes unstable

  // create and set up wave-equation solver object:
  rsWaveEquation1D<double> wvEq;
  wvEq.setNumGridPoints(numGridPoints);  // grid spacing is 1 / (numGridPoints-1)

  wvEq.setWaveSpeed(1.0);

  // create arrays for string state and set up initial conditions for the solver:
  int Ng = numGridPoints;  // for convenience
  std::vector<double> u(Ng), v(Ng);
  RAPT::rsArrayTools::fillWithZeros(&u[0], Ng);
  //u[Ng/2] = 1.0;
  rsWindowFunction::hanning(&u[Ng/2-width/2], width);
  RAPT::rsArrayTools::fillWithZeros(&v[0], Ng);
  wvEq.setInitialConditions(&u[0], &v[0], Ng, timeStep);

  // go through a couple of rounds of updates:
  double** plotMatrix;
  rsMatrixTools::allocateMatrix(plotMatrix, numTimeSteps, numGridPoints);
  for(int n = 0; n < numTimeSteps; n++) {
    wvEq.getState(&u[0], Ng);
    wvEq.getState(plotMatrix[n], Ng);
    //rsPlotVector(u);
    wvEq.updateState(timeStep);
  }

  // plot results:
  GNUPlotter plt;
  std::vector<double> t(numTimeSteps), x(numGridPoints);
  RAPT::rsArrayTools::fillWithIndex(&t[0], numTimeSteps);
  RAPT::rsArrayTools::fillWithIndex(&x[0], numGridPoints);
  plt.addDataMatrix(numTimeSteps, numGridPoints, &t[0],  &x[0],plotMatrix);
  // make convenience-function addDataMatrix that doesn't require x,y axes to be passed

  plt.addGraph("i 0 nonuniform matrix w image notitle");   
  plt.addCommand("set palette gray");                   // maximum is white
  plt.plot();
  // we need a bipolar palette for this

  //plt.addCommand("set hidden3d");
  //plt.plot3D();
  // todo: set hidden 3D and/or plot as heatmap

  rsMatrixTools::deallocateMatrix(plotMatrix, numTimeSteps, numGridPoints);

  // Observations:
  // -using a hanning window of width 10, it looks good - the impulses move to the boundary and get 
  //  reflected (with sign inversion), 
  // -using a single impulse-spike, the point where the spike was set oscillates at the Nyquist 
  //  freq (it alternates betwen +1 and -1) and a wave at the spatial Nyquist freq spreads out
  // -width = 4 gives a pixelated appearance
  // -waveSpeed = 0.5 shows numerical dispersion
}





// move all this plotting stuff into rs_testing module
/*
void plotMatrix(rsMatrix<double>& A, bool asHeatMap)  // use const
{
  GNUPlotter plt;
  //plt.addDataMatrixFlat( A.getNumRows(), A.getNumColumns(), A.getDataPointerConst());
  plt.addDataMatrixFlat( A.getNumRows(), A.getNumColumns(), A.getRowPointer(0));
  if(asHeatMap) {
    plt.addCommand("set size square");  // make optional
    plt.addGraph("i 0 nonuniform matrix w image notitle");   
    plt.addCommand("set palette gray");
    //plt.setRange(0, A.getNumRows()-1, 0, A.getNumColumns()-1, -1.0, +1.0); // doesn't work
    plt.plot();
  }
  else
    plt.plot3D();
}
*/
// move to rs_testing, maybe have an option to plot it as image/heatmap
// factor out a function that takes a plotter reference as argument, so we can do some setup calls
// before plotting - such as setting plotting ranges

void normalizeMatrix(rsMatrix<double>& A, bool removeMean)
{
  double *a = A.getDataPointer();
  int N = A.getSize();
  RAPT::rsArrayTools::normalize(a, N, 1.0, removeMean); // allow values other than 1
}

void plotMatricesAnimated(std::vector<rsMatrix<double>>& frames)
{
  if(frames.size() == 0)
    return;
  int numRows = frames[0].getNumRows();
  int numCols = frames[0].getNumColumns();
  GNUPlotter plt;;

  // write matrices into the datafile:
  bool normalize = true; // make parameter
  plt.setDataPrecision(2);
  for(size_t i = 0; i < frames.size(); i++) {
    rsAssert(frames[i].getNumRows() == numRows && frames[i].getNumColumns() == numCols, 
      "all matrices must have the same dimensions");
    if(!normalize)
      plt.addDataMatrixFlat(
        frames[i].getNumRows(), frames[i].getNumColumns(), frames[i].getRowPointer(0));
    else {
      rsMatrix<double> A = frames[i];
      normalizeMatrix(A, true);
      plt.addDataMatrixFlat(A.getNumRows(), A.getNumColumns(), A.getRowPointer(0));
    }
  }

  // set up gnuplot for creating animated gif:
  std::string datafile = plt.getDataPath();

  plt.addCommand("set palette gray");

  //plt.addCommand("set palette defined(-1 'green', 0 'black', 1 'red')");
  // we need a better palette - maybe something with nonlinear saturation  - the range around zero
  // needs to expanded, the ends may be compressed
  // or: optionally normalize the matrices for each frame - that would use the full contrast range
  // at the expense of giving false values

  // see here for suitable color maps:
  // https://www.kennethmoreland.com/color-maps/

  // try to get rid of the coordinate axes, color-bar, etc - just show the pure content

  plt.addCommand("set cbrange [-1:1]"); // should be set by caller - make parameters

  plt.addCommand("set terminal gif animate delay 4 optimize");
  //plt.addCommand("set terminal gif animate delay 20 optimize");
  //plt.addCommand("set terminal gif animate delay 5");
  plt.addCommand("set output 'gnuplotOutput.gif'"); 
  plt.addCommand("stats '" + datafile + "' nooutput");

  // let gnuplot loop over the frames:
  plt.addCommand("do for [i=1:int(STATS_blocks-1)] {"); 
  plt.addCommand("  plot '" + datafile + "' index (i-1) nonuniform matrix w image notitle");
  plt.addCommand("}");
  plt.invokeGNUPlot();
}
// -datafiles tend to get large - maybe we can reduce the precision - 5 significant digits should be 
//  enough for line plots, for color-coded images, 3 is actually already enough
// -with 200 frames, gnuplot takes a *really* long time to produce the gif - and yes, it's really 
//  gnuplot, not our code here
// -can we just write the images into .ppm files and use some other tool to create the animated gif 
//  from them? ...maybe even a commandline tool, so we can automate it?
//  maybe this: http://www.imagemagick.org/script/command-line-processing.php
// -try how long it takes without optimization turned on - doesn't seem to help very much - it 
//  still takes ages - it even seems as if the time increases superlinearly with the number of 
//  frames - WTF?
// -can we also produce mp4 files?

void rectangularMembrane()
{
  int numGridPoints = 65;    // using powers of two for timeStep also an (inverse)-power-of-2 / (numGridPoints-1)
  int numTimeSteps  = 200;   // with more than 200, it takes ridiculously long :-(
  int width         = 6;     // width of initial impulse/excursion
  int xPos          = 15;    // x-coordinate of initial displacement
  int yPos          = 25;    // y-coordinate of initial displacement

  double timeStep   = 1.0 / (numGridPoints-1);  
  timeStep /= sqrt(2.0);  // C = 1/sqrt(2) is the stability limit
  //timeStep /= 2;

  // set up PDE solver:
  rsRectangularMembrane<double> membrane;
  membrane.setGridDimensions(numGridPoints, numGridPoints);
  membrane.setWaveSpeed(1.0);
  membrane.setTimeStep(timeStep);


  int Ng = numGridPoints;
  rsMatrix<double> u(Ng, Ng), v(Ng, Ng);
  u.setAllValues(0);
  v.setAllValues(0);
  for(int i = xPos-width/2+1; i < xPos+width/2; i++) {    // why "+1" for start-index?
    for(int j = yPos-width/2+1; j < yPos+width/2; j++) {
      double dx = xPos - i;
      double dy = yPos - j;
      double r  = sqrt(dx*dx + dy*dy);
      u(i, j)   = 0.5 * (1 + cos(2*PI*r/width));  
      // factor out into exciteWithRaisedCosine(strength, width, x, y) 
      // or addRaisedCosine, addBellInput, addBellExcitation...maybe we should initialize with all
      // zeros and use a sort of addExcitation function - this can be used during realtime 
      // operation as well...addBellExcitation, addTriangularExcitation, addSpikeExcitation, etc
    }
  }
  membrane.setInitialConditions(u, v);

  std::vector<rsMatrix<double>> frames;
  for(int n = 0; n < numTimeSteps; n++) {
    membrane.getState(u); // maybe getState should return a matrix - but no - that would enforce allocs
    //plotMatrix(u, true);
    frames.push_back(u);
    membrane.updateState();
  }
  plotMatricesAnimated(frames);  // does nothing yet - is under construction

  // -maybe try the simplified scheme for special case lambda = 1/sqrt(2) - Eq. 11.12 and then
  //  compare with general case for a setting that would allow for simplified computations

  // -it can actually generate nice patterns when we place the initial bump at the center 
  //  -> experiment with other symmetric initial conditions - maybe make a .js version with p5.js
  //  no need to wait ages for the gif to be rendered
  //  should we use an even or odd number of datapoints for that? 
  //  -> odd (tried with 16/8/8 vs 17/8/8) - so 65/32/32 is nice
  // 

  // when rendering the gif, gnuplot constantly has disk i/o and its memory usage may be less than 
  // the size of the datafile - is this a hint that it doesn't actually read in larger datafiles at
  // once and keeps the data in ram? maybe we can set a cache-size somewhere

  // maybe we can use python to produce the video - include the PDE solver into the rsPy python 
  // module, call it from there and do visualization in python
}


void rectangularRoom()
{
  // Simulates propagation of waves in a rectangular room....

  // grid resolutions along the 3 coordinates:
  int Nx = 11;
  int Ny = 11;
  int Nz = 11;
  float dt = 0.0025f;  // time-step between two samples

  // room lengths in the 3 coordiniates (length, width, height)
  float Lx = 1.f;
  float Ly = 1.0f;
  float Lz = 1.0f;
  // NSS (pg.292) says that for isotropic problems, it's best to choose Lx,Nx etc. such that 
  // hx=hy for 2D, so in 3D isotropic, we should probably use hx=hy=hz



  rsRectangularRoom<float> room;
  room.setGridDimensions(Nx, Ny, Nz);
  room.setRoomDimensions(Lx, Ly, Lz); // maybe rename to setShape, setSizes
  room.setTimeStep(dt);

  // maybe for a user of a room-reverb, it's better to parametrize it via size, xy-ratio, xz-ratio
  // because size is a more intuitive parameter

  // todo: 
  // -initialize the room with a pressure impulse somewhere...do we need to init the u_t and 
  //  u_tt arrays, too or can they be zero
  // -run an update loop and record the pressure over time at various points of interest and plot
  //  the time series ..i hope to see something that resembles a stylized impulse response of a 
  //  room
  // -to make it more realistic, add damping...maybe frequency dependent damping... how?

  int Nt = 5000; // number of time-steps

  // indices where to put the initial impulse and to read out the signal:
  int ix = 2;
  int iy = 3;
  int iz = 4;



  std::vector<float> E_kin(Nt), E_pot(Nt), E_sec(Nt);
  std::vector<float> h(Nt), h1(Nt), h2(Nt);    // impulse response at location of excitation


  room.reset();
  room.injectPressureAt(ix, iy, iz, 1.f);
  for(int n = 0; n < Nt; n++)
  {
    E_kin[n] = room.getKineticEnergy();
    E_pot[n] = room.getPotentialEnergy() * 600;

    E_sec[n] = room.getSecondDerivativeEnergy();


    h[n]     = room.getPressureAt(ix, iy, iz);
    h1[n]    = room.getPressureDerivativeAt(ix, iy, iz) * dt;
    h2[n]    = room.getPressureSecondDerivativeAt(ix, iy, iz) * (dt*dt);


    room.updateState();
    //const rsMultiArray<float>& u = room.getState();

    // somehow visuallize the state...maybe we could plot Nz surfaces
    // maybe plot the total energy in the room (sum-of-squares of u plus sum-of-squares
    // of u_t = potential + kinetic)?
    // plot potential and kinetic energies and total energy

    int dummy = 0;
  }

  // todo: write the impulse response to a wave file

  rsPlotVectors(h, h1, h2);
  //rsPlotVectors(E_pot, E_kin, E_pot+E_kin, E_sec);
  rsPlotVectors(E_pot, E_kin, E_pot+E_kin);
  // E_pot and E_kin seem to be on a vastly different scale - figure out the scale factors from 
  // physical considerations - i think, we need to take into account the spatial and temporal
  // sampling intervals ...maybe E_kin should be multiplied by the temporal interval?
  // E_kin wiggles around 150, E_pot around 0.5 - try to figure out, how these averages behave
  // as functions of spatial and temporal sampling interval

  // Observations
  // -increasing the timeStep makes the wiggles in the energies faster (as expected) and also 
  //  increases the mean of the kinectic energy - it also seems to decrease the mean of the
  //  potential energy ..well, at least the increase from 0.01 to 0.02 increased the mean
  //  between 0.005 and 0.01, there seems to be no such increase
  //  -dt = 0.0025 gives good temporral resolution when Nx=Ny=Nz=11 and Lx=Ly=Lz=1
  //  -scaling E_pot by 600 makes both averages equal in this case - figure out where that factor
  //   comes from...maybe that factor is not yet exact - the sum of the energies still wiggles
  //  -the stability limit (for 11,11,11; 1,1,1) seems to be somewhere between dt = 0.05 and 0.06
  //   ->try to figure out Courant numbers
  // -with dt=0.0025, Nx=Ny=Nz=11, Lx=Ly=Lz=1
  //  -increasing Nx to 21 scales E_pot by factor 2 (likewise for Ny,Nz)
  //  -increasing Lx to 4 reduces E_pot's mean from 150 to 100
  //  -increasing Nx to 21 and Lx to 2 seems to reduce the amplitude of the wiggle
  //  -increasing Nx,Ny,Nz all to 21 let's the computation take a really long time
  // -the energy in the 2nd derivative is on the order of 100000 ...maybe it scales with 1/h^2 where
  //  the 1st-derivative energy scales only with 1/h

  // from Möser, page 31 - the energy density of a sound field:
  // E = (1/2) * ( p^2/(rho_0 * c^2) + rho*v^2)
  // p is the pressure, rho_0 the air-density(?), v the speed
  // rho the sound-density(

  // the 3D wave equation for the sound pressure p = p(x,y,z,t) is:
  // p_tt = c^2 * L(p) = c^2 * (p_xx + p_yy + p_zz) where L(p) is defined as the Laplacian of p
  // and c is the speed of sound
  // the Laplacian measures, by how much a value at a particular position differs from the average
  // of its neighborhood - so it makes intuively sense that p changes this way...or well, maybe one
  // could also expce p_t and not p_tt to be proportional to the Laplacian? ...but that's not how 
  // it works

  // -check out, how numerical energy is computed in Bilbao's book - but he uses schemes instead of
  //  physical variables

  // https://en.wikipedia.org/wiki/Sound_power

  // i guess, we get a lot of cache-misses when the sizes grow larger...we should use 
  // Nx <= Ny <= Nz for best cache locality...maybe with more complex data-layout using a 
  // Hilbert curve, cache locality could be improved - but that would be really complicated!
 


  int dummy = 0;
}




// maybe to really challenge the blep/blamp class, try to hardsync a sinewave and try to anti-alias
// more higher order derivatives with bladratics, blubics, blartics, etc.

void particleForceDistanceLaw()
{
  // Plots the force vs the distance of the rsPartcielSystem class for various choices of the
  // parameters.

  rsParticleSystem<float> ps(1);
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
  RAPT::rsArrayTools::fillWithRangeLinear(d, N, dMin, dMax);
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

void getTwoParticleTrajectories(rsParticleSystem<float>& ps, int N, float* x1, float* y1, float* z1,
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
  rsParticleSystem<float> ps(2);

  // both particles have unit mass, charge and size:
  ps.particles[0].mass   = 1.f;
  ps.particles[0].charge = 1.f;
  ps.particles[0].size   = 1.f;
  ps.particles[1].mass   = 1.f;
  ps.particles[1].charge = 1.f;
  ps.particles[1].size   = 1.f;


  // place them at (-1,0,0) and (+1,0,0) with zero velocity initially:
  ps.initialPositions[0]  = rsVector3D<float>(-0.5f, +0.0f,  +0.0f);
  ps.initialPositions[1]  = rsVector3D<float>(+0.5f, -0.0f,  +0.0f);
  ps.initialVelocities[0] = rsVector3D<float>( 0.0f, -0.01f, -0.0f);
  ps.initialVelocities[1] = rsVector3D<float>( 0.0f, +0.01f, +0.0f);

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
  rsVector3D<float> test = cross(rsVector3D<float>(2,3,5), rsVector3D<float>(7,11,13)); 
  // bool result = test == rsVector3D<float>(-16,9,1); // ok - cross-product is correct


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
// todo:
// -create particle systems with more complex interactions, like every particle has a number of 
//  generic "features" A,B,C,D,... represented as real numbers (think: mass, charge, etc.)
// -for each such feature, there's a specific force-distance law
// -let there be different types of particles, like type "red": A = 0.5, B = -2, C = 1.5, etc. - 
//  i.e. each type has a particular set of numeric values assigned for the different features
//  (these particle types could be seen as "atoms" of the various chemical elements)
// -one of these features could be "size" or "radius" and the corresponding force between two
//  particles should be zero if their distance is less than the sum of the radii and 
//  increase sharply, when the distance gets smaller than that -> avoid interpenetration
// -let there be a drag force that tends to slow down motion and a temperature that gives the 
//  particles random bumps into random directions at each time step
//  ...or maybe drag and temperature will also emerge if we don't define any globally?
// -with the right set of force-laws, i hope that we may get interesting emergent behavior, like 
//  particles assembling to structures which move around, grow, replicate, eat, shit, etc. - i.e. 
//  life-like behavior
// -some simple structures that would be interesting to observe (and/or create):
//  -dipoles: exert a force of certain kind into one direction and the same force with opposite 
//   sign into the opposite direction - can we assemble something like that from the given atoms?
//   -would they tend to rotate in the right circumstances (i.e. feel a torque?)
//   -or should dipoles and torques be built into the basic particles ("atoms")
// -this can be done in 2D or 3D - in 2D, it would be easy to visualize - atoms of different kind 
//  can be given different colors
// -as a further complication, the particles/atoms could have a state and behave differently in 
//  different states (i.e. the force-laws could be state dependent...or the response to forces 
//  could be state dependent)
// -it would be interesting, if energy (in the form of temperature) could be absorbed or liberated 
//  when structures assemble or disassemble (chemical bonds form or break up)
// -the force-distance-law F(d) for each feature should be a function that goes to minus inf for
//  d -> 0, cross the x-axis at some specified "preferred distance", then show a bump of positive 
//  values (height and position adjustable) and then fall down to zero. for example, a function 
//  like F(d) = a/(b + d^n) - c/d^m, for example F(d) = 2/(1+d^3) - 0.1/d^2, see
//    https://www.desmos.com/calculator/42glxvmavv 
//  the 5 free parameters a,b,c,n,m should be adjusted in terms of asymptotic repulsion m, 
//  preferred distance, bump position, bump height, asymptotic attraction n
//  -as a variation, the asymptotic atttraction could also go down exponentially
//  -the repulsion is there to model that no two particles can occupy the same position
// -see: https://www.youtube.com/watch?v=Z_zmZ23grXE


//template<class T>
//inline double rsAbs(rsVector2D<T> v)
//{
//  return (double) sqrt(v.x*v.x + v.y*v.y); // abs(vector) is defined as length(vector)
//}
// that's a bit weird naming ("absolute value of a vector") but is needed in isCloseTo for vector
// arguments - find a better way

template<class TArg, class TTol>
inline bool isCloseTo(TArg x, TArg y, TTol tol)
{
  if(rsAbs(x - y) <= tol)
    return true;
  else
    return false;
}
// move to rapt or test utilities

inline bool isCloseTo(
  rsVector2D<std::complex<double>> v, rsVector2D<std::complex<double>> w, double tol)
{
  bool r = true;
  rsVector2D<std::complex<double>> d = w-v;
  r &= d.x.real() <= tol;
  r &= d.x.imag() <= tol;
  r &= d.y.real() <= tol;
  r &= d.y.imag() <= tol;
  return r;
}

inline bool isCloseTo(
  rsMatrix2x2<std::complex<double>> A, rsMatrix2x2<std::complex<double>> B, double tol)
{
  bool r = true;
  rsMatrix2x2<std::complex<double>> D = A-B;
  r &= D.a.real() <= tol;
  r &= D.a.imag() <= tol;
  r &= D.b.real() <= tol;
  r &= D.b.imag() <= tol;
  r &= D.c.real() <= tol;
  r &= D.c.imag() <= tol;
  r &= D.d.real() <= tol;
  r &= D.d.imag() <= tol;
  return r;
}
// move to test utilities





bool quantumSpinMeasurement()
{
  // This should be turned into a unit test...
  bool pass = true;   // move to unit tests


  typedef std::complex<double> Complex;
  typedef rsQuantumSpin<double> QS;
  typedef rsVector2D<std::complex<double>>  Vec;
  typedef rsMatrix2x2<std::complex<double>> Mat;
  typedef rsNoiseGenerator<double> PRNG;

  // some work variables:
  Vec A, B, C;            // for generic states
  std::complex<double> p; // for inner products
  double P;               // for probabilities
  double tol = 1.e-13;    // tolerance for rounding errors
  double r1, r2;          // results of 1st and 2nd measurement
  std::complex<double> one(1,0), zero(0,0), two(2,0), half(.5,0), s(1/sqrt(2.0),0);
  PRNG prng;  // randum number genertor for measurements
  prng.setRange(0, 1);  // this is super important and a source for trouble if forgotten
                        // todo: allow ony a specific kind of PRNG that only has range 0..1

  // create some qubits in pure states:
  Vec u, d, l, r, i, o; // maybe use capital letters
  QS::prepareUpState(u);
  QS::prepareDownState(d);
  QS::prepareLeftState(l);
  QS::prepareRightState(r);
  QS::prepareInState(i);
  QS::prepareOutState(o);


  // check normalization of pure states:
  pass &= (p=QS::bracket(u,u)) == 1.0;
  pass &= (p=QS::bracket(d,d)) == 1.0;
  pass &= isCloseTo(p=QS::bracket(l,l), one, tol);
  pass &= isCloseTo(p=QS::bracket(r,r), one, tol);
  pass &= isCloseTo(p=QS::bracket(i,i), one, tol);  // fails
  pass &= isCloseTo(p=QS::bracket(o,o), one, tol);  // fails

  // check orthogonality:
  pass &= (p=QS::bracket(u,d)) == 0.0;  // (1) Eq 2.3
  pass &= (p=QS::bracket(d,u)) == 0.0;
  pass &= (p=QS::bracket(r,l)) == 0.0;
  pass &= (p=QS::bracket(l,r)) == 0.0;
  pass &= (p=QS::bracket(i,o)) == 0.0;  // fails
  pass &= (p=QS::bracket(o,i)) == 0.0;  // fails

  // check, if Eq 2.8 and 2.9 are satisfied:
  pass &= isCloseTo(p = QS::bracket(o,u) * QS::bracket(u,o), half, tol);
  pass &= isCloseTo(p = QS::bracket(o,d) * QS::bracket(d,o), half, tol);
  pass &= isCloseTo(p = QS::bracket(i,u) * QS::bracket(u,i), half, tol);
  pass &= isCloseTo(p = QS::bracket(i,d) * QS::bracket(d,i), half, tol);
  pass &= isCloseTo(p = QS::bracket(o,r) * QS::bracket(r,o), half, tol);
  pass &= isCloseTo(p = QS::bracket(o,l) * QS::bracket(l,o), half, tol);
  pass &= isCloseTo(p = QS::bracket(i,r) * QS::bracket(r,i), half, tol);
  pass &= isCloseTo(p = QS::bracket(i,l) * QS::bracket(l,i), half, tol);

  // check up-spin probabilities of the various pure spin states:
  pass &= isCloseTo(P = QS::getUpProbability(u), 1.0, tol); // pure up-spin   has P(up) = 1
  pass &= isCloseTo(P = QS::getUpProbability(d), 0.0, tol); // pure down-spin has P(up) = 0
  pass &= isCloseTo(P = QS::getUpProbability(r), 0.5, tol); // all other pure spin states (left, 
  pass &= isCloseTo(P = QS::getUpProbability(l), 0.5, tol); // right, in, out) have up-spin 
  pass &= isCloseTo(P = QS::getUpProbability(i), 0.5, tol); // probability of 1/2
  pass &= isCloseTo(P = QS::getUpProbability(o), 0.5, tol);

  pass &= isCloseTo(P = QS::getRightProbability(u), 0.5, tol);
  pass &= isCloseTo(P = QS::getRightProbability(d), 0.5, tol);
  pass &= isCloseTo(P = QS::getRightProbability(r), 1.0, tol);
  pass &= isCloseTo(P = QS::getRightProbability(l), 0.0, tol);
  pass &= isCloseTo(P = QS::getRightProbability(i), 0.5, tol);
  pass &= isCloseTo(P = QS::getRightProbability(o), 0.5, tol);

  pass &= isCloseTo(P = QS::getInProbability(u), 0.5, tol);
  pass &= isCloseTo(P = QS::getInProbability(d), 0.5, tol);
  pass &= isCloseTo(P = QS::getInProbability(r), 0.5, tol);
  pass &= isCloseTo(P = QS::getInProbability(l), 0.5, tol);
  pass &= isCloseTo(P = QS::getInProbability(i), 1.0, tol);
  pass &= isCloseTo(P = QS::getInProbability(o), 0.0, tol);

  // test some stuff:
  A = Vec(one, one);              // this is an invalid state because..
  P = QS::getTotalProbability(A); // ..it has total probability 2
  QS::normalizeState(A);          // this call normalizes the total probability
  pass &= isCloseTo(P = QS::getTotalProbability(A), 1.0, tol);
  A = u;
  B = d;
  C = s*u + s*d;
  pass &= (C == r);               // (1) Eq 2.5 
  pass &= (s*u - s*d == l);
  QS::randomizeState(A, &prng);
  QS::randomizeState(B, &prng);
  pass &= isCloseTo(P = QS::getTotalProbability(A), 1.0, tol);
  pass &= isCloseTo(P = QS::getTotalProbability(B), 1.0, tol);
  pass &= QS::bracket(A,B) == conj(QS::bracket(B,A));


  std::complex<double> e1, e2;
  Vec E1, E2;

  // some test with spin operators:
  Mat pauliZ, pauliY, pauliX;
  QS::setToPauliZ(pauliZ);
  QS::setToPauliY(pauliY);
  QS::setToPauliX(pauliX);
  // (1) pg 138 says that any spin operator can be written a linear combination of the 3 Pauli
  // matrices and the identity matrix

  e1 = pauliZ.getEigenvalue1();  pass &= e1 == -1.0;
  e2 = pauliZ.getEigenvalue2();  pass &= e2 == +1.0;
  E1 = pauliZ.getEigenvector1(); pass &= isCloseTo(E1, d, tol); // "down"
  E2 = pauliZ.getEigenvector2(); pass &= isCloseTo(E2, u, tol); // "up"
  // E1 ane E2 are swapped - why?  ...not true anymore?

  e1 = pauliX.getEigenvalue1();  pass &= e1 == -1.0;
  e2 = pauliX.getEigenvalue2();  pass &= e2 == +1.0;
  E1 = pauliX.getEigenvector1(); pass &= isCloseTo(E1, l, tol); // "left" - wrong - not normalized
  E2 = pauliX.getEigenvector2(); pass &= isCloseTo(E2, r, tol); // "right"

  e1 = pauliY.getEigenvalue1();  pass &= e1 == -1.0;
  e2 = pauliY.getEigenvalue2();  pass &= e2 == +1.0;
  E1 = pauliY.getEigenvector1(); pass &= isCloseTo(E1, o, tol); // "out"
  E2 = pauliY.getEigenvector2(); pass &= isCloseTo(E2, i, tol); // "in"

  // test eigenvalue and eigenvector compuation:
  Mat op;
  op.setValues(one, two, two, one);
  e1 = op.getEigenvalue1(); pass &= e1 == -1.0;
  e2 = op.getEigenvalue2(); pass &= e2 == +3.0;

  //E1 = op.getEigenvector1(); // (1, 0)     -> wrong result
  //E2 = op.getEigenvector2(); // (1,-1) * s
  //E1 = pauliZ.getEigenvector1(); 
  //E2 = pauliZ.getEigenvector2();


  // test spin measurements via Pauli matrices:
  QS::prepareDownState(A);
  p = QS::measureObservable(A, pauliZ, &prng); pass &= p == -1.0;
  p = QS::measureObservable(A, pauliZ, &prng); pass &= p == -1.0;
  P = QS::getStateProbability(A, d);           pass &= P ==  1.0;
  QS::prepareUpState(A);
  p = QS::measureObservable(A, pauliZ, &prng); pass &= p == +1.0;
  p = QS::measureObservable(A, pauliZ, &prng); pass &= p == +1.0;
  P = QS::getStateProbability(A, u);           pass &= P ==  1.0;

  QS::prepareLeftState(A);
  p = QS::measureObservable(A, pauliX, &prng); pass &= p == -1.0;
  p = QS::measureObservable(A, pauliX, &prng); pass &= p == -1.0;
  P = QS::getStateProbability(A, l);           pass &= isCloseTo(P, 1.0, tol);
  QS::prepareRightState(A);
  p = QS::measureObservable(A, pauliX, &prng); pass &= p == +1.0;
  p = QS::measureObservable(A, pauliX, &prng); pass &= p == +1.0;
  P = QS::getStateProbability(A, r);           pass &= isCloseTo(P, 1.0, tol);

  QS::prepareOutState(A);
  p = QS::measureObservable(A, pauliY, &prng); pass &= p == -1.0;
  p = QS::measureObservable(A, pauliY, &prng); pass &= p == -1.0;
  P = QS::getStateProbability(A, o);           pass &= isCloseTo(P, 1.0, tol);
  QS::prepareInState(A);
  p = QS::measureObservable(A, pauliY, &prng); pass &= p == +1.0;
  p = QS::measureObservable(A, pauliY, &prng); pass &= p == +1.0;  // fails
  P = QS::getStateProbability(A, i);           pass &= isCloseTo(P, 1.0, tol);

  // test spin measurements via dedicated functions:
  QS::prepareDownState(A);
  p = QS::measureSpinZ(A, &prng);    pass &= p == -1.0;
  p = QS::measureSpinZ(A, &prng);    pass &= p == -1.0;
  P = QS::getStateProbability(A, d); pass &= P ==  1.0;
  QS::prepareUpState(A);
  p = QS::measureSpinZ(A, &prng);    pass &= p == +1.0;
  p = QS::measureSpinZ(A, &prng);    pass &= p == +1.0;
  P = QS::getStateProbability(A, u); pass &= P ==  1.0;

  QS::prepareLeftState(A);
  p = QS::measureSpinX(A, &prng);    pass &= p == -1.0;
  p = QS::measureSpinX(A, &prng);    pass &= p == -1.0;
  P = QS::getStateProbability(A, l); pass &= isCloseTo(P, 1.0, tol);
  QS::prepareRightState(A);
  p = QS::measureSpinX(A, &prng);    pass &= p == +1.0;
  p = QS::measureSpinX(A, &prng);    pass &= p == +1.0;
  P = QS::getStateProbability(A, r); pass &= isCloseTo(P, 1.0, tol);

  QS::prepareOutState(A);
  p = QS::measureSpinY(A, &prng);    pass &= p == -1.0;
  p = QS::measureSpinY(A, &prng);    pass &= p == -1.0;
  P = QS::getStateProbability(A, o); pass &= isCloseTo(P, 1.0, tol);
  QS::prepareInState(A);
  p = QS::measureSpinY(A, &prng);    pass &= p == +1.0;
  p = QS::measureSpinY(A, &prng);    pass &= p == +1.0;
  P = QS::getStateProbability(A, i); pass &= isCloseTo(P, 1.0, tol);


  // test with random states, if matrix based and dedicated function based measurements do the
  // same thing:
  int n;
  int N = 1000;
  rsNoiseGenerator<double> prng1, prng2;
  prng1.setRange(0.0, 1.0);
  prng1.reset();
  prng2.setRange(0.0, 1.0);
  prng2.reset();
  for(n = 0; n < N; n++)
  {
    QS::randomizeState(B, &prng);

    A = B; r1 = QS::measureObservable(A, pauliZ, &prng1);
    A = B; r2 = QS::measureSpinZ(A, &prng2);
    pass &= (r1 == r2);

    A = B; r1 = QS::measureObservable(A, pauliX, &prng1);
    A = B; r2 = QS::measureSpinX(A, &prng2);
    pass &= (r1 == r2);

    A = B; r1 = QS::measureObservable(A, pauliY, &prng1);
    A = B; r2 = QS::measureSpinY(A, &prng2);
    pass &= (r1 == r2);
  }
  // todo: make additional meausrements after the collapse - they should always give the same 
  // results

  // test, if the statistical distribution is as desired - set it into a state
  // au = sqrt(0.8), ad = sqrt(0.2) - we should see roughly 80% "up" measurements and 20% "down"
  B = Vec(std::complex<double>(sqrt(0.8), 0), std::complex<double>(sqrt(0.2), 0));
  std::vector<double> spins(N);
  for(n = 0; n < N; n++) {
    A = B;
    spins[n] = QS::measureSpinZ(A, &prng);
  }
  double mean = rsMean(spins);
  P = (mean+1)/2;  // -1..+1 -> 0..1  P = 0.81 (with N=1000) - close enough to 0.8

  // compute expectation values for the spin component measurements in state B:
  double Ez, Ex, Ey;  
  Ez = QS::getExpectedMeasurement(pauliZ, B); // should be 0.6 = (0.8*2)-1 and close to mean1
  Ex = QS::getExpectedMeasurement(pauliX, B); // 0.8
  Ey = QS::getExpectedMeasurement(pauliY, B); // 0.0
  P  = Ex*Ex + Ey*Ey + Ez*Ez;  // should be one according to (1) Eq 3.27
  pass &= isCloseTo(P,  1.0, tol);
  pass &= isCloseTo(Ez, 0.6, tol);
  // Ex = 0.8 and Ey = 0 - maybe try to compute these manually


  // todo:
  // test uncertainty computation:
  // -if the system is in an eigenstate of an observable (say PauliZ), this should be zero for
  //  that observable, i.e. uncertainty(pauliZ, A) = 0, iff A is eigenstate of pauliZ
  // -if the system is in eigenstate pauliX, uncertainty of pauliZ should be maximal
  // -is the uncertainty product of pauliX,pauliZ independent of the state? i think so

  // verify (1) Eq 7.11 for our 3 sets of basis vectors:
  Mat M, M2, L;
  Mat I = Mat::identity();
  M = QS::projector(u) + QS::projector(d); pass &= M == I;
  M = QS::projector(r) + QS::projector(l); pass &= isCloseTo(M, I, tol);
  M = QS::projector(i) + QS::projector(o); pass &= isCloseTo(M, I, tol);

  // test for a random (normalized) state A, if A is an eigenvector of its projector |A><A| with
  // eigenvalue 1:
  QS::randomizeState(A, &prng);
  //QS::normalizeState(A); // not necessarry
  M  = QS::projector(A);
  e1 = M.getEigenvalue1();   // 0
  e2 = M.getEigenvalue2();   // 1
  E1 = M.getEigenvector1();  // 
  E2 = M.getEigenvector2();  // is this = k * A for some (complex) constant k
  Complex c1, c2;         // move decalaration to top
  c1 = E2.x / A.x;
  c2 = E2.y / A.y;
  pass &= isCloseTo(c1, c2, tol);   // if c1 == c2, E2 is a multiple of A, so A is eigenvetor of M
  c1 = QS::bracket(E1, E2);         
  pass &= isCloseTo(c1, zero, tol); // E1 should be orthogonal to E2
  M2 = M*M;
  pass &= isCloseTo(M, M2, tol);    // M should be equal to its square
  c1 = QS::sandwich(A, pauliZ, A); // todo: use some more general observable instead of pauliZ
  c2 = (M*pauliZ).getTrace();
  pass &= isCloseTo(c1, c2, tol);  // should be equal (for any observable) by (1) Eq 7.12
  QS::randomizeState(B, &prng);
  L = QS::projector(B);            // since projectors are Hermitian, L is an observable (right?)
  c1 = QS::sandwich(A, L, A);
  c2 = (M*L).getTrace();
  pass &= isCloseTo(c1, c2, tol);  // should be equal (for any observable) by (1) Eq 7.12



  // density matrix stuff:



  // ...
  // todo: test Eq. 7.13, page 198: 
  // -create two random states psi, phi and their projectors projPsi, projPhi
  // -assign probabilities probPsi, probPhi and create a denstity matrix 
  //  rho = probPsi*projPsi + probPhi*projPhi
  // -create a random observable L (maybe by using the projector of yet another random state)
  // -compute both sides of 7.13
  Mat PA, PB;  // projectors of A and B
  QS::randomizeState(A, &prng);   // state psi (created by Alice with probability 0.7)
  QS::randomizeState(B, &prng);   // state phi (created by Alice with probability 0.3)
  QS::randomizeState(C, &prng);   // used as observable by Bob
  PA = QS::projector(A);
  PB = QS::projector(B);
  L  = pauliZ;
  Complex EA, EB, EL;
  Mat rho = QS::densityMatrix(std::vector<double>({ 0.75, 0.25 }), std::vector<Vec>({ A, B }));
  EA = QS::sandwich(A,L,A); pass &= EA == (PA*L).getTrace();  // 7.12
  EB = QS::sandwich(B,L,B); pass &= EB == (PB*L).getTrace();  // 7.12
  EL = (rho * L).getTrace();   pass &= isCloseTo(EL, 0.75*EA + 0.25*EB, tol);



  // ...oh - wait - there is no other way to evaluate the lhs of eq 7.13 - hmm - scrap that test

  // test 7.10, pg 192: (A (x) B) * (a (x) b) = A*a (x) B*b where (x) is the Kronecker product and
  // * is the matrix product - does this relation hold for matrices of any size or just for NxN 
  // matrices A,B and Nx1 column vectors a,b?

  // what about the eq in the middel of pg 195: Tr(L) = sum_i sandwich(i, L, i)
  // or pg 197: Tr(projector(phi) * L) = sandwich(phi, L, phi)



// -what about eq 71. pg 185 - is this an implicit equation for any observable M?

  rsAssert(pass);
  return pass;
}

//bool quantumSpinMeasurement()

/** Checks, if the vector v is an eigenvector of the matrix A by verifying if A*v equals some 
multiple of v within a given tolerance. */
template<class TElem, class TTol>
bool isEigenVector(const rsMatrix<TElem>& A, const rsMatrix<TElem>& v, TTol tol)
{
  rsAssert(v.isColumnVector());

  // figure out the element in v with the largest absolute value - if it is zero, then the vector
  // v must be the zero vector which is an eigenvector of any matrix:
  int N = v.getNumRows(); // dimensionality of v
  int n = 0;
  TTol max = 0; 
  for(int i = 0; i < N; i++) {
    if(rsAbs(v.at(i, 0)) > max) {
      max = rsAbs(v.at(i,0));
      n = i; }}
  if(v.at(n,0) == TElem(0))
    return true;

  // the largest (absolute) element in v (which is at at index n) in nonzero - compute the matrix
  // vector product w = A*v and the ratio of w[n]/v[n] - if v is an eigenvector, this ratio should
  // be the corresponding eigenvalue (we call it lambda) and all elements w[i] should be in the 
  // same ratio, i.e. w[i] should be lambda*v[i] for all i:
  rsMatrix<TElem> w = A*v;
  TElem lambda = w.at(n,0) / v.at(n,0);
  for(n = 0; n < N; n++) {
    if(rsAbs(lambda*v.at(n,0) - w.at(n,0)) > tol)
      return false; }
  return true;
}
// maybe rename to isRightEigenVector, make a similar function isLeftEigenVector
// move into matrix class
// maybe, it could return the computed eigenvalue?

// todo: implement equations for entanglement - this becomes a new function
bool quantumSpinEntanglement()
{
  bool pass = true;  // move to unit tests

  typedef std::complex<double> Complex;
  typedef rsMatrix<Complex> Mat;
  typedef std::vector<Complex> Vec;
  typedef rsQuantumState<double> QS;

  double tol = 1.e-13;


  //typedef rsQuantumSpin<double> QS;
  //typedef rsVector2D<std::complex<double>>  Vec;
  //typedef rsMatrix2x2<std::complex<double>> Mat;
  //Complex S[4];  // our state (or wave function), consisting of amplitudes for the 4 basis vectors
  // Note that the fact that we need 4 complex numbers (amplitudes) to represent the state is 
  // because 2^2 = 4 and not because 2*2 = 4. A state made from 3 spins would need 2^3 = 8 
  // amplitudes and we can't represent that as an array of 3 spins
  // see 136...
  // the measurement of one of the two spins can be simulated by figuring out the probabilities of 
  // the 4 eigenstates of the 4x4 matrix representing the measurement of the respective spin (it's
  // constructed by taking the kronecker product of the 2x2 1-spin operator and the 1x1 identity 
  // matrix), selecting one at random according to the probabilities and setting the whole 4D 
  // system into the corrsponding eigenstate - which happens to be a product state which implies
  // that it will also fix what will be mesured for the corresponding observable in the other spin

  //Vec S[2];  // our state, consisting of 2 spins

  // create base states |uu>, |ud>, |du>, |dd> of the cobined system of two spins (see pg 189):
  Complex one(1,0), zero(0,0), i(0,1);
  Mat u(2, 1, Vec({ one, zero }));            // |u> - spin up
  Mat d(2, 1, Vec({ zero, one }));            // |d> - spin down
  Mat uu = Mat::getKroneckerProduct(u, u);    // |uu> - up/up
  Mat ud = Mat::getKroneckerProduct(u, d);    // |ud> - up/down
  Mat du = Mat::getKroneckerProduct(d, u);    // |du> - down/up
  Mat dd = Mat::getKroneckerProduct(d, d);    // |dd> - down/down

  // create observables sz, tz (simga-z, tau-z) and the product observables (pg 190)
  Mat pauliZ(2, 2, Vec({ one,  zero, zero, -one  }));
  Mat pauliX(2, 2, Vec({ zero, one,  one,   zero }));
  Mat pauliY(2, 2, Vec({ zero, -i,   i,     zero }));
  Mat id2x2( 2, 2, Vec({ one,  zero, zero,  one  }));     // 2x2 identity matrix
  Mat sztx = Mat::getKroneckerProduct(pauliZ, pauliX);    // sigma_z  (x) tau_x
  Mat sxtz = Mat::getKroneckerProduct(pauliX, pauliZ);    // sigma_x  (x) tau_z
  Mat tzsz = Mat::getKroneckerProduct(pauliZ, pauliZ);    // tau_z (x) sigma_z
  Mat txsx = Mat::getKroneckerProduct(pauliX, pauliX);    // tau_x (x) sigma_x
  Mat tysy = Mat::getKroneckerProduct(pauliY, pauliY);    // tau_y (x) sigma_y
  Mat st   = txsx + tysy + tzsz;                          // sigma * tau, pg 180
  Mat szI  = Mat::getKroneckerProduct(pauliZ,  id2x2);    // sigma_z (x) identity, Eq 7.4, pg 187

  Mat sxI  = Mat::getKroneckerProduct(pauliX,  id2x2);
  Mat syI  = Mat::getKroneckerProduct(pauliY,  id2x2);

  Mat Itx  = Mat::getKroneckerProduct(id2x2,  pauliX);    // identity (x) tau_x
    // see page 170 bottom "if we were being pedantic,..." - yes, we are!

  // check some of the relations on page 350 (relations that hold for the Alice- and 
  // Bob-observables):
  pass &= szI * uu ==  uu;
  pass &= szI * ud ==  ud;
  pass &= szI * du == -du;
  pass &= szI * dd == -dd;
  pass &= Itx * uu ==  ud;
  pass &= Itx * ud ==  uu;
  pass &= Itx * du ==  dd;
  pass &= Itx * dd ==  du;
  // ...maybe check more - maybe calso relations on 351

  // check effect of our mixed observables on some states:
  pass &= sztx * ud == uu;

  // Check, if base states are eigenvectors of observables sigma_z and tau_x:
  pass &= isEigenVector(szI, uu, tol); // +1 (eigenvalue)
  pass &= isEigenVector(szI, ud, tol); // +1
  pass &= isEigenVector(szI, du, tol); // -1
  pass &= isEigenVector(szI, dd, tol); // -1
  // maybe, we should also check the eigenvalues?

  // oh - for these, it fails - aren't they supposed to be eigenvectors of tau_x, too?
  // or wait - no - they are already checked above - the equation on page 350
  //pass &= isEigenVector(Itx, uu, tol);
  //pass &= isEigenVector(Itx, ud, tol);
  //pass &= isEigenVector(Itx, du, tol);
  //pass &= isEigenVector(Itx, dd, tol);


  // what are the eigenvectors of sztx, sxtz? are these the maximally entangled states?

  // create singlet state and the 3 triplet states (pg 166):
  Complex s = Complex(1. / sqrt(2), 0.0);
  Mat sing = s * (ud - du);   // singlet state
  Mat trp1 = s * (ud + du);   // triplet state 1
  Mat trp2 = s * (uu + dd);   // triplet state 2
  Mat trp3 = s * (uu - dd);   // triplet state 3

  // check some relations for these states
  pass &= tzsz * sing == -sing;  // pg 177, |sing>| is eigenvector with eigenvalue -1
  pass &= txsx * sing == -sing;  // pg 178

  pass &= isEigenVector(st, sing, tol);  // pg 181
  pass &= isEigenVector(st, trp1, tol);  
  pass &= isEigenVector(st, trp2, tol);  
  pass &= isEigenVector(st, trp3, tol);  
  //

  // check 6.11 - the sum of expectation values of sx, sy, sz should be one for product states
  //Mat v = uu; 
  Complex Ex, Ey, Ez;
  Ex = QS::sandwich(uu, sxI, uu);  // 0
  Ey = QS::sandwich(uu, syI, uu);  // 0
  Ez = QS::sandwich(uu, szI, uu);  // 1
  pass &= Ex == zero && Ey == zero && Ez == one;
  Ex = QS::sandwich(ud, sxI, ud);  // 0
  Ey = QS::sandwich(ud, syI, ud);  // 0
  Ez = QS::sandwich(ud, szI, ud);  // 1
  pass &= Ex == zero && Ey == zero && Ez == one;
  Ex = QS::sandwich(du, sxI, du);  // 0
  Ey = QS::sandwich(du, syI, du);  // 0
  Ez = QS::sandwich(du, szI, du);  // -1
  pass &= Ex == zero && Ey == zero && Ez == -one;
  Ex = QS::sandwich(dd, sxI, dd);  // 0
  Ey = QS::sandwich(dd, syI, dd);  // 0
  Ez = QS::sandwich(dd, szI, dd);  // -1
  pass &= Ex == zero && Ey == zero && Ez == -one;

  // For a more general product state, the relation Ex + Ey + Ez = 1 should still hold (Eq 6.11). 
  // This is not true for entangled states. For a maximally entangled state, this sum should become
  // zero, indicating that the measurement of z-spin is totally unpredictable when the system
  // is in a maximally entangeld state - on the other hand, if we do make the measurment, we can
  // predict with 100% probability, that the z-spin of the other half of the state (Bob's half, 
  // Itz, tau-z) will be the opposite from our "Alice" measurement

  // pg 173/174:
  Ex = QS::sandwich(sing, sxI, sing);  // 0
  Ey = QS::sandwich(sing, syI, sing);  // 0
  Ez = QS::sandwich(sing, szI, sing);  // 0
  pass &= Ez == zero;



  // todo: check desnity matrices for entangled states


  // check 7.10


  rsAssert(pass);
  return pass;
}

// -check pg 166, singlet and triple states
// -pg 188, eq 7.7: c kronekcer product of matrices












// gates:
// this here:
// https://www.youtube.com/watch?v=ZN0lhYU1f5Q
// says: measure, hadamard, phase, T (rotate |1> by pi/4), cnot

// https://homepages.cwi.nl/~rdewolf/qcnotes.pdf
// http://mmrc.amss.cas.cn/tlb/201702/W020170224608150244118.pdf


template<class T>
class rsQuantumGate
{

public:

  typedef rsVector2D<std::complex<T>> QBit;

  /** Acts as a logical "not" on a single qubit. */
  static void pauliX(QBit& q) 
  {
    QBit t = q;  // temporary
    q.x = t.y;
    q.y = t.x;
  }

  static void hadamard(QBit& q)
  {
    T s = T(1)/sqrt(T(2));  // maybe make static member, if used often
    QBit t = q;  // temporary
    q.x = s*t.x + s*t.y;
    q.y = s*t.x - s*t.y;
  }

  /*
  static void toffoli(Vec& v1, Vec& v2, Vec& v3)
  {

  }
  */

};
// https://en.wikipedia.org/wiki/Quantum_logic_gate
// https://medium.com/@jonathan_hui/qc-programming-with-quantum-gates-8996b667d256
// https://medium.com/@jonathan_hui/qc-programming-with-quantum-gates-2-qubit-operator-871528d136db

bool quantumGates()
{  
  bool pass = true;  // move to unit tests

  typedef std::complex<double> Complex;
  typedef rsVector2D<Complex> QBit;
  typedef rsQuantumSpin<double> QS;
  typedef rsQuantumGate<double> QG;
  typedef rsMatrix2x2<Complex> Mat2;  // for comparing actions of QG to matrix multiplications

  QBit q1, q2, q3, q4;  // a couple of qbits

  // test Pauli-X:
  QS::prepareUpState(q1);    // q1 = |1>  (up)
  QS::prepareDownState(q2);  // q2 = |0>  (down)
  q3 = q1;                   // q3 = q1 = |1)
  QG::pauliX(q3);            // q3 = |0> = q2, bcs Pauli-X acts as logical "not"
  pass &= q3 == q2;
  //...or is |up> used as |0> and |down> as |1>? look up what is conventional
  // seems so: https://en.wikipedia.org/wiki/Bloch_sphere


  Mat2 pauliX; QS::setToPauliX(pauliX); // not gate


  rsQuantumComputer<double> cmp;
  cmp.applyGate(pauliX, 0);
  cmp.applyGate(pauliX, 1);
  cmp.applyGate(pauliX, 2);
  cmp.applyGate(pauliX, 3);




  rsAssert(pass);
  return pass;
}
// i don't really understand, how the interaction between teh states via the gates works - how
// do staes get entangled and how can we simulate this entanglement when applying gates? say
// q1 is entengled with q2 and we apply pauli-x ("not") to q1. something should automatically
// happen with q2, too - i guess....
// https://www.monoidal.net/papers/tutorialqpl-1.pdf
// https://www.monoidal.net/papers/tutorialqpl-2.pdf

// see here - Quirk - a toy quantum computer simulator
// https://www.youtube.com/watch?v=aloFwlBUwsQ
// https://algassert.com/quirk
// https://github.com/Strilanc/Quirk

// this is a really good paper on how to simulate a quantum computer
// https://arxiv.org/pdf/1805.00988.pdf



template<class T>
void plotQuantumSpinStateTrajectory(
  const std::vector<rsVector2D<std::complex<T>>>& Psi, T dt)
{
  typedef rsQuantumSpin<double> QS;
  int N = (int) Psi.size();
  std::vector<T> t(N), dr(N), di(N), ur(N), ui(N);
  for(int n = 0; n < N; n++) {
    t[n] = n*dt;   // time axis
    ur[n] = QS::getUpAmplitude(  Psi[n]).real();
    ui[n] = QS::getUpAmplitude(  Psi[n]).imag();
    dr[n] = QS::getDownAmplitude(Psi[n]).real();
    di[n] = QS::getDownAmplitude(Psi[n]).imag();
  }
  GNUPlotter plt;
  plt.addDataArrays(N, &t[0],  &ur[0], &ui[0], &dr[0], &di[0]);
  plt.plot();
}

bool quantumSpinEvolution()
{
  bool pass = true;   // move to unit tests - or maybe not

  // Simulates a spin (i.e. a spinning electron) in a magnetic field as explained in (1) pg 116 ff. 
  // The magnetic field is aligned with the z-axis. For a classical spinning charge, the energy of 
  // the system is proportional to the dot product of the spin and the field. In quantum mechanics,
  // the Hamiltonian H is proportional to the z-spin operator (the observable associated with the
  // Pauli matrix sigma_z). We solve the time dependent Schrödinger equation both numerically  and
  // analytically and plot the trajectories of the real and imaginary parts of the probability 
  // amplitudes of |u> and |d> ("up" and "down"). From the numeric solution, we also plot the 
  // expectation values of the spin measurements along the 3 axes. The expectation value for the 
  // z-spin is constant while the other two oscillate. We also plot the difference between analytic
  // an numeric solution representig the error of the numerical solver algorithm (which is forward 
  // Euler)
  //
  // References:
  //  (1) The Theoretical Minimum - Quantum Mechanics (Leonard Susskind, Art Friedman)

  int n, N    = 2000;  // loop index and number of samples
  double w    = 1;     // angular frequency?
  double hBar = 1;     // we use https://en.wikipedia.org/wiki/Planck_units

  // todo: let user define a B-vector (for the magnitic field) and use that for the Hamiltonian by
  // forming the scalar product sigma * B where B is the 3 vector defining the magnetic field
  // ((1) page 116)

  typedef std::complex<double> Complex;
  typedef rsQuantumSpin<double> QS;
  typedef rsVector2D<Complex>  Vec;
  typedef rsMatrix2x2<Complex> Mat;

  // set up the random number generator to be used for measurements:
  rsNoiseGenerator<double> prng;
  prng.setSeed(420);
  prng.setRange(0.0, 1.0);

  // Create an initial spin state:
  Vec Psi0;
  QS::randomizeState(Psi0, &prng);
  //QS::prepareUpState(Psi0);
  //QS::prepareDownState(Psi0);
  //QS::prepareRightState(Psi0);
  //QS::prepareLeftState(Psi0);
  //QS::prepareInState(Psi0);
  //QS::prepareOutState(Psi0);

  double p = 2;             // phase - tweaking this is intersting
  double k = 1 / sqrt(2.);  // to normalize the state
  double s = k * sin(p);
  double c = k * cos(p);
  //QS::prepareState(Psi0, c, s, s, -c); 

  // create Pauli matrices and Hamiltonian:
  rsVector3D<rsMatrix2x2<Complex>> sigma = QS::pauliVector();
  Mat H = Complex(hBar*w/2) * sigma.z; // (1) Eq 4.23


  // Solve the time dependent Schrödinger equation numerically using the forward Euler method:
  Complex i(0, 1);  // imaginary unit
  std::vector<Vec> stateTrajectory(N);     // records the trajectory of our spin state Psi
  std::vector<double> ex(N), ey(N), ez(N); // expectation values of spin component measurements
  double step = 0.01;                      // integration step size
  Complex cStep = Complex(step);           // ...needs to be complexified 
  Vec Psi = Psi0;                          // initialize state vector Psi(t) = Psi(0)
  for(n = 0; n < N; n++) 
  {
    // record time series:
    stateTrajectory[n] = Psi;
    ex[n] = QS::getExpectedMeasurement(sigma.x, Psi);
    ey[n] = QS::getExpectedMeasurement(sigma.y, Psi);
    ez[n] = QS::getExpectedMeasurement(sigma.z, Psi);

    // update state:
    Vec dPsi = -i * H * Psi;   // (1) Eq 4.9 or 4.10 (time dependent Schrödinger equation)
    Psi = Psi + cStep * dPsi;  // forward Euler step
    QS::normalizeState(Psi);   // avoid divergence due to error build up
  }
  plotQuantumSpinStateTrajectory(stateTrajectory, step);
  rsPlotVectors(ex, ey, ez);


  // Solve the time dependent Schrödinger equation analytically using the "Recipe for a 
  // Schrödinger Ket" in (1):
  Vec E1 = H.getEigenvector1();
  Vec E2 = H.getEigenvector2();
  Complex e1  = H.getEigenvalue1();
  Complex e2  = H.getEigenvalue2();
  Complex a10 = QS::bracket(E1, Psi0);   // initial multiplier for E1, (1) Eq 4.31
  Complex a20 = QS::bracket(E2, Psi0);   // initial multiplier for E2, (1) Eq 4.31
  std::vector<Vec> stateTrajectory2(N);
  Psi = Psi0;                            // re-initialize time dependent state vector
  Complex a1, a2;                        // time dependent multipliers
  for(n = 0; n < N; n++)
  {
    // compute state analytically:
    double t = step * n;
    a1  = a10 * exp(-i/hBar * e1 * t);  // (1) Eq 4.30
    a2  = a20 * exp(-i/hBar * e2 * t);
    Psi = a1*E1 + a2*E2;                // (1) Eq 4.29 

    // record time series:
    stateTrajectory2[n] = Psi;
  }
  plotQuantumSpinStateTrajectory(stateTrajectory2, step);

  // Compute error of numeric solution and plot it:
  std::vector<Vec> error = stateTrajectory2 - stateTrajectory;
  plotQuantumSpinStateTrajectory(error, step);



  // Verify experimentally some formulas in (1):

  // compute commutators of the observables represented by the Pauli matrices X,Y,Z with the 
  // Hamiltonian H (see (1) Eq. 4.24):
  Mat ZH = Mat::commutator(sigma.z, H);  // 0
  Mat XH = Mat::commutator(sigma.x, H);
  Mat YH = Mat::commutator(sigma.y, H);
  // ZH is the zero matrix. This implies that (d/dt) <Z> = (-i/hBar) * <[Z,H]> = 0 (Eq 4.18). This 
  // means that the change (with time) of the expectation value is zero, i.e. the expectation value 
  // of Z is constant over time. In general, the expectation value of observables that commute with
  // the Hamiltonian are conserved (see (1) pg 115)
  // compare results to (1) Eq 4.27

  // compute commutators of Pauli matrices (see (1) Eq. 4.26)
  Mat XY = Mat::commutator(sigma.x, sigma.y);  // 2*i*Z
  Mat YZ = Mat::commutator(sigma.y, sigma.z);  // 2*i*X
  Mat ZX = Mat::commutator(sigma.z, sigma.x);  // 2*i*Y
 
  // Observations:
  // -The real and imaginary parts of au and ad move sinusoidally with a frequency determined by w
  //  ->figure out the formula for the frequency (scale the time axis appropriately)
  // -when starting it |r> state, the sinusoidal amplitudes are all equal
  // -when starting in |u> state, u has amplitude 1 and d amplitude u, u.re is cosine and u.im is
  //  negative sine
  //  -ez is constant 1 and ex, ey are constant 0
  // -when starting in |d> state, d.r = cos, d.i = sin, u.r = u.i = 0
  // -when starting in |i> state, u and d pairs actually are not phase-shifted with respect to one
  //  another
  // -when starting in |o> state, u and d seem to be shfted 180° degrees 
  //  ur = cos, ui = -sin, dr = sin, di = -cos - so we see the full qudruplet of phase shifts - a 
  //  rather nice state to start with
  // -blue/black seem to be alway 90° shifted with respect to each other, same for grren/red
  // -Q: what's the phase-shift between the blue/black pair and the green/red pair? is this fixed
  //  or does it depend on the seed?
  // -with a random start state, we my get any constant value for ez (around 0.4 for seed=420) and 
  //  ex,ey vary sinusodially - but the frequency is higher than that of re/im components of u/d 
  //  (twice as high? maybe that's because an absolute value is involved which itself involves 
  //  taking squares? -> figure out)
  // -ez always being a constant
  // todo:
  // -figure out, how the phase relationships depend on the initial state
  // -maybe try to find a "nice" initial state where 
  //  u.re = cos(phi), u.im = sin(phi), d.re = -cos(phi), d.im = -sin(phi)

  // Notes:
  // -The time dependent Schrödinger equation dPsi = -i * H * Psi actually encapsulates a system of
  //  4 differential equations - the state vector Psi has two elements and each is a complex number 
  //  (todo: write the 4 equations out in full). If the state vector would be N-dimensional, it 
  //  would be a system of 2*N (first order) differential equations. If the state "vector" would be 
  //  a continuous function, it would be a partial differential equation and H would be a continuous
  //  differential operator (right?). Maybe it's also worthwhile to think about the trivial quantum
  //  system that has only a single possible value to be measured, i.e. N=1. I think, in this case,
  //  the Schrödinger equation just rotates a single complex number around the unit circle? Such a 
  //  motion can be described by a system of two (real) differential equations 
  //  (iirc: x' = -y, y' = x or something), so everything fits.
  // -It seems, we are getting some sort of rotation in 4D space

  // -According to Eq 4.17, for any observable represented by a matrix L, we have
  //  (d/dt) <L> = (i/hBar) * <[H,L]>, where <L> is the expectation value of observable L and 
  //  [H,L] is the commutator of H and L i.e. H*L-L*H (and <[H,L]> is its expectation value)
  //  i*[H,L] is Hermitatian (i.ie an observable), if H and L are both Hermitian - the i is 
  //  important, [H,L] itself is not Hermitian (says the book)
  //  ->test this experimentally

  // todo: take as observable not the spin along the z-axis but along an arbitrary axis (i think, along 
  // the z-axis, the spin is constant anyway?

  // todo: factor out the computation of a state trajectory, given an initial state, a Hamiltonian
  // a number of steps and a step-size
  // move the plotting into a different function

  // plot not only the state but also the expectation value(s) of some observable(s) as function(s)
  // of time, see if this matches hat Eq 4.17 predicts

  // todo: can we accumulate all the steps into a single big step represented by a matrix U(t) 
  // (the time-development operator, see Eq 4.1) that has accumulated all the little steps up to t?


  // todo:
  // Maybe this can be turned into an oscillatore ...like rsQuantumSpinOsc or something? it creates
  // 4 sines that have a certain pahses relationship between them - can this be useful?

  // see also:
  // https://en.wikipedia.org/wiki/Quantum_state#Pure_states

  return pass;
}


// maybe combine a 2-state system with a 3-state system (for an example 3-state system, see


void quantum3StateSystem()
{
  // we implemenent a numerical soultion to:
  // http://folk.uio.no/jaakko/FYS3110/Griffiths_3.38.pdf

  typedef complex<double> Complex;
  typedef rsMatrix<Complex> Mat;
  typedef vector<Complex> Vec;

  Mat H(3, 3, Vec({1,0,0, 0,2,0, 0,0,2}));
  Complex lambda = 2, mu = 3;

  Mat A = lambda * Mat(3, 3, Vec({0,1,0, 1,0,0, 0,0,2}));
  Mat B = mu     * Mat(3, 3, Vec({2,0,0, 0,0,1, 0,1,0}));


  int dummy = 0;

}


// this is just some quick throw-away code and can be deleted soon:
double hermitePolynomial(double x, int n)
{
  if(n == 0) return 1;
  if(n == 1) return 2*x;
  double h0 = 1;
  double h1 = 2*x;
  for(int i = 1; i < n; i++)
  {
    double tmp = 2 * (x*h1 - i*h0); // H[n+1](x) = 2*x*H[n](x) - 2*n*H[n-1](x)
    h0 = h1;
    h1 = tmp;
  }
  return h1;
}
// this is the physicist's version of Hermite polynomials, the probabilist's version would use the
// recursion: tmp = x*h1 - i*h0. ...without the factor two
double hermitePolynomial2(double x, int n)  // just for test
{
  if(n == 0) return  1;
  if(n == 1) return  2*x;
  if(n == 2) return  4*x*x - 2;
  if(n == 3) return  8*x*x*x - 12*x;
  if(n == 4) return 16*x*x*x*x - 48*x*x + 12;
  if(n == 5) return 32*x*x*x*x*x - 160*x*x*x + 120*x;
  return 0;
}
void testHermitePoly()
{
  GNUPlotter plt;
  int N = 501;
  std::vector<double> x(N), y(N);
  RAPT::rsArrayTools::fillWithRangeLinear(&x[0], N, -2.0, +2.0);
  double tmp1, tmp2;
  tmp1 = hermitePolynomial( 5.0, 2);
  tmp2 = hermitePolynomial2(5.0, 2);
  for(int i = 0; i < 5; i++) {
    for(int n = 0; n < N; n++)
      y[n] = rsPolynomial<double>::evaluateHermite(x[n], i);
    plt.addDataArrays(N, &x[0], &y[0]);
  }
  plt.plot();
}

template<class T>
void normalize(std::complex<T>* x, int N, T amplitude = 1.0)
{
  T max = 0.0;
  for(int n = 0; n < N; n++) {
    if(abs(x[n].real()) > max) max = abs(x[n].real());
    if(abs(x[n].imag()) > max) max = abs(x[n].imag());
  }
  T s = amplitude / max;
  for(int n = 0; n < N; n++)
    x[n] *= s;
}

void quantumParticle()
{
  //testHermitePoly();  // throw away soon

  int    Nx    = 128;    // number of spatial samples
  double xMin  = -1.0;   // minimum x-coordinate
  double xMax  =  1.0;   // maximum x-coordinate
  double dt    = 1.e-5;  // time step for Euler solver

  double w     = 40;     // angular frequency (?) - controls actually width of initial distribution
  int    E     = 0;      // energy level
  double scl   = 5;      // scales x in the Hermite polynomial

  int numCycles = 400;  // number of cycles to record
  int skipRatio = 200;   // number of iterations to per recorded cycle (i.e. oversampling)
                         // ...needs more

  double dx = (xMax-xMin)/(Nx-1);  // spatial sampling rate
  double r  = dt / (dx*dx);        // should be <= 1/2 for stability of explicit method


  // the numerical solving scheme may get unstable when choosing higher number of spatial samples
  // -> counteract with smaller dt

  // define the potential function:
  typedef std::complex<double> Complex;
  Complex i(0,1);       // imaginary unit
  std::function<Complex (double x)> V;
  V = [=](double x)  { return 0.5 * w*w * x*x; }; // from 10.13
  // use 10.13 together with 10.15 for the initial state
  // old: return c * 0.5*m*w*w * xs*xs;  // https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator


  // create the initial wavefunction:
  double f = 0;
  std::vector<double>  x(Nx);
  std::vector<Complex> Psi_0(Nx);
  RAPT::rsArrayTools::fillWithRangeLinear(&x[0], Nx, xMin, xMax);
  for(int k = 0; k < Nx; k++)  {
    double gauss = exp(-0.5 * w * x[k]*x[k]);  // Eq 10.15
    Psi_0[k] = gauss;

    //double gauss = exp(-(x[k]-mu)*(x[k]-mu) / (sigma*sigma));
    //Psi_0[k]  = 0.125 * (1.0+i) * gauss;
    //Psi_0[k] *= rsPolynomial<double>::evaluateHermite(scl * x[k], E); // energy level
    //Psi_0[k] *= exp(i*(2.*PI*f*x[k]/(xMax-xMin)));                    // initial momentum/velocity
  }
  // todo: 
  //  -maybe multiply by +i or -i or maybe exp(i*phi) for phi being an arbitrary angle
  //  -is this how we give it an initial velocity? nope! but how do we? i think, we need to shift 
  //   the Fourier trafe to a frequency other than 0 - by multiplying it with exp(i*lamda*x)? - we 
  //   probably should not multiply by a real sinusoid because that would give a symmetric 
  //   magnitude spectrum for the momentum
  // -maybe multiply by a Hermite polynomial this gives the eigenfunctions of the harmonic 
  //  oscillator - the order of the polynomial is the energy level 
  //    https://en.wikipedia.org/wiki/Hermite_polynomials
  // maybe normalize to unit mean
  // see also here:
  // https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation#One-dimensional_examples
  // it says, what scl should actually be


  // create and set up the particle object:
  rsQuantumParticle<double> p;
  p.setAxisX(x);
  p.initializeWaveFunction(Psi_0);
  p.setPotentialFunction(V);



  int numSamples = numCycles * Nx;
  std::vector<Complex> xOut;
  xOut.reserve(numSamples);
  // simulate the time evolution of the wavefunction:
  int recordedCycles = 0;
  int iterations = 0;
  while(recordedCycles < numCycles)
  {
    if(iterations % skipRatio == 0) {
      GNUPlotter::plotComplexArrayReIm(&x[0], &(p.getWaveFunction())[0], Nx); // ugly syntax
      RAPT::rsAppend(xOut, p.getWaveFunction());
      recordedCycles++;
    }
    p.updateWaveFunction(dt);
    iterations++;
  }
  normalize(&xOut[0], (int) xOut.size());
  writeToWaveFile("SchroedingerParticle.wav", xOut, 44100);
  std::cout << "Done";
  
 
  /*
  int plotInterval = 500;
  int Nt = 100000;
  for(int n = 0; n < Nt; n++)
  {
    if(n % plotInterval == 0)
      GNUPlotter::plotComplexArrayReIm(&x[0], &(p.getWaveFunction())[0], Nx); // ugly syntax
    p.updateWaveFunction(dt);
  }
  */


  // Ovservations:
  // -k controls the "LFO" frequency
  // -sigma controls the initial width - and the shape of the amp-envelope is different for each 
  //  sigma
  // -it sounds like a wah-wah - increasing sigma turns up the "resonance"
  // -maybe the amp-envelope is more interesting as waveform than the signal itself?
  // -maybe try other initial shapes
  //  -gaussian bell multiplied by Hermite polynomials (these are the other eigenfunctions to the
  //   different energy levels)
  //  -anything
  // -with increasing sigma, the signal gets louder - we may need to normalize somehow
  //  -it also seems get more unstable
  //

  // Ideas: 
  // -stabilize the algorithm by filtering out the Nyquist freq after each iteration and/or
  //  normalizing the total probability over the whole buffer
  // -the filter could be a two-sample moving average
  // -see: https://en.wikipedia.org/wiki/Finite_difference_method#Explicit_method - it says,
  //  it's stable for: r := dt / dx^2 <= 1/2  ...at leats in the cae of the ehat equation - but
  //  maybe that applies more generally?
  // -maybe try the implicit and the crank-nicholson method (tridiagonal system?)
  // -ok, with Nx = 128, dt = 5.e-5, the explicit method seems to be stable, too - nevertheless
  //  -> try crank-nicolson - mayb we can choose a bigger time-step with it and get an overall
  //  faster algorithm -  oh no - it just stays stable for longer but not indefinitely

  //GNUPlotter::plotComplexArrayReIm(&x[0], &Psi_0[0], Nx);
  //GNUPlotter plt;
  //plt.addDataArrays(Nx, &x[0], 

  int dummy = 0;
}


void tennisRacket()
{
  // Numerically integrates the system of differential equations describing the angular velocities 
  // of a rigid object about its 3 principal axes of rotation. Rotation around the axis with the 
  // intermediate moment of inertia is unstable and flips over periodically whenever there is a 
  // small amount of initial angular velocity along any of the other two axes. It's an unstable
  // equilibrium of the Euler equations.
  // see:
  // https://en.wikipedia.org/wiki/Tennis_racket_theorem
  // https://en.wikipedia.org/wiki/Euler%27s_equations_(rigid_body_dynamics)
  // https://www.youtube.com/watch?v=1VPfZ_XzisU
  // https://arxiv.org/pdf/1606.08237.pdf

  // User parameters:
  int N = 5000;       // number of samples
  double h = 0.01;    // step-size ("delta-t")
  double I1, I2, I3;  // the 3 moments inertia
  double c = 2;       // ratio of I1/I2 and I2/I3
  I1 =  c;
  I2 =  1;
  I3 =  1/c;
  double w1, w2, w3;  // the 3 (initial) angular velocities (with respect to the principal axes)
  w1 = 0.01;
  w2 = 1;
  w3 = 0.0;
  bool renormalize = true; // renormalize rotational energy in each step (counteract numeric drift)

  // coefficients for the creative extra terms:
  double p = 0.0;   // a small negative number (-0.01) leads to some asymmetry (and decay, which 
                    // can probably be compensated by enforcing constant rotational energy and/or
                    // angular momentum)
  // move this into a 2nd implementation with more creative add-ons

  // create time axis and vectors to hold the results (w1,w2,w3 as functions of time):
  std::vector<double> t(N), W1(N), W2(N), W3(N);
  double a1, a2, a3;  // the 3 angular accelerations
  double E, E0;       // rotational energy and its initial value
  double k1, k2, k3;  // multipliers resulting from the moments of inertia
  k1 = (I2-I3)/I1;
  k2 = (I3-I1)/I2;
  k3 = (I1-I2)/I3;
  E0 = (I1*w1*w1 + I2*w2*w2 + I3*w3*w3)/2; // https://en.wikipedia.org/wiki/Moment_of_inertia#Kinetic_energy_2

  // Use forward Euler method to integrate the ODE system of the Euler rotation equations 
  // (...Euler is every-fucking-where):
  for(int i = 0; i < N; i++)
  {
    // compute absolute time and record angular velocities:
    t[i]  = i*h;
    W1[i] = w1;
    W2[i] = w2;
    W3[i] = w3;

    // compute angular accelerations:
    a1 = k1*w2*w3;
    a2 = k2*w3*w1 + p*w2;
    a3 = k3*w1*w2;

    // update angular velocities:
    w1 += h*a1;
    w2 += h*a2;
    w3 += h*a3;

    // optionally renormalize rotational energy:
    E = (I1*w1*w1 + I2*w2*w2 + I3*w3*w3)/2;
    if(renormalize) {
      double r = sqrt(E0/E);
      w1 *= r;
      w2 *= r;
      w3 *= r;
    }
  }

  // Plot w1,w2,w3 as functions of time:
  GNUPlotter plt;
  plt.addDataArrays(N, &t[0], &W1[0], &W2[0], &W3[0]);
  plt.plot();

  //rosic::writeToMonoWaveFile("TennisRacket.wav", &W2[0], N, 44100);

  // Observations:
  // -if there's a small initial w1 or w3 component, the instability lets the w2 component (blue) 
  //  flip back and forth periodically
  // -the excursions of w3 (green) are greater than those of w1 (black)
  // -the difference in the strength of these excursions depends of the ratio of I1 and I3
  // -in case of an initial w1 component, the w1 excursions are always positive and the w3 
  //  excursions alternate
  // -in case of an initial w3 component, the w1 excursions alternate and the w3 excursions are 
  //  always positive
  // -choosing I1=4,I2=5,I3=1 (i.e. the middle I is greatest), we get a stable rotation, as the 
  //  theory predicts (w2 stays constant) 
  //  -there is some (sinusoidal?) wobble between w1 and w3 - they seem to exchange energy
  // -choosing I1=4,I2=0.5,I3=1 (i.e. the middle I is smallest), w2 wobbles a little bit, w1 and w3
  //  also wobble but with unequal amounts (w3 wobbles more)
  // -i guess, the instability around the intermediate axis is due to k2 being negative while k1,k3
  //  are positive?
  // -using I=(8,4,1) instead of I=(4,2,1) results in an increase of frequency of the flips
  //  -guess: is the frequency proportional to the magnitude to the vector valued moment of inertia 
  //   sqrt(I1^2 + I2^2 + I3^2)
  // -the frequency of the flips seems to go up with the length of the initial angular velocity 
  //  vector but also with the ratio of the initial disturbance to w2 (figure out details)



  // todo:
  // -try to figure out, how the frequency of the flips depends on the various parameters - can we 
  //  find a formula - or even better, a sort of inverse formula that let's us dial in the 
  //  frequency as user parameter?
  //  -nope: I=(4,2,1) has the same frequency as I=(8,4,2), so the absolute numbers seem irrelevant
  //  -I=(16,4,2) has much higher frequency than I=(8,4,2)
  //  -try I = (c,1,1/c) and figure out frequency f as function of c (make plots)
  //   -c=1 -> f=0

  //  ...maybe it depends on the total mass or the norm of the inertia vector?
  // -figure out effects of having initial nonzero values for both, w1 and w3 
  // -figure out effects of the sign(s) of the initial angular velocities
  // -compute angular momentum and rotational energy (as functions of time) - they should remain
  //  constant
  // -figure out the formula for the angular momentum in the fixed "lab" reference frame - maybe 
  //  enforce this to stay constant

  // -according to this https://arxiv.org/pdf/1606.08237.pdf (section 3), the system has as analytic 
  //  solutions the jacobi elliptic functions cn,sn,dn - plot them, too
}

// see also:
// https://en.wikipedia.org/wiki/Moment_of_inertia
// https://en.wikipedia.org/wiki/Poinsot%27s_ellipsoid
// https://en.wikipedia.org/wiki/Polhode

template<class T>
T rotationalEnergy(const rsVector3D<T>& I, const rsVector3D<T>& w)
{
  return T(0.5) * (I.x*w.x*w.x + I.y*w.y*w.y + I.z*w.z*w.z);
  // https://en.wikipedia.org/wiki/Moment_of_inertia#Kinetic_energy_2
}

template<class T>
void plotVectorComponents(std::vector<T>& t, std::vector<rsVector3D<T>>& v)
{
  int N = (int) t.size();
  rsAssert((int)v.size() == N);
  std::vector<T> x(N), y(N), z(N);
  for(int n = 0; n < N; n++) {
    x[n] = v[n].x;
    y[n] = v[n].y;
    z[n] = v[n].z;
  }
  GNUPlotter plt;
  plt.addDataArrays(N, &t[0], &x[0], &y[0], &z[0]);
  plt.setPixelSize(1000, 300);
  plt.plot();
}
// move to rapt, make input vectors const (requires changes to GNUPlotter)


void tennisRacket2()
{
  // In contrast to tennisRacket, we use vectors here and introduce some additional tweaks in an 
  // attempt to provide a musically viable set of user parameters so this thing can be used as a 
  // waveform generator or resonator.

  // User parameters:
  typedef rsVector3D<double> Vec;
  int    N = 8000;          // number of samples
  double h = 0.01;          // step-size ("delta-t")
  Vec I(   4, 2, 1);        // moments of inertia along principa axes
  Vec w(0.01, 1, 0);        // initial angular velocities
  bool renormalize = false; // keep energy constant

  // experimental parameters: 
  Vec d(0.0, 0.02, 0.0); // damping/dissipation - works only when renormalization is off


  // integrate ODE system via forward Euler method :
  Vec a;                     // angular acceleration
  Vec M(0, 0, 0);            // external torque vector (zero at the moment)
  std::vector<double> t(N);  // time axis
  std::vector<Vec>    W(N);  // vectorial angular velocity as function of time
  double E;                  // rotational energy...
  double E0 = rotationalEnergy(I, w); //...and its initial value
  for(int i = 0; i < N; i++)
  {
    // compute absolute time and record angular velocity vector:
    t[i] = i*h;
    W[i] = w;

    // compute angular acceleration vector:
    a.x = (M.x - (I.z - I.y)*w.y*w.z) / I.x;
    a.y = (M.y - (I.x - I.z)*w.z*w.x) / I.y;
    a.z = (M.z - (I.y - I.x)*w.x*w.y) / I.z;
    // formula from https://en.wikipedia.org/wiki/Euler%27s_equations_(rigid_body_dynamics)

    // todo: add experimental extra terms to the angular acceleration...
    a.x -= d.x*w.x / I.x;
    a.y -= d.y*w.y / I.y;
    a.z -= d.z*w.z / I.z;


    // update angular velocity vector:
    w += h*a;

    // optionally renormalize rotational energy:
    E = rotationalEnergy(I, w);
    if(renormalize)
      w *= sqrt(E0/E);
  }

  // -when only the 1st axis is damped, it tends towards rotating only around the 3rd axis
  // -when only the 3rd axis is damped, it tends towards rotating only around the 1st axis
  // -when only the 2nd axis is damped, it dies out

  // todo: 
  // -try asymmetrical damping that affects the both half-waves of the "square-wave" differently
  //  ->goal: some sort of pulse-width control

  plotVectorComponents(t, W);
}

void tennisRacket3()
{
  // This produces the same plot as the functions above but uses the rsTennisRacket class

  int N = 8000; 
  double h = 0.01;
  rsTennisRacket<double> tr;
  tr.setInertiaRatio(2);
  tr.setState(0.01, 1, 0);  // initial state
  tr.setStepSize(h);

  std::vector<double> t(N), w1(N), w2(N), w3(N);
  for(int n = 0; n < N; n++) {
    t[n]  = n*h;
    w1[n] = tr.getW1();
    w2[n] = tr.getW2();
    w3[n] = tr.getW3();
    tr.updateState(0, 0, 0);
  }

  GNUPlotter plt;
  plt.plotFunctionTables(N, &t[0], &w1[0], &w2[0], &w3[0]);
}

// todo: make a function tennisRacketFreqs that plots the frequency or period of the flipping
// against the inertia ratio - goal: figure out how the inertia ratio controls the freq

double tennisRacketPeriod(rsTennisRacket<double>& tr)
{
  // todo: estimate period of the tennis racket by measuring the distance between two upward
  // zero crossings...

  int first = -1, second = -1;  // sample instants of zero crossings
  double x1 = tr.getSample(0);
  int n = 0;
  while(second == -1) {
    double x  = tr.getSample(0);
    if(x1 <= 0 && x > 0) // zero crossing found
    {
      if(first == -1)
        first = n;
      else
        second = n;
    }
    x1 = x;
    n++;
  }
  return second - first;

  // we could refine this to estimate zeros with subsample precision - but that may be overkill 
  // here
}

void tennisRacketFreq()
{
  // We plot the frequency of the flipping as function of the inertia ratio


  rsTennisRacket<double> tr;
  tr.setStepSize(0.01);

  std::vector<double> r = { 1.01, 1.125,1.25,1.5,1.75,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20 };
  //std::vector<double> r = { 2,4,8,16,32,64 };
  int N = (int) r.size();
  std::vector<double> p(N), f(N); // period and frequency

  for(int i = 0; i < N; i++)
  {
    tr.setInertiaRatio(r[i]);
    tr.setState(0.01, 1, 0);  // initial state
    p[i] = tennisRacketPeriod(tr);
    f[i] = 1/p[i];
  }

  GNUPlotter plt;
  //plt.plotFunctionTables(N, &r[0], &p[0]);
  plt.plotFunctionTables(N, &r[0], &f[0]);

  // It looks similar to a straight line, but not quite...
  // maybe some variant of this?
  // https://en.wikipedia.org/wiki/Logarithmic_integral_function
  // https://de.wikipedia.org/wiki/Integrallogarithmus#/media/Datei:Li10.png

  // -maybe plot the numerical derivative
  // -maybe use the logarithm of the ratio for the x-axis - let it extend also to negative values,
  //  i.e. r < 1 - maybe it's symmetrical?
}