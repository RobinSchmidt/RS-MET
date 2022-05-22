using namespace rotes;
using namespace rosic;
using namespace RAPT;

void rotes::testOscillatorStereo()
{
  // We test some functionality the stereo wavetable oscillator used in Straightliner. We pass it
  // a sawtooth wave, set up some of the oscillator's waveform modification parameters and retrieve
  // the resulting waveform in its decimated version (which is supposed to be used for plots on the
  // GUI) and plot it.

  // set up the WaveTable:
  static const int prototypeLength = 2048;
  double prototypeWaveFlat[prototypeLength];
  double *prototypeWave[2];
  prototypeWave[0] = &prototypeWaveFlat[0];
  prototypeWave[1] = &prototypeWaveFlat[0];
  int n;
  for(n=0; n<prototypeLength; n++)
    prototypeWaveFlat[n] = RAPT::rsSawWave( 2*PI*n / prototypeLength);
  MipMappedWaveTableStereo waveTable;
  waveTable.setWaveform(prototypeWave, prototypeLength);

  // set up the oscillator:
  OscillatorStereo osc;
  osc.setWaveTableToUse(&waveTable);
  osc.setStereoPhaseShift(90.0);
  osc.setStartPhase(95.0);

  // obtain the data for the plot:
  static const int plotLength = 160;
  double plotIndices[plotLength];
  RAPT::rsArrayTools::fillWithIndex(plotIndices, plotLength);

  double plotDataFlat1[2*plotLength];
  double *plotData1[2];
  plotData1[0] = &plotDataFlat1[0];
  plotData1[1] = &plotDataFlat1[plotLength];
  osc.getWaveformForDisplay(plotData1, plotLength);

  osc.setStartPhase(96.0);
  double plotDataFlat2[2*plotLength];
  double *plotData2[2];
  plotData2[0] = &plotDataFlat2[0];
  plotData2[1] = &plotDataFlat2[plotLength];
  osc.getWaveformForDisplay(plotData2, plotLength);

  plotData(plotLength, plotIndices, plotData1[0], plotData2[0]);

  // Observations:
  // The waveform plot looks a bit ugly, jaggy and also, both channels look almost the same but not
  // quite. Is that supposed to be so?
}


void renderLorenzSoundToFile(int numSamples)
{
  rosic::LorentzSystem lorentzSystem;
  lorentzSystem.setPseudoFrequency(500); // make this a parameter

  int N = numSamples;
  std::vector<double> x(N), y(N), z(N);
  lorentzSystem.getState(&x[0], &y[0], &z[0]);
  for(int n = 1; n < N; n++)
  {
    lorentzSystem.iterateState();
    lorentzSystem.getState(&x[n], &y[n], &z[n]);
  }

  RAPT::rsArrayTools::normalize(&x[0], N, 1.0, true);
  RAPT::rsArrayTools::normalize(&y[0], N, 1.0, true);
  RAPT::rsArrayTools::normalize(&z[0], N, 1.0, true);
  writeToMonoWaveFile("LorenzX.wav", &x[0], N, 44100, 16);
  writeToMonoWaveFile("LorenzY.wav", &y[0], N, 44100, 16);
  writeToMonoWaveFile("LorenzZ.wav", &z[0], N, 44100, 16);
}

void rotes::testLorentzSystem()
{
  // Produces and plots output of the 3 variables x,y,z as functions fo time t  of 
  // rosic::LorentzSystem.

  //renderLorenzSoundToFile(100000);

  static const int N = 2000;

  rosic::LorentzSystem lorentzSystem;
  lorentzSystem.setPseudoFrequency(1000);

  double t[N], x[N], y[N], z[N];
  lorentzSystem.getState(&x[0], &y[0], &z[0]);
  for(int n = 1; n < N; n++)
  {
    lorentzSystem.iterateState();
    lorentzSystem.getState(&x[n], &y[n], &z[n]);
  }

  RAPT::rsArrayTools::fillWithIndex(t, N);
  plotData(N, t, x, y, z);
}

bool rotes::testSnowflake()
{
  bool testResult = true;

  rosic::Snowflake sf;
  sf.setSampleRate(1.0);
  sf.setFrequency(0.35 / 4); // inc = frequency*numLines/sampleRate

  int N = 1000;  // number of samples
  int n;         // sample index
  std::vector<double> xt(N), yt(N), xf(N), yf(N); // table-based and on-the-fly generated outputs

  // create on-the-fly output:
  sf.setUseTable(false);
  sf.reset();
  for(n = 0; n < N; n++)
    sf.getSampleFrameStereo(&xf[n], &yf[n]);

  // create table-based output:
  sf.setUseTable(true);
  sf.reset();
  for(n = 0; n < N; n++)
    sf.getSampleFrameStereo(&xt[n], &yt[n]);

  // xt should equal xf, same for yt,yf - check this and use it as return value
  double xError = 0, yError = 0;
  for(n = 0; n < N; n++)
  {
    xError = RAPT::rsMax(xError, fabs(xt[n]-xf[n]));
    yError = RAPT::rsMax(yError, fabs(yt[n]-yf[n]));
    testResult &= xError == 0 && yError == 0;
    rsAssert(testResult);
  }

  // compare outputs with and without reset in on-the-fly mode (there are some audible buzzy 
  // artifacts when reset is off - it should be equal up to numerical roundoff - but isn't):

  std::vector<double> x0(N), y0(N), x1(N), y1(N); // reset modes 0 (never) and 1 (each cycle)

  // create a unit square:
  sf.clearRules();
  sf.setAxiom("F+F+F+F+");
  sf.setNumIterations(0);
  sf.setTurnAngle(90);
  sf.updateAllInternals();

  sf.setUseTable(false);

  sf.reset();
  sf.setResetRatio(0, 1);  // reset after each cycle
  for(n = 0; n < N; n++) 
    sf.getSampleFrameStereo(&x1[n], &y1[n]);

  sf.reset();
  sf.setResetRatio(0, 0);  // never reset
  for(n = 0; n < N; n++) 
    sf.getSampleFrameStereo(&x0[n], &y0[n]);


  GNUPlotter plt;
  plt.setRange(-1.1, +1.1, -1.1, +1.1);
  plt.setPixelSize(400, 400);
  plt.addCommand("set size square");

  //plt.addDataArrays(N, &xt[0], &yt[0]);
  //plt.addDataArrays(N, &xf[0], &yf[0]);

  plt.addDataArrays(N, &x0[0], &y0[0]);
  plt.addDataArrays(N, &x1[0], &y1[0]);
  plt.plot();


  return testResult;
}

void fillResetInstantArray(rosic::ResetCounter& rc, int* r, int N)
{
  int i = 0;      // counts iterations/ticks
  int j = 0;      // counts resets
  rc.reset();
  while(true) {
    if(rc.tick()) {
      r[j] = i;
      j++;
      if(j == N)
        break;
    }
    i++;
  }
}

bool rotes::testResetter()
{
  bool result = true;

  static const int N = 15;
  int r[N];   // reset instants
  rosic::ResetCounter rc;

  rc.setInterval(5.0);
  fillResetInstantArray(rc, r, N);

  // check, if resets occur at expected times:
  int i;
  for(i = 0; i < N; i++)
  {

  }

  rc.setInterval(5.5);
  fillResetInstantArray(rc, r, N);
  // difference between entries should alternate between 5 and 6: 5,6,5,6,5,6,...

  rc.setInterval(5.25);  
  fillResetInstantArray(rc, r, N);
  // difference between entries tould follow the pattern: 5,5,5,6,5,5,5,6,...

  return result;
}


void rotes::testTurtleReverse()
{
  std::string commands = "F+F+F+F+F+F+F-F-F-F+F+F-F+F+F-F";

  TurtleGraphics tg;
  tg.setAngle(30);
  tg.init();

  std::vector<double> x, y;

  // forward part:
  x.push_back(tg.getX());
  y.push_back(tg.getY());
  size_t i = 0;
  while(i < commands.size()) {
    if(tg.interpretCharacter(commands[i])) {
      x.push_back(tg.getX());
      y.push_back(tg.getY()); }
    i++; }

  // reversed part:
  tg.setReverseMode(true);
  i = (int) commands.size()-1;
  while(i >= 0) {
    if(tg.interpretCharacter(commands[i])) {
      x.push_back(tg.getX());
      y.push_back(tg.getY()); }
    i--; }

  // plot:
  int N = (int) x.size();
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0]);
  plt.addDataArrays(N, &y[0]);
  plt.plot();
}


void setupTurtleSource(rosic::TurtleSource& ts)
{
  // used to set up a TurtleSource to specific settings (the idea is to call it for two objects to
  // make sure, they have the same settings)

  //ts.setTurtleCommands("F");
  //ts.setTurnAngle(0);

  //ts.setTurtleCommands("F+");   // just forward/backward
  //ts.setTurnAngle(180);

  //ts.setTurtleCommands("F+F+");   // just forward/backward
  //ts.setTurnAngle(180);

  //ts.setTurtleCommands("F+F+F+");   // triangle
  //ts.setTurnAngle(120);

  //ts.setTurtleCommands("F+F+F+F+");   // square
  //ts.setTurnAngle(90);

  //ts.setTurtleCommands("F+F+F+F+F+");  // pentagon
  //ts.setTurnAngle(72);


  ts.setTurtleCommands("F+F+F+F+F+F+F+F+");   // octagon
  ts.setTurnAngle(45);


  //ts.setTurtleCommands("FF");
  //ts.setTurnAngle(0);
  //// a simple "F" doesn't work right in forward mode, "FF" works right in forward but not reverse 
  //// mode

  ts.setSampleRate(1);
  //ts.setFrequency(0.05);
  ts.setFrequency(0.125);
  //ts.setFrequency(0.25);
  //ts.setFrequency(0.5);
  //ts.setFrequency(1.0);
  ts.setResetRatio(0, 0); // no reset
  //ts.setResetRatio(0, 2); // reset after 2 cycles
}

bool runTurtleTest(rosic::TurtleSource& ts1, rosic::TurtleSource& ts2, int numSamples, 
  bool plot = false)
{
  // compares the output of two TurtleSource objects and returns true, if they are almost equal

  int N = numSamples;
  std::vector<double> x1(N), y1(N), x2(N), y2(N);
  int n;
  double err;
  double tol = 1.e-14;
  bool result = true;

  for(n = 0; n < N; n++)
  {
    ts1.getSampleFrameStereo(&x1[n], &y1[n]);
    ts2.getSampleFrameStereo(&x2[n], &y2[n]);
    err = RAPT::rsMax(fabs(x2[n]-x1[n]), fabs(y2[n]-y1[n]));
    //rsAssert(ts2.checkIndexConsistency());

    result &= err <= tol;
    //rsAssert(result);
  }

  if(plot)
  {
    GNUPlotter plt;

    //plt.addDataArrays(N, &x1[0], &y1[0]);
    //plt.addDataArrays(N, &x2[0], &y2[0]);
    //plt.setRange(-2, 2, -2, 2);

    plt.addDataArrays(N, &x1[0]);
    plt.addDataArrays(N, &y1[0]);
    plt.addDataArrays(N, &x2[0]);
    plt.addDataArrays(N, &y2[0]);

    plt.plot();
  }

  return result;
}

void rotes::testTurtleSource()
{
  int N = 100; // number of samples
  bool result = true;

  // set up two TurtleSource objects with the same settings, the only difference being that one 
  // uses the table and the other one doesn't:
  rosic::TurtleSource ts1, ts2;
  ts1.setUseTable(true);

  setupTurtleSource(ts1);
  setupTurtleSource(ts2);

  ts1.reset(); result &= ts1.isInInitialState();
  ts2.reset(); result &= ts2.isInInitialState();

  result &= runTurtleTest(ts1, ts2, N, true);
  ts1.setFrequencyScaler(-1.0);
  ts2.setFrequencyScaler(-1.0);
  result &= runTurtleTest(ts1, ts2, N, true);

  ts1.reset(); result &= ts1.isInInitialState();
  ts2.reset(); result &= ts2.isInInitialState();
  result &= runTurtleTest(ts1, ts2, N, true);
    // the result ist shifted (i think, only the x-values - maybe bcs dy is initially 0?)
    // maybe in resetTurtle(), we need to init the turtle differently when in reverse mode

  double reverseRatio = 1.0;
  ts1.setReverseRatio(0, reverseRatio);
  ts2.setReverseRatio(0, reverseRatio);
  ts1.reset(); result &= ts1.isInInitialState();
  ts2.reset(); result &= ts2.isInInitialState();
  result &= runTurtleTest(ts1, ts2, N, true);

}

void samplerPatchTest_BandpassSaw()
{
  // Create the sfz-strings defining the instruments. Maybe factor out into separate functions - 
  // one for each patch - maybe as static functions in a class. Maybe this example patch stuff 
  // should go into an Experiment. Or maybe we should create repo for example patches and this repo
  // could also contain some unit tests. Maybe it should contain some example sfz patches along 
  // with their expected output. Maybe have two repos - one for the patches and one for the tests.
  // Or maybe the tests could go into the existing RS-MET-Tests repo
  // see also here: https://github.com/sfz/tests/

  // A Sawtooth with resonant lowpass and filter envelope:
  std::string sfz1 = "\
<group>\n\
<region>\n\
sample=Saw2048.wav\n\
loop_start=0 loop_end=2048 loop_mode=loop_continuous pitch_keycenter=21\n\
cutoff=2000 resonance=15 fil_type=lpf_2p\n\
fileg_attack=0.2 fileg_decay=0.4 fileg_sustain=0.5 fileg_release=0.5 fileg_depth=600\n\
volume=-10\n\
";


  std::string sfz2 = "\
<group>\n\
<region>\n\
sample=Saw2048.wav\n\
loop_start=0 loop_end=2048 loop_mode=loop_continuous pitch_keycenter=21\n\
cutoff=200 resonance=5 fil_type=hpf_2p\n\
cutoff2=1500 resonance2=15 fil2_type=lpf_2p\n\
fileg_attack=0.2 fileg_decay=0.4 fileg_sustain=50 fileg_release=0.2 fileg_depth=1200\n\
fillfo_freq=7 fillfo_depth=300\n\
volume=-15\n\
ampeg_attack=0.1 ampeg_decay=0.5 ampeg_sustain=50 ampeg_release=0.4\n\
amplfo_freq=11 amplfo_depth=3\n\
";
  // there are two amplifiers in this patch in the se - i think, one should be sufficient
  // todo: add an fillfo, amplfo

  // A sawtooth with amp env using shape parameters:
  std::string sfz3 = "\
<group>\n\
<region>\n\
sample=Saw2048.wav\n\
loop_start=0 loop_end=2048 loop_mode=loop_continuous pitch_keycenter=21\n\
ampeg_attack=0.1 ampeg_decay=0.5 ampeg_sustain=50 ampeg_release=0.4\n\
ampeg_attack_shape=0.5\n\
";
  // triggers assert due to ampeg_attack_shape

  // Create the playback data:
  float fs = 44100;
  int   N  = 90000;
  int   v  = 64;    // velocity

  using Note = rsTestNoteEvent;
  using NoteList = std::vector<Note>;
  NoteList notes = { Note{45, v, 0, 50000} };
  //NoteList notes = { Note{45, v, 0, 50000},  Note{52, v, 10000, 50000} };

  // Create a sampler engine, set it up from the sfz string and let it produce the output according
  // to our sequence of notes
  rosic::Sampler::rsSamplerEngine2 se;
  se.setSampleRate(fs);
  se.setFromSFZ(sfz2);
  using Vec = std::vector<float>;
  Vec outL(N), outR(N);
  getSamplerNotes(&se, notes, outL, outR);
  //rsPlotVectors(outL, outR);

  rosic::writeToStereoWaveFile("BandpassSaw1.wav", &outL[0], &outR[0], N, (int)fs, 16);


  // ToDo:
  // -use the Saw2048 sample, apply a lowpass and a highpass (both 1st order) with an ADSR 
  //  envelope, maybe use keytracking at least for the highpass, maybe both
  // -Make a patch featuring 3 eq bands at 500, 1000, 2000 Hz, the outer ones narrow dips, the 
  //  inner a broad peak, Modulate all frequencies by an LFO. Emulates smallstone phaser.
  //  Maybe use a stereo-shift
  // -Make another similar patch but this time with the phaser on the group level in a mode where
  //  the group acts as a sub-bus
}


void rotes::testSamplerEngine()
{

  generateTestSamples();
  samplerPatchTest_BandpassSaw();


  int dummy = 0;
}