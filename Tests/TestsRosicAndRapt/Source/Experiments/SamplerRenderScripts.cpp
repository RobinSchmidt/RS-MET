

// get rid - there's a copy of it now in RenderScriptTools.cpp which should be used
std::vector<double> randomizePhases(const std::vector<double>& x, int seed, double amount)
{
  int N = (int)x.size();
  RAPT::rsAssert(RAPT::rsIsPowerOfTwo(N), "Works only for power of 2 lengths");

  // Do forward transform:
  std::vector<double> y(N), a(N), p(N);
  using FT = RAPT::rsFourierTransformerRadix2<double>;   // todo: use Bluestein later
  FT ft;
  ft.setBlockSize(N);
  ft.getRealSignalMagnitudesAndPhases(&x[0], &a[0], &p[0]);

  // Randomize phases:
  RAPT::rsNoiseGenerator<double> ng;
  ng.setSeed(seed);
  ng.setRange(0.0, amount * 2*PI);
  for(int k = 1; k < N; k++) {    // DC is unaffected
    p[k] += ng.getSample();
    p[k] =  RAPT::rsWrapToInterval(p[k], -PI, +PI); }

  // Do inverse transform and return result:
  ft.getRealSignalFromMagnitudesAndPhases(&a[0], &p[0], &y[0]);
  return y;
}

void createBassdrumPsy1Sample(double freqScale = 1.0, bool plot = false)
{
  // Create a bassdrum sample using a weighted sum of exponential envelopes with different decay
  // times for a sine-sweepdown. Create also overtones at twice and thrice the frequency, maybe let
  // them have an attack/decay and an attack envelope respectively

  // Rendering parameters:
  int fs = 44100;
  //int N  = 44100*2;
  int N = rsPowInt(2, 16); // maybe have special function rsPowerOf2
  //bool plot = false;  // if true, we plot stuff, if false, we render a wavefile

  // Frequency envelope parameters:
  double freqDecay1  =   10;    // in ms
  double freqDecay2  =   80;
  double freqDecay3  =  240;
  double freqWeight1 =  6.0;
  double freqWeight2 = -1.0;
  double freqWeight3 =  1.0;
  double freqFloor   =  0.0;
  double freqCeil    =  800;  // an excursion/depth would be better
                              //double freqScale   = 1.0;    // 1: fundamental, 2,3,4,etc: overtones

  // Amplitude envelope parameters:
  double ampDecay1  =   20;
  double ampDecay2  =  100;
  double ampDecay3  =  400;
  double ampWeight1 =  0.75;
  double ampWeight2 = -1.0;
  double ampWeight3 =  1.0;


  int plotDecimate = 32;
  if(plot)
  {
    fs /= plotDecimate;
    N  /= plotDecimate;
  }

  // Create frequency envelope:
  using AT = RAPT::rsArrayTools;
  RAPT::rsAttackDecayEnvelope<double> eg1, eg2, eg3; // a simple decay env would suffice

  eg1.setAttackSamples(0);
  eg2.setAttackSamples(0);
  eg3.setAttackSamples(0);
  eg1.setDecaySamples(freqDecay1 * 0.001 * fs);
  eg2.setDecaySamples(freqDecay2 * 0.001 * fs);
  eg3.setDecaySamples(freqDecay3 * 0.001 * fs);
  using Vec = std::vector<double>;
  Vec env1(N), env2(N), env3(N), env(N);
  eg1.noteOn(60, 100);
  eg2.noteOn(60, 100);
  eg3.noteOn(60, 100);
  for(int n = 0; n < N; n++)
  {
    env1[n] = eg1.getSample();
    env2[n] = eg2.getSample();
    env3[n] = eg3.getSample();
    env[n]  =  freqWeight1*env1[n] + freqWeight2*env2[n] + freqWeight3*env3[n];
  }
  if(plot) rsPlotVectors(env1, env2, env3, env);

  // Create the raw sine-sweep:
  AT::normalize(&env[0], N);  // makes it inconvenient to port to realtime
  Vec xL(N), xR(N);
  double phi = 0;  // phase
  for(int n = 0; n < N; n++)
  {
    xL[n] = sin(phi - PI/4);
    xR[n] = sin(phi + PI/4);
    double f = freqFloor + (freqCeil - freqFloor) * env[n];
    double w = freqScale * 2*PI*f/fs;
    phi += w;
  }
  if(plot) rsPlotVectors(xL, xR);

  // Create amplitude envelope:
  eg1.setDecaySamples(ampDecay1 * 0.001 * fs);
  eg2.setDecaySamples(ampDecay2 * 0.001 * fs);
  eg3.setDecaySamples(ampDecay3 * 0.001 * fs);
  eg1.reset(); eg1.noteOn(60, 100);
  eg2.reset(); eg2.noteOn(60, 100);
  eg3.reset(); eg3.noteOn(60, 100);
  for(int n = 0; n < N; n++)
  {
    env1[n] = eg1.getSample();
    env2[n] = eg2.getSample();
    env3[n] = eg3.getSample();
    env[n]  =  ampWeight1*env1[n] + ampWeight2*env2[n] + ampWeight3*env3[n];
  }
  if(plot) rsPlotVectors(env1, env2, env3, env);

  // Apply amplitude envelope and fade-out:
  AT::normalize(&env[0], N);
  for(int n = 0; n < N; n++)
  {
    xL[n] *= env[n];
    xR[n] *= env[n];
  }
  RAPT::rsFadeOut(&xL[0], N-N/8, N-1);
  RAPT::rsFadeOut(&xR[0], N-N/8, N-1);

  // Create ambience sample:
  Vec amb = randomizePhases(xL+xR, 2, 1.0);
  RAPT::rsArrayTools::normalize(&amb[0], N, 1.0);

  // Plot final result or write to wvaefile: 
  if(plot)  rsPlotVectors(xL, xR);
  if(!plot)
  {
    std::string fileName = "BassdrumPsy1";
    if(freqScale == 1.0) fileName += "_Prime";
    if(freqScale == 1.5) fileName += "_Fifth1";
    if(freqScale == 2.0) fileName += "_Octave1";
    if(freqScale == 3.0) fileName += "_Fifth2";
    if(freqScale == 4.0) fileName += "_Octave2";
    std::string nameWithExt = fileName + ".wav";
    rosic::writeToStereoWaveFile(nameWithExt.c_str(), &xL[0], &xR[0], N, (int)fs);

    fileName += "_Ambience";
    nameWithExt = fileName + ".wav";
    rosic::writeToMonoWaveFile(nameWithExt.c_str(), &amb[0], N, (int)fs);
  }

  // ToDo:
  // -Write this as a jupyter notebook - we really want quck REPL evaluation for this
  // -Maybe write this as an APE script
  //  -Give the user just a single choice parameter to select the preset, settings are encoded in 
  //   the code...but maybe give some macro-parameters to the user
  //  -but APE does not yet support midi - maybe trigger the drum with an input impulse
  // -Have an envelop for the waveshape. it should control the amount of feedback PM:
  //  y[n] = sin(phi + fb*y[n])  -> nonlinear -> needs implicit solver maybe using y[n-1] as 
  //  initial guess, Newton iteration: f(y) = y - sin(phi + fb*y), f'(y) = 1 - fb*cos(phi + fb*y)
  //  -maybe have also PM by a 2nd independent signal...maybe a sine at twice the 
  //  freq? I want something that turns a sine into a square - the self-feedback leads to a 
  //  sawtooth shape...or maybe make it squarish by waveshaping maybe use something like
  //  tanh(d*x) / d except when d == 0 which is treated as special case returning just x...have a
  //  drive envelope for d
  // -create a difference between a signal with and without the waveshaping and phasemodulation 
  //  applied - the difference gives us an "overtones" signal that we can mix with the fundamental
  //  signal to taste
  // -Maybe remove the normalization of the freq-env. it leads to the effect that increasing 
  //  freqWeight1 increases the overall freq but we want it to affect only the initial freq
  // -Maybe produce stereo output by using +-45° phase shift for left/right
  // -Mix in a very short noise-burst to the attack (might have settings for: decay, HP, LP, color, seed)
  //  -maybe it should have a filter env
  //  -maybe clip or saturate sweep+noise
  // -Apply reverb - maybe using a phase-randomization algorithm: do one big FFT (maybe 
  //  length = 2^16 = 65536), randomize phases, iFFT ...this need its own amp-env
  // -Add an envelope generator to ToolChain that implements a weighted sum of exponential decays
  // -Maybe give ToolChain a rendering functionality

  // -try to raise the amp-env to a time-varying power - with 2 , it should look more like Gaussian
  // -apply distortion, maybe bitcrushing in a time variant manner to add some noisiness to the
  //  transient
  // -for such effects, obtain a difference signal to the original and store the pure effect 
  //  seperately to be mixed in later
  // -use a more sawtooth-like waveform to give it some overtones...however, maybe the 1st and 2nd
  //  should be notched out because 100 or 150 hz are ugly in bassdrums. or: maybe also subtract
  //  from a sample using a sine, obtain difference and pass that through a highpass before mixing
  // -the (gaussian) envelope should be applied after reverb/convolution -> sort of gated reverb


  // ToDo:
  // -render some goasque creaking sounds and make an sfz (use pan-envelopes)...but how? maybe using a 
  //  percussive sound with echo using a delay-length corresponding to a note and a lot of feedback?
  // -render some multiplicative synthesis lead sounds
  // -render some swooshes and reversed sounds ...maybe a noise with filter env
  // -create a bouncy goa-bass that harmonizes well with this bassdrum...maby they need to also to 
  //  the one_shot thing? but the we must take care of interference...maybe only the amp-envs 
  //  should get added on to of each other - the oscs just keep playing
  // -maybe wrap bassdrum and bass into a single instrument using keyranges - the wen can also 
  //  apply some dynamics processing to their mix (if bus_mode is active)
  // -maybe convolve it with an expoentially decaying white noise sound, 
  //  -maybe that noise should have a color envelope from white to brown, say


  // -For dubstep sounds, try to use a percussive sample like a snare, set up a loop and modify its 
  //  location via midi-cc, maybe have a second one an octave below going on. maybe change the loop 
  //  length (while adjusting the increment accordingly), maybe also when moving loop_start and 
  //  loop_end simultaneously, adapt the increment to account for the fact that the loop point is
  //  receding or approaching ...but maybe not
  // -it has often very formantish stuff going on, try also comb-filters and sync-stuff
  // -LFOs are important for modulating the position - their frequency should be midi controlled
  // -use MSEGs, attacks should be normal, decays should use "anti-analog" shape
  // -speed should speed up or down in step of simple ratios (like doubling, tripling, halving, ..)
  // -distortion and bitcrushing also helps
  // -Maybe a multiplicative synthesis sample could also be used as basis
  // https://www.youtube.com/watch?v=dknbNrr4EDo&list=RDQM6d8Vubj1zMg&start_radio=1


  int dummy = 0;
}

void createSweepDrummerSamples(int sampleRate, double freqScale)
{
  // We do the same as the function above but this time with less code using the new rsSweepDrummer
  // class. When this is finished, the function above will be obsolete...




  int dummy = 0;
}

// move this to RenderScripts:
void createBassdrumPsy1Samples()
{
  //createBassdrumPsy1Sample(1.0, true); return;  // for development


  bool plot = false;
  createBassdrumPsy1Sample(1.0, plot);
  createBassdrumPsy1Sample(1.5, plot);
  createBassdrumPsy1Sample(2.0, plot);
  createBassdrumPsy1Sample(3.0, plot);
  createBassdrumPsy1Sample(4.0, plot);
}


void createMiscSamples()
{
  // Create miscelanneous other samples that are useful as raw material in the sampler engine.

  using Vec = std::vector<double>;
  using AT  = RAPT::rsArrayTools;

  int fs = 44100;
  int N  = 0;       // number of samples

                    // Generate a unit impulse:
  double one = 1.0;
  rosic::writeToMonoWaveFile("UnitImpulse.wav", &one, 1, fs);
  // -Maybe we should make it a few samples long. Having a wavfile containing just a single value
  //  may be a corner case that some sampler engines won't like? But it actually makes for a nice
  //  unit test to see what an engine does in such an extreme case. OK, my engine handles it 
  //  nicely: it just produces a unit impulse indeed. 
  // -But: we can't create filter blips this way because the region player immediately stops, 
  //  giving the filter no ringtout time. We need a couple of samples of silence after the impulse
  //  and set up a loop over this silence portion.
  //  

  // Generate 5 seconds of white noise with uniform amplitude distribution:
  N = 5*fs;
  Vec x(N);
  RAPT::rsNoiseGenerator<double> ng;
  for(int n = 0; n < N; n++)
    x[n] = ng.getSample();
  rosic::writeToMonoWaveFile("UniformWhiteNoise.wav", &x[0], N, fs);
  // Notes:
  // -Layering this sample with itself with various offsets can be used to obtain Irwin-Hall
  //  distributed white noise.
  // -5 seconds should be long enough to have no noticable repetition pattern and/or 
  //  comb-filtering artifacts when layering several shifted copies.


  int dummy = 0;
}