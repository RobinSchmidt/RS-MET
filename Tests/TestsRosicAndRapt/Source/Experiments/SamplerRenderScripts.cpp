
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

std::vector<double> renderAlternatingSquareSweep(int N=2048)
{
  // Renders the sequence 1,-1,1,1,-1,-1,1,1,1,1,-1,-1,-1,-1,....
  // Its intended to be used as transient for drum sounds

  RAPT::rsAssert(RAPT::rsIsPowerOfTwo(N), "Works only for power of 2 lengths");

  std::vector<double> y(N);
  int M = 2;  // current sub-length
  int n = 0;  // current index in sub-sequence
  int s = 0;  // current start
  while(M < N)
  {
    n = 0;

    while(n < M/2)
    {
      y[s+n] = 1;
      n++;
    }
    while(n < M)
    {
      y[s+n] = -1;
      n++;
    }

    s += M;
    M *= 2;
  }

  return y;
}

void renderSweepBassdrums(int sampleRate)
{
  // We do the same as the function above but this time with less code using the new rsSweepDrummer
  // class. When this is finished, the function above will be obsolete...

  /*
  // test to render some transient to mix in:
  int N = 2048;
  auto test = renderAlternatingSquareSweep(N);
  rsPlotVector(test);
  rosic::writeToMonoWaveFile("SquareSweep2048.wav", &test[0], N, sampleRate);
  // using this as input to a reverb and then cutting off the initial section and further shaping
  // the response can give intersting sounds. I think, we should generally explore the possibility
  // to render drumsounds via FDNs
  */


  rsSweepDrummer<double> sd;
  sd.setSampleRate(sampleRate);
  //sd.setFreqScale(1.0);
  //sd.setTimeScale(1.0);
  //sd.setFreqFloor(0.0);
  //sd.setFreqDepth(800);

  // Helper function to apply time variying quantization (experimental - doesn't seem to sound 
  // good):
  auto quantize = [](std::vector<double>& xL, std::vector<double>& xR, 
    double amount, double attack, double decay)
  {
    int N = (int) xL.size();
    for(int n = 0; n < N; n++)
    {
      xL[n] = rsQuant(xL[n], amount);
      xR[n] = rsQuant(xR[n], amount);
    }
  };


  // Helper function to render output and write to file:
  auto render = [&](const char* fileName, double length, double fadeLength)
  {
  // render signal:
    int N = (int)ceil(length * sampleRate);
    std::vector<double> xL(N), xR(N);
    sd.noteOn(33, 64);                      // A1, 55Hz, mid velocity
    for(int n = 0; n < N; n++)
      sd.processFrame(&xL[n], &xR[n]);

    // experimental:
    //quantize(xL, xR, 0.025, 0.0, 100.0);
    // ...nah - that doesn't sound good

    //using AT = RAPT::rsArrayTools;
    rsNormalizeJointly(xL, xR);
    int fadeStart = (int)ceil((length-fadeLength)*sampleRate);
    rsFadeOut(&xL[0], fadeStart, N);
    rsFadeOut(&xR[0], fadeStart, N);

    rosic::writeToStereoWaveFile(fileName, &xL[0], &xR[0], N, sampleRate);

    // todo: apply dithering, use the new implementation and save to 24 or 32 bit (without dither)
  };


  // Freq envelope:
  sd.setFreqDecays( 10.0,   80, 240);   // in ms
  sd.setFreqWeights( 6.0,   -1,   1);   // raw factors
  sd.setFreqDepth(800/6.0);             // in Hz

  // Amp Enevlope:
  sd.setAmpDecays(  20.0,  100, 400);
  sd.setAmpWeights(  0.75,  -1,   1); 

  // The freq-env first goes down quickly, then almost stabilizes, and then goes down further more
  // slowly:
  render("SweepBassdrum1.wav", 1.2, 0.2);
  // The freq excursion seems to go far above 800/6 = 13.333. the first cycle is 65 samples long
  // which translates to 48000/65 = 738.46.. Hz. -> figure out what makes the freq-excursion go so
  // high ...ah! the freq starts out at 800 due to freqWeights[0] == 6. Maybe alway choose the 1st
  // weight = 1 such that the freq-depth directly gives the excursion

  // This has a freq-env that goes like down-up-down (might be interesting for toms):
  sd.setFreqDecays( 10.0, 120, 240);
  sd.setFreqWeights( 6.0,  -2,   2); 
  render("SweepBassdrum2.wav", 1.2, 0.2);

  // This has more tonal character, stabilizing at A1 = 55Hz
  sd.setFreqDecays(  5.0,  15,  60);
  sd.setFreqWeights(10.0,  -1,   1);
  sd.setFreqDepth(1000/10.0);
  sd.setAmpDecays(  20.0,  80, 300);
  sd.setFreqFloor(55.0);
  render("SweepBassdrum3_1.wav", 1.2, 0.2);
  // maybe this could use some overtones to emphasize the tonal character

  // Reaching the 55Hz more quickly - it's almost the same as previous but a tad punchier. Maybe 
  // get rid of the previous...but maybe the previous is nice for a lower velocity
  sd.setFreqDecays(4.0,  10,  40);       // different
  sd.setFreqWeights(10.0,  -1,   1);     // same
  sd.setFreqDepth(1300/10.0);            // different
  sd.setFreqFloor(55.0);                 // same
  sd.setAmpDecays(  20.0,  80, 300);     // same
  render("SweepBassdrum3_2.wav", 2.0, 0.2);


  int dummy = 0;

  // ToDo:
  // -Add an envelope for the quantization, use functionrsQuantize (or was it rsQuant)...hmmm...
  //  don't know - maybe quantization is not a good idea. It doesn't sound good. The idea is to 
  //  somehow introduce some noisy element into the attack transient. Simply mixing in decaying
  //  noise may work too, especially in combination with reverb. But that's somehow dissatisfying.
  //  it would be nice to "derive" some noisy signal from the raw bassdrum signal. Maybe we can
  //  drive a chaos generator whose chaoticity is controlled by some aspect of the bassdrum? Maybe
  //  take the derivative, square it, smooth it

  // -allow rendering of the envelopes for inspection
  // -maybe make a class rsSweepDrumPresets that we can use like: preset.setup(sd, 2); where
  //  2 is the index of the preset
  // -Maybe render the raw signal without amp-env. The amp-env can be applied later in the 
  //  sampler.
  // -apply waveshape envelopes - maybe it should go just up and down. at the start, the pitch
  //  is high, so we may not want as many overtones. at the end, it's fading away amplitude-wise,
  //  so the brightness should perhaps follow along. check out zeroDelayfeedbackPhaseMod, maybe
  //  use asin
  // -Split off transient by applying a short amp env at the beginning and subtracting the so 
  //  enveloped signal - save transient and body in seperate files to be recombined with weights
  //  in the sampler - weights can be velocity dependent with stronger dependency for the transient
  //  ...but maybe that can be realized in the sampler itself by using one short decay and one
  //  with an equal attack? then we could remic transient and body in realtime - try in rgc:sfz
  //  to load the same sample twice, for one apply a decay-only (50 ms) and to the other an 
  //  attack-only (also 50 ms) and add the outputs (maybe the release needs to equal the decay). 
  //  I hope that the sum will give back the original signal. If this is the case, we can use this
  //  to apply the amp-envelopes in realtime with dynamic response to velocity. Generally, we want 
  //  the velocity to control the "transientness" and spectral brightness (and also volume). Maybe 
  //  velocity could control a mix between more or less overtones (and/or maybe a filter's cutoff).
  //  Maybe velocity could also control an additional pitch envelope. and 
}

void createBassdrumSamplesOld()  
{
  // rename to createSweepBassdrumSamples or remove and integrate into createMiscSamples

  //createBassdrumPsy1Sample(1.0, true); return;  // for development

  // old - may become obsolete as soon as the new renderSweepBassdrums is finished and re-creates 
  // all samples we create here:
  bool plot = false;
  createBassdrumPsy1Sample(1.0, plot);
  createBassdrumPsy1Sample(1.5, plot);
  createBassdrumPsy1Sample(2.0, plot);
  createBassdrumPsy1Sample(3.0, plot);
  createBassdrumPsy1Sample(4.0, plot);
}

void createNoiseBursts(int sampleRate)
{
  rsNoiseBurst<double> nb;
  nb.setSampleRate(sampleRate);
  nb.setIrwinHallOrder(8);
  nb.setSpectralSlopeStart(0.0);   // start out white

  // Helper function to render output and write to file:
  auto render = [&](const char* fileName, double length, double fadeLength)
  {
    // render signal:
    int N = (int)ceil(length * sampleRate);
    std::vector<double> xL(N), xR(N);
    nb.noteOn(60, 64);                      // A1, 55Hz, mid velocity
    for(int n = 0; n < N; n++)
      nb.processFrame(&xL[n], &xR[n]);

    //using AT = RAPT::rsArrayTools;
    rsNormalizeJointly(xL, xR);
    int fadeStart = (int)ceil((length-fadeLength)*sampleRate);
    rsFadeOut(&xL[0], fadeStart, N);
    rsFadeOut(&xR[0], fadeStart, N);

    rosic::writeToStereoWaveFile(fileName, &xL[0], &xR[0], N, sampleRate);

    // todo: apply dithering, use the new implementation and save to 24 or 32 bit (without dither)
  };

  nb.setAmpAttack(100.0);
  nb.setAmpDecay(2000.0);  // maybe set up decay in terms of T60, maybe call it setReverbTime60
  nb.setSpectralSlopeChange(-5.0);  // in (dB/oct) / sec
  render("NoiseBurst_100_2000_5.wav", 5.0, 0.2);

  nb.setAmpAttack(50.0);
  nb.setAmpDecay(1000.0);
  nb.setSpectralSlopeChange(-10.0);
  render("NoiseBurst_50_1000_10.wav", 2.6, 0.2);

  nb.setAmpAttack(30.0);
  nb.setAmpDecay(800.0);
  nb.setSpectralSlopeChange(-15.0);
  render("NoiseBurst_30_800_15.wav", 1.8, 0.2);

  nb.setAmpAttack(20.0);
  nb.setAmpDecay(500.0);
  nb.setSpectralSlopeChange(-25.0);
  render("NoiseBurst_20_500_25.wav", 1.0, 0.2);

  nb.setAmpAttack(10.0);
  nb.setAmpDecay(300.0);
  nb.setSpectralSlopeChange(-30.0);
  render("NoiseBurst_10_300_30.wav", 0.6, 0.2);

  nb.setAmpAttack(10.0);
  nb.setAmpDecay(300.0);
  nb.setSpectralSlopeChange(-60.0);
  render("NoiseBurst_10_300_60.wav", 0.6, 0.2);



  // settings befor introducing the normalizer for the slope-filter: 
  // 50/100/-20; 50,200,-10; 10,50,-25


  // Observations:
  // -with 50,200,-10, there is a strange noise artifact at the end (after 4.5 secs). With a faster
  //  slope-change, this happens earlier. seems like when the slope hits -45.220833333333331, the
  //  computed DC gain of the slope filter becomes infinite. At -45.707916666666669, it becomes 
  //  NaN. Maybe we need to limit the slope to -45, Or: Maybe use 2 slope filters in series each with 
  //  only half of the slope. Or: maybe make a SlopeFilterChain class that let's the use select the 
  //  number M of filters and each partial filter realizes a slope of slope/M.
  // -When we choose a too high setting for slope-change, the later part of the signal gets 
  //  boosted. I think, it's because the slope filter doesn't normalize to unit magnitude at DC
  //  but at 1 kHz (the pivot frequency). We need to normalize the slope-filter at DC. We can 
  //  counteract the low-freq boost by choosing a shorter decay for the amp-env - but that's 
  //  clunky. We need a variant of SlopeFilter that normalizes at DC - and maybe another one that 
  //  normalizes at Nyquist when we want to introduce a 2nd filter to attenuate the bass over time.
  //  But maybe for that, a time-variant highpass is more appropriate. A 2nd order Butterworth
  //  highpass adjusted around 400 Hz seems to be good
  // -With shorter settings, it could be useful to create snare-like sounds
  //
  // ToDo:
  // -Use a state-variable and/or state-vector based implementation for better time-variant 
  //  behavior
  // -Maybe try other shapes for the slope envelope. It's currently just linear.
  // -Try higher order Irwin-Hall distributions for the input noise. Maybe try (time-varying)
  //  bimodal and trimodal distributions. Maybe it should start out bimodal and over time, the
  //  modes should merge into one, maybe exponentially - i.e. use a bimodal distribution with modes
  //  a,-a for some a and let that a decay to zero.
  // -Try to fix the problem in the slope filter when the slope gets really large. Actually, we 
  //  should plot the magnitude responses anyway. Maybe try to use a filter-bank based approach
  //  instead of the slope filter. Maybe Linkwitz-Riley filters could be good for this

  // Sound Design ideas:
  // -The noise burst could probably be used to create nice pizzicato sounds. Maybe use some sort
  //  of filter-blip as input. Maybe using modal synthesis.
}


void createMiscSamples()
{
  createBassdrumSamplesOld();  // becomes obsolete when renderSweepBassdrums renders all of them
  renderSweepBassdrums(48000); 
  createNoiseBursts(48000);

  rsConvolveFiles("SweepBassdrum1.wav", "NoiseBurst_100_2000_5.wav");
  //rsConvolveFiles("SweepBassdrum1.wav", "NoiseBurst_50_1000_10.wav");
  //rsConvolveFiles("SweepBassdrum1.wav", "NoiseBurst_30_800_15.wav");
  //rsConvolveFiles("SweepBassdrum3_2.wav", "NoiseBurst_30_800_15.wav");
  // The right channel is a bit quiet at the start. maybe try another seed for right. But maybe it
  // has to do the phases of left and right input signal and not with the seed. Maybe try to swap 
  // input channels of one of the signals only. The convolution wet result could use some highpass
  // but this is something we can do in the sampler engine. For the more tonal bassdrum, the reverb
  // does not sound as good. It takes away some of the tonal character (turns everything into a big
  // "tock" sound). Seems like the reverb works better for atonal drums. Maybe try a longer attack
  // for the reverb to leave the transient of the drum more intact.

  return;

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




// This function and its callers may be obsolete now. But the callers still contain some comments 
// that should be preserved:
void applyAllpassChain(const double* x, int N, double* y, double freq, double quality, 
  int numStages, double sampleRate)
{
  // As a convention, we switch to a first order allpass chain when the quality fcator Q is zero:
  bool useBiquad = quality != 0.0;

  // Create and set up the allpass chain filter:
  rosic::AllpassChain apc;
  apc.setSampleRate(sampleRate);
  apc.setFrequency(freq);
  apc.setNumStages(numStages);
  apc.setQ(quality);
  if(useBiquad)
    apc.setMode(apc.SECOND_ORDER_ALLPASS);
  else
    apc.setMode(apc.FIRST_ORDER_ALLPASS);

  // Apply the allpass to x and write the result into y:
  for(int n = 0; n < N; n++)
    y[n] = apc.getSample(x[n]);

  // ToDo: 
  // -Allow for arbitrary number of stages. Currently rosic::AllpassChain has a limit of 24.
}
void applyAllpassChain(std::vector<double>& inOut, double freq, double quality, int numStages,
  double sampleRate)
{
  applyAllpassChain(&inOut[0], (int) inOut.size(), &inOut[0], freq, quality, numStages, sampleRate);
}
void createAllpassBassdrum1()
{
  // This is the first attempt to render a sine-sweepdown based bassdrum sample.

  // User parameters:
  int    sampleRate = 48000;  // Sample rate in Hz
  double length     = 0.5;    // Length in seconds
  double Q          = 1.0;    // Using 0 will switch to 1st order allpass
  double S          = 5;      // Number of stages per allpass chain

  // Render sample:
  int N = ceil(length * sampleRate);  // Number of samples to render
  using Vec = std::vector<double>;
  Vec x(N);
  x[0] = 1;
  applyAllpassChain(x, 16000.0,  Q, S, sampleRate);
  applyAllpassChain(x,  8000.0,  Q, S, sampleRate);
  applyAllpassChain(x,  4000.0,  Q, S, sampleRate);
  applyAllpassChain(x,  2000.0,  Q, S, sampleRate);
  applyAllpassChain(x,  1000.0,  Q, S, sampleRate);
  applyAllpassChain(x,   500.0,  Q, S, sampleRate);
  applyAllpassChain(x,   250.0,  Q, S, sampleRate);
  applyAllpassChain(x,   125.0,  Q, S, sampleRate);
  applyAllpassChain(x,    62.5,  Q, S, sampleRate);
  applyAllpassChain(x,    31.25, Q, S, sampleRate);

  // Normalize and write it to a wavefile:
  RAPT::rsArrayTools::normalize(&x[0], N);
  rosic::writeToMonoWaveFile("AllpassBassdrum1.wav", &x[0], N, sampleRate, 16);
  //rsPlotVector(x);
}
void createAllpassBassdrum2()
{
  // User parameters:
  int    sampleRate = 48000;  // Sample rate in Hz
  double length     = 0.5;    // Length in seconds
  double loQ        = 0.5;
  double hiQ        = 1.0;
  double loF        = 27.5;
  double hiF        = 14080;
  int    numStages  = 1;       // Number of stages per allpass chain
  int    numChains  = 50;      // With 10, we get exact octaves

  // Render sample:
  int N = ceil(length * sampleRate);  // Number of samples to render
  using Vec = std::vector<double>;
  Vec x(N);
  x[0] = 1;
  for(int i = 1; i <= numChains; i++)
  {
    double f = rsLinToExp(double(i), 1.0, double(numChains), loF, hiF);
    double q = rsLinToExp(double(i), 1.0, double(numChains), loQ, hiQ);
    applyAllpassChain(x, f, q, numStages, sampleRate);
  }

  // Normalize and write it to a wavefile:
  RAPT::rsArrayTools::normalize(&x[0], N);
  rosic::writeToMonoWaveFile("AllpassBassdrum2.wav", &x[0], N, sampleRate, 16);

  // Observations:
  // -When loQ=4 and hiQ=1, then the sound becomes kinda warbly. It's also warbly when having it the 
  //  other way around
  // -Using more stages per chain (and reducing the number of chains to keep the product constant),
  //  the amp-envelope seems to become more steppy. There seem to be plateaus. With just 1 stage 
  //  per chain, the amp-env does not feature such plateaus
  //
  // Conclusions:
  // -We may scrap the numStages parameter and just always use 1 stage
  //
  // ToDo:
  // -Maybe let the frequencies and Q-values pass through a nonlinear function - maybe the linfrac
  //  warping map should be applied to the normalized values
}
void createAllpassBassdrum3()
{
  // User parameters:
  int    sampleRate = 48000;  // Sample rate in Hz
  double length     = 0.5;    // Length in seconds
  double loQ        = 1.0;
  double hiQ        = 1.0;
  double shQ        = 0.0;    // Shape parameter for Q
  double loF        = 27.5;
  double hiF        = 14080;
  double shF        = 0.0;    // Shape parameter for frequency
  int    numFilters = 50;     // Number of allpass filters

  // Render sample:
  int N = ceil(length * sampleRate);  // Number of samples to render
  using Vec = std::vector<double>;
  Vec x(N);
  x[0] = 1;
  auto shape = [](double x, double s) { return RAPT::rsRationalMap_01(x, s); };
  for(int i = 0; i < numFilters; i++)
  {
    double p = double(i) / double(numFilters-1);               // Goes from 0 to 1
    double f = rsLinToExp(shape(p, shF), 0.0, 1.0, loF, hiF);
    double q = rsLinToExp(shape(p, shQ), 0.0, 1.0, loQ, hiQ);
    applyAllpassChain(x, f, q, 1, sampleRate);
  }

  // Normalize and write it to a wavefile:
  RAPT::rsArrayTools::normalize(&x[0], N);
  rosic::writeToMonoWaveFile("AllpassBassdrum3.wav", &x[0], N, sampleRate, 16);

  // Observations:
  // -Positive values of shF make the sweep initially faster, I think. It sounds more snappy.
  // -For Q=0, we get a DC output. That's wrong!
}

//=================================================================================================

/** This is a class a for post-processing the rendered samples. It includes functions for common
post-filtering operations for spectral shaping (for example white-to-pink) or clean-up (for example 
highpass below desired fundamental) or boosting/cutting certain frequencies with peaking filters to
impose a tonality. It also has functionality to shorten samples, apply fade-out, etc. ...TBC... */

class rsSamplePostProcessor
{

public:

  //-----------------------------------------------------------------------------------------------
  // \Setup

  void setSampleRate(double newSampleRate) { sampleRate = newSampleRate; }


  //-----------------------------------------------------------------------------------------------
  // \Processing

  void applyOnePoleLowpass(double* x, int N, double cutoff);
  // Maybe add optional parameters: bool bidirectional = false, int numPasses = 1


  void applyOnePoleHighpass(double* x, int N, double cutoff);


  // ToDo: applyTilt, applyOnePoleHighpass, applyButterworthHighpass, applyEllipticHighpass, 
  // shortenTail(cutThresholdDb, fadeTime) etc.


  //-----------------------------------------------------------------------------------------------
  // \Convenience

  void applyOnePoleLowpass(std::vector<double>& x, double cutoff)
  { applyOnePoleLowpass(&x[0], (int) x.size(), cutoff); }

  void applyOnePoleHighpass(std::vector<double>& x, double cutoff)
  { applyOnePoleHighpass(&x[0], (int) x.size(), cutoff); }

  /** Shortens the tail of the signal x. Cuts off everything from the end that falls below a given 
  threshold in dB. The threshold is interpreted as being relative to the maximum sample value. 
  After shortening the length, it then applies a smooth fade-out envelope to the new end of given 
  length in seconds. A good value for the threshold is -60 dB and for the fade out time 0.02 
  seconds, i.e. 20 milliseconds. */
  void shortenTail(std::vector<double>& x, double thresholdDb, double fadeOutTime, 
    double releaseTime);


protected:

  double sampleRate = 1.0;

};
// Maybe move to somewhere else - maybe rosic/rendering
// See also rsBiDirectionalFilter in MiscUnfinished.h in rapt/Unfininished/Sampling. It has similar
// functionality - maybe merge.

void rsSamplePostProcessor::applyOnePoleLowpass(double* x, int N, double cutoff)
{
  RAPT::rsOnePoleFilter<double, double> flt;
  flt.setSampleRate(sampleRate);
  flt.setMode(flt.LOWPASS_IIT);
  flt.setCutoff(cutoff);
  for(int n = 0; n < N; n++)
    x[n] = flt.getSample(x[n]);
}

void rsSamplePostProcessor::applyOnePoleHighpass(double* x, int N, double cutoff)
{
  RAPT::rsOnePoleFilter<double, double> flt;
  flt.setSampleRate(sampleRate);
  flt.setMode(flt.HIGHPASS_MZT);
  flt.setCutoff(cutoff);
  for(int n = 0; n < N; n++)
    x[n] = flt.getSample(x[n]);
}

void rsSamplePostProcessor::shortenTail(std::vector<double>& x, double thresholdDb, 
  double fadeOutTime, double releaseTime)
{
  RAPT::rsEnvelopeFollower<double, double> ef;
  ef.setSampleRate(sampleRate);
  ef.setAttackTime(0.0);
  //ef.setReleaseTime(1000 * 0.25 * 1.0/lowFreq);
  ef.setReleaseTime(1000 * releaseTime);
  int N = (int) x.size();
  std::vector<double> env(N);
  for(int n = 0; n < N; n++)
    env[n] = ef.getSample(x[n]);
  //rsPlotVectors(x, env);
  // I tried to use the more advanced rsEnvelopeFollower2 but this produced total garbage results
  // for this sort of signal. The simpler works much better.

  // Find last sample that exceeds the threshold:
  double envMax   = RAPT::rsArrayTools::maxValue(&env[0], N);
  double thresh   = RAPT::rsDbToAmp(thresholdDb);
  thresh *= envMax;  // threshold should be relative
  int nCut = N-1;
  while(nCut > 0)
  {
    if(env[nCut] >= thresh)
      break;
    nCut--;
  }

  // Shorten the signal:
  N = nCut+1;
  x.resize(N);

  // Apply a smooth fade out:
  int fadeSamples = sampleRate * fadeOutTime;
  rsFadeOut(&x[0], N-fadeSamples-1, N-1);
  //rsPlotVectors(x, env);
}


// Refactor this into:
// -Generation of exciter signal x
//  -Maybe make a class rsWhiteExciter with various options:
//   -impulse, impulse-train, allpass-delay-chain, chirpUp, chirpDown, noise, allpass-FDN, 
//    noise-burst
// -Application of the allpass zapper
// -Post-processing
//  -Browning lowpass
//  -Cleanup highpass
//  -Shortening
//  -Normalization
// -Filename generation
// -File writing
void createBrownZap(int numStages, double lowFreq = 15, double highFreq = 8000, double freqShape = 0.0, 
  double lowQ = 1.0, double highQ = 1.0, double qShape = 0.0, double maxLength = 1.0, int sampleRate = 48000)
{
  // The result is like in the createAllpassBassdrumN() functions above but here, we use the class
  // rsWhiteZapper which encapsulates the allpass based algorithm which the other functions 
  // implement manually.

  // Create and set up the zapper object:
  rosic::rsWhiteZapper wz;
  wz.setSampleRate(sampleRate);
  wz.setNumStages(numStages);
  wz.setLowFreq(lowFreq);
  wz.setHighFreq(highFreq);
  wz.setFreqShape(freqShape);
  wz.setLowQ(lowQ);
  wz.setHighQ(highQ);
  wz.setQShape(qShape);

  // Render sample:
  int N = ceil(maxLength * sampleRate);  // Number of samples to render
  using Vec = std::vector<double>;
  Vec x(N);
  x[0] = 1;

  // Test - with noise-burst:
  //x = createNoise(N, -1.0, +1.0, 0);
  //Vec e = attackDecayEnvelope(N, 50, 150);  // Decay = 100-200 seems nice
  //x = x*e;
  //rsPlotVector(x);
  // Using a noise-burst makes it sound more acoustic and natural, less electronic. But maybe this
  // doesn't count as raw-material anymore because it can be recreated from the impulse responses 
  // using convolution. Actually, the user could just play the kick sample through a convolution
  // reverb that uses a noise-burst as impulse-response

  for(int n = 0; n < N; n++)
    x[n] = wz.getSample(x[n]);
  //rsPlotVector(x);


  // Optionally plot a phase-spectrum of the allpass impulse response:
  bool plotPhase = true;  // Maybe make this a function parameter
  if(plotPhase)
  {
    SpectrumPlotter<double> sp;
    sp.setFftSize(65536);
    sp.setLogFreqAxis(true);
    sp.setSampleRate(sampleRate);
    sp.setFreqAxisUnit(SpectrumPlotter<double>::FreqAxisUnits::hertz);
    //sp.plotDecibelSpectra(N, &x[0]);  // Should be flat. Yep - it is.
    sp.plotPhaseSpectra(N, &x[0]);
    // ToDo: Maybe alternatively plot phase-delay or group-delay
  }


  // Post-process the white zap. First we turn the spectrum from white to brown by applying a first
  // order lowpass tuned somewhere below the lowest allpass tuning freq. The resulting -6 dB/oct 
  // magnitude response is ideal in the sense that the amplitude stays constant during the sweep. 
  // Then we do a bit of cleanup by applying a highpass. This removes some bumpiness that is 
  // introduced by the lowpass (it looks like we get some sort of undulating DC towards the end in
  // the lowpass - we want to remove that). Finally, the sample is shortened to remove the silence
  //  in the tail:
  rsSamplePostProcessor pp;
  pp.setSampleRate(sampleRate);
  pp.applyOnePoleLowpass( x, 0.5*lowFreq);
  pp.applyOnePoleHighpass(x, 1.0*lowFreq);
  pp.applyOnePoleHighpass(x, 1.0*lowFreq);
  pp.applyOnePoleHighpass(x, 1.0*lowFreq);
  pp.shortenTail(x, -60.0, 0.02, 0.25/lowFreq); 
  N = x.size();
  RAPT::rsArrayTools::normalize(&x[0], N);  // implement and use sp.normalize(x);
  // Maybe a 3rd order Butterworth highpass would be better than applying a 1st order highpass 3
  // times? Try it! Maybe also try a somwhat lower cutoff for the highpass like 0.75*lowFreq






  // Create filename from the parameters (maybe factor out):
  std::string name = "ZappyKick"; // Nah - not all possible settings lead to bassdrums
  //std::string name = "AllpassZap";
  name += "_NS=" + std::to_string(numStages);
  if(lowFreq  != 15)   name += "_FL=" + rosic::rsToString(lowFreq);
  if(highFreq != 8000) name += "_FH=" + rosic::rsToString(highFreq);
  name += "_FS=" + rosic::rsToString(freqShape);
  if(lowQ   != 1.0)  name += "_QL=" + rosic::rsToString(lowQ);
  if(highQ  != 1.0)  name += "_QH=" + rosic::rsToString(highQ);
  if(qShape != 0.0)  name += "_QS=" + rosic::rsToString(qShape);
  // ToDo: add mode
  name += ".wav";
  // Maybe let the caller pass a basic name - don't hardcode the "ZappyKick" name here. It may
  // default to AllpassZap. The function may be eventually used to render other allpass based 
  // sounds such as toms and snares.


  // Write it to a wavefile:
  rosic::writeToMonoWaveFile(name, &x[0], N, sampleRate, 16);





  // ToDo: factor this out:
  // Create a spectrogram:
  //rsSpectrogramProcessor<double> specProc;
  //specProc.setBlockAndTrafoSize(256, 2048);
  //specProc.setHopSize(64);
  //specProc.setAnalysisWindowType(...);


  // Plot a spectrogam:
  //SpectrogramPlotter<double> sp;
  //GNUPlotter plt;
  // I think, the way the SpectrumPlotter class works is to call appropriate setup functions on an
  // existing GNUPlotter object.


  // Observations:
  // -When applying the zapper twice, the result is similar to using a single zapper with twice the
  //  number of stages, I think. So, we don't really get anything interesting new by chaining 
  //  multiple zappers.
  // -The Q factor does not seem to affect the sweep-speed but the waveshape. For low Q, it's 
  //  sinusoidal. For high Q, it looks like harmonics appear. Higher Q-values should be used only
  //  at high frequencies, i.e. highQ can be somewhat higher, but lowQ should always be low-ish. I
  //  think, with lower Q, we get a cleaner time-separation of the different frequencies. With 
  //  higher Q, they get mixed up more - or something.
  // -Positive values of freqShape (freqShape) like +0.8 have the effect of cramming the signal more 
  //  together. But the attack also gets longer - it gets crammed around some center that is not at
  //  the start of the sample. Negative values like -0.8 seem to lengthen the thump/body portion of
  //  the zap and shorten the transient.
  // -Lowering highFreq makes the transient more clicky. Raising highFreq makes the transient more 
  //  sweepy/zappy.
  // -With a slope of -6, we seem to get a constant amplitude. Maybe try using a leaky integrator,
  //  i.e. lowpass with 6dB/oct tuned to somewhere below lowFreq
  // -When increasing lowFreq, it does not sweep down as much and the whole sweep happens faster. I 
  //  think cranking up lowFreq by an octave makes the zap happen in half of the time
  // -Playing the generated sample back at different speeds does not change the sonic impression 
  //  very much
  // -The highest frequency used affects how much the transient is smeared over time. When using 
  //  1kHz as highes allpass freq, we get a strong impulse at the start. When usiing 16 kHz as 
  //  highest freq, the initial impulse is more smeared out over time. I think, The higher we go,
  //  the more smeared out the higher frequencies will be. That means, we can use that parameter to
  //  determine the initial sweep speed? Going higher makes the sweep initially slower or 
  //  something? Or to put it the other was: when the highest freq is rather low, then the high 
  //  frequencies are more crammed into the initial section.
  //  -I think, the sweet spot is around 8000. Going up to 16 doesn't seem to make much of an 
  //   audible difference - but maybe that's just my ears
  // -Previously, I was using rosic::AllpasChain which implements a chain of up to 24 allpass 
  //  filters that all use the same set of coefficients. With this setup, using more stages per 
  //  chain (and reducing the number of chains to keep the product constant), the amp-envelope 
  //  seems to become more steppy. There seem to be plateaus. With just 1 stage per chain, the 
  //  amp-env did not feature such plateaus. We could simulate this old setup by letting
  //  rsWhiteZapper use groups of allpasses that have the same settings. For example, for a total
  //  number of 50 allpasses, we could use 10 groups of 5 or 5 groups of 10 or 2 groups of 25. I 
  //  don't think, that's desirable, though.
  // -With lower Q, its seems to sweep further down at the end.
  // -Using a sawtooth burst as input is not so interesting (maybe try an impulse-train burst?)
  // -The phase response always goes down from 0° to -(numStages*360°) in a sigmoid shape.
  // -The almost linear portion of the sigmoid is confined between lowFreq and highFreq.
  // -The midpoint of the sigmoid occurs at the geometric mean of lowfreq and highFreq (when 
  //  freqShape = 0) - tested with loFreq = 500, hiFreq = 2000 and loFreq = 250, hiFreq = 4000. in
  //  both cases, the phase response's symmetry center goes through f = 1000.
  // -Higher values for Q make the sigmoid sharper (i.e. more resembling a hard-clipper) at the 
  //  expense of making the linear section wiggly (try Q = 16 with numStages = 10). The number of
  //  wiggles seem to equal numStages. At around Q = 4, the wiggles become barely visible in the 
  //  plot, with Q = 2 completely invisible, with Q = 1 we may actually already be in the 
  //  oversmoothed range - although, it's still fine. (ToDo: check for other values of numStages).
  //  Try createBrownZap(50, 100, 1000, 0.0, Q, Q, 0.0); and tweak Q
  // -By the way, the allpass chain *closely* approximates a linear phase response.
  // -Shape parameters other than 0 skew the sigmoid and make it asymmetric - maybe a bit like
  //  the Gompertz function.
  // -One octave sweeps from around 200 to 100 Hz make for nice electric toms:
  //  createBrownZap(45, 100, 200, -3.0); Pulling the lowfreq and highFreq even closer together
  //  does not change the sound very much. Even using 149 and 151, it sounds still rather similar.
  // -It actually seems to work to have lowFreq == highFreq
  // -for createBrownZap(45, 25, 200, -3.0); the phase plot looks wrong. Maybe try increasing 
  //  maxLength - maybe it's because of truncation
  //
  // Conclusions:
  // -Overall length is proportional to numStages and inversely proportional to lowFreq
  // -It makes sense to let it sweep down to around 15 Hz. If a sweep is desired that sweeps faster
  //  and only down to 30 Hz, we can just play the sample back at twice the speed. This will give a
  //  similar sound.
  // -Giving the user control over how Q changes with frequency seems overkill. Maybe a single Q
  //  parameter is good enough. ...but maybe experiment a bit more... Or maybe I have just not set 
  //  used it in the right range - I tried 1..16 - but maybe we should stay the smooth-oversmooth 
  //  range like 0.125..4 or something. Maybe it makes sense to have a sharper transition at one 
  //  end of the sigmoid.
  // -I think, the most important sound-shaping feature to be added is to give more flexibility to
  //  the curve that distributes the allpass tuning frequencies. Using a fixed Q fo all allpasses
  //  that migght even be hardcoded seems good enough. But maybe when we have more flexibility for
  //  the freqs, more flexibility for the Qs may become more desirable? We'll see. However, the
  //  current implementation of the Q-shape should be kept - but it doesn't need to be a user 
  //  parameter on a GUI, when I make a module for ToolChain from it.
  //
  // ToDo:
  // -Maybe instead of writing the sample to disk, return it a std::vector. Or maybe factor out a 
  //  getWhiteZapBassdrum function
  // -We currently need numStages >= 2. For 0 and 1, it produces garbage (inf, nan). Fix this!
  // -Test it with extreme settings like having the highest freq at fs/2. What if the lowest freq 
  //  is 0? Does the allpass design admit this? However, the exponential scaling will not admit it
  //  anyway.
  // -Other interesting post-processing effects could be a peak-eq tuned to the desired bassdrum
  //  fundamental frequency.
  // -Let the function take the parameters as function arguments
  // -Try a first order allpass chain
  // -Plot a spectrogram
  // -Plot a group-delay response
  // -Plot a ringing response
  // -Use an EQ to boost the bass frequencies
  // -Try to make a phaser from that
  // -Check out, if we can get a nice stereor signal by using noise-bursts with different seeds for
  //  left and right (or for mid and side)
  // -In a synth based on that, allow to mix noise-burst and impulse inputs - and maybe other types
  //  of inputs as well
  // -Try applying some of the allpasses in forward mode and some in backward mode. Maybe use lower
  //  Q for backward mode. When using the exact same allpass-chain for forward and backward pass,
  //  the phase-responses shuld cancel and the original impulse should be reconstructed
  // -Figure out the shape of the pitch env by using the zero-crossing based pitch detector. Or:
  //  use class rsSingleSineModeler. This may give better time resolution. or maybe try both 
  //  approaches
  //
  // Ideas:
  // -Maybe the allpass action cannot only be imagined as delaying frequencies but also as 
  //  lengthening the event at which that frequency-band occurs by simulatneously lowering its
  //  amplitude to keep the energy constant. It would be nice to be able to set up an allpass in 
  //  terms of its "lengthening" response. Maybe that's what group delay is, after all?
}


void createAllpassDrums()
{
  // The idea is to use the impulse responses of allpass filters to create percussive sounds. These 
  // sounds will occupy the full spectrum, i.e. they will be white. That's a good starting point 
  // for a drum sound. Later, it may be further shaped by using lopwass/highpass/bandpass/peak/etc.
  // filters.

  // Obsolete:
  //createAllpassBassdrum1();
  //createAllpassBassdrum2();
  //createAllpassBassdrum3();

  // Experimental:
  //createBrownZap(30,  100,  500, 0.0);   // Tom?
  //createBrownZap(50,  100, 2000, 0.0);   // Laser Zap
  //createBrownZap(50,  250, 4000, 0.0);   // 4 octaves
  //createBrownZap(50,  250, 500, 0.0);      // 1 octave

  // For plotting tests:
  createBrownZap(50, 10, 1000, 0.0);


  // Nice soft bassdrum or low tom:
  createBrownZap(45, 25, 200, -3.0);  // phase-plot looks wrong - I think this is a bug in the plotter
  createBrownZap(45, 25, 200,  0.0);
  createBrownZap(45, 25, 200, +3.0);

  // One octave sweeps from around 200 to 100 Hz make for nice electrnic toms:
  createBrownZap(45, 100, 200, -3.0);
  createBrownZap(45, 120, 180, -3.0);
  createBrownZap(45, 140, 160, -3.0);
  createBrownZap(45, 149, 151, -3.0);
  // They sound all quite similar.





  // Create the ZappyKickXXX.wav bassdrums:
  int numStagesLo  = 20;
  int numStagesHi  = 80;
  int numStagesInc = 10;
  double shapeLo   = -5;
  double shapeHi   = +5;
  for(int numStages = numStagesLo; numStages <= numStagesHi; numStages += numStagesInc)
  {
    for(int k = shapeLo; k <= shapeHi; k++)
    {
      double freqShape = k;
      createBrownZap(numStages, 15, 8000, freqShape);
    }
  }
  // My favorites: 20/-2, 20/0, 30/-3, 40/-3, 40/1..2, 50/-3, 50/1..2, 60/-4, 60/2, 70/-4, 70/2..3,
  // 80/-4, 80/3 where 2..3 means the perfect setting might be in between 2 and 3 etc..

  // ToDo:
  // -Allow for more flexible shaping like in the linear fractional interpolation scheme. We want 
  //  to determine the slope at 0 as s0 = log2(param1), s1 = pow(s0, param2-1), shape = param3.
  //  param2 is some sort of SigmoidVsSpikey parameter and param3 is a shift up/down parameter, 
  //  I think
  // -Make a function createNoiseBurstDrums
  // -Instead of using noise-bursts as excitatition, try using impulse responses of allpass-delay
  //  series. Maybe these clicks can actually be used for another sample pack.
  // -Generally, to synthesize drums, we may use layers of such allpass-based "zap" sounds and 
  //  noise-burst based sounds. The noise bursts may also pass through a time-varying filter and
  //  waveshaper(s). The overall mix, too.
  // -It would be nice to have two user parameters: length, snap which would adjust the freqShape 
  //  and numStages
  // -Maybe provide 16/44.1 versions for free and sell 24/96 versions through one of these stores:
  //  https://toneisland.com/sites-to-sell-sample-packs/
  //  Loopmasters seems interesting...or maybe use Patreon
  //  Maybe let the big, commercial pack also have more samples, i.e. sample the parameters more 
  //  densely. From the higher quality alone, the commercial pack would be larger by a factor of
  //  3.2653 = (96/44.1) * (24/16). If we sample the parameters with a density that higher by a 
  //  factor of k, we'd get an additional factor of k^numParameters. But maybe sample only the 
  //  one or two most important parameters more densely.
  // -For normalized sample packs, we may need to store the applied gain factor somewhere. Maybe in
  //  metadata? When sfz files are generated, the gain can be stored there.
}

void createSamplerWaveforms()
{
  // ToDo: factor out the mip-map creation and move it to RenderScriptTools.h/cpp

  // Create mip-mapped multisamples for single-cycle waveforms that can be used in the sampler 
  // engine. Most tables are of length 2048 except the bottom one which is of length 4096. We want
  // a cycle length that is a power of two and we want the samples to have a rootKey that is 
  // exactly a midi note. All midi notes except for the As have irrational frequencies. So let's 
  // choose an A as rootKey. A4 at 440Hz is midi-key 69. When we use key=21 (27.5 Hz) for the table
  // of length 2048, we need to choose a sample-rate of 56320 to make it all work out. The formula 
  // for the frequency is: sampleRate/cycleLength, so we get 56320/2048 = 27.5

  // ToDo:
  // -Create prototype wave of length 8192
  // -Create various filtered versions moving average filters...the goal is to preserve the 
  //  time domain waveform as closely as possible ...maybe use a 
  //  filter -> decimate -> interpolate procedure for the higher tables
  //  ...maybe to avoid boundary artifacts from the MA filters, we should use 3 cycles for the
  //  filtering and then extract the middle one

  using Vec = std::vector<double>;
  using SWR = rosic::StandardWaveformRenderer;
  using namespace RAPT;
  using namespace rosic;

  std::string waveName = "Saw";  // name of the waveform as it appears in the .wav files

  // Render prototype sawtooth wave of length 8192 and obtain the first two decimated versions of
  // it by just suing averages of successive samples:
  Vec w8192(8192);
  SWR::renderSawWaveform(&w8192[0], 8192); //rsPlotVector(w8192);
  rsScale(w8192, 0.8125); //
  // Some nicely representable approximation to the scaling by the Wilbraham-Gibbs constant to 
  // avoid clipping due to Gibb's overshoot 
  // https://en.wikipedia.org/wiki/Gibbs_phenomenon
  // https://mathworld.wolfram.com/Wilbraham-GibbsConstant.html


  Vec w4096 = rsDecimateViaMean(w8192, 2); //rsPlotVector(w4096);
  Vec w2048 = rsDecimateViaMean(w4096, 2); //rsPlotVector(w2048);

  std::string fileName;
  fileName = waveName + "_K9.wav";
  rosic::writeToMonoWaveFile(fileName.c_str(), &w4096[0], 4096, 56320, 16);

  //fileName = waveName + "_K21.wav";
  //rosic::writeToMonoWaveFile(fileName.c_str(), &w2048[0], 2048, 56320, 16);
  // seems like we don't need the 8192 version..


  // From the waveform of length 2048, we create the mip-map using the FFT/iFFT technique:
  MipMappedWaveTableStereo mipMap;
  double* pWave[2];
  pWave[0] = pWave[1] = &w2048[0];
  mipMap.setWaveform(pWave, 2048);
  Vec tmp(mipMap.getTableLength());  // temp buffer for the succesive mip-map levels
  int key = 21;
  for(int i = 0; i < mipMap.getNumLevels(); i++)
  {
    mipMap.copyDataTo(&tmp[0], 0, i);
    fileName = waveName + "_K" + std::to_string(key) +  ".wav";
    rosic::writeToMonoWaveFile(fileName.c_str(), &tmp[0], 2048, 56320, 16);
    key += 12;
    //rsPlotVector(tmp);
  }

  // Observations:
  // -In the 0th mip-map, it seems like the Nyquist freq is missing?
  // -Allow the user to specify a tapering function such that we do not necessarily have to use
  //  hard brickwall filters with all of their strong ringing
  //
  // ToDo:
  // -Figure out (and maybe fix) the missing Nyquist freq in the 0th mip-map level.
  //  -Commenting out the "Truncate the spectrum for the next iteration" part doesn't fix it.
  //  -Maybe it has to do with the weird encoding of DC and Nyquist gain in the 0th spectral bin?
  // -Looks like we could have one more level...but it will consume more memory
  // -[Done?] Write the different rendered levels to wavefiles with a meaningfully formatted way, 
  //  such as Saw_K21 for the w2048 wave.

  int dummy = 0;
}