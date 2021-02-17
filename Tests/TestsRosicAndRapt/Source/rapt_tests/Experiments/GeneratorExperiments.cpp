using namespace RAPT;

void fractalize(double* x, double* y, int N, int a = 2, int b = 0)
{
  int inc = 1;
  while(inc <= N/2)  // is this the right condition?
  {
    double amp = 1.0 / double(inc);  // amplitude - maybe raise this to a power p

    for(int n = 0; n < N; n++)
      y[n] += amp * x[(inc*n)%N];


    inc = a*inc + b;
  }
}
// todo: make anothenr version that doesn't multiply by 2 in each iteration but adds 1...or maybe
// in general: have a factor a and offset b and the the next freq is fNew = a*fOld + b. That would
// cover both cases and many more
// -apply the function recursively to the fractalized wave

void waveformFractalization()
{
  static const int N = 2048;  // period length of prototype wave in samples
  int a = 3;  // must be at least 1
  int b = 0;

  double x[N], y1[N], y2[N];

  // render prototype waveform:
  double w = 2*PI/N;
  for(int n = 0; n < N; n++)
  {
    //x[n] = sin(w*n);   // todo: allow other prototype waveforms later
    //x[n] = rsSawWave(w*n);
    x[n] = rsSqrWave(w*n);  // fractalized with a=2,b=0 gives saw wave
    //x[n] = rsTriWave(w*n);
  }

  // fractalize it:
  rsArrayTools::fillWithZeros(y1, N); // maybe rename to clear - or make alias
  fractalize(x, y1, N, a, b);

  // fractalize the fractalized wave again:
  rsArrayTools::fillWithZeros(y2, N);
  fractalize(y1, y2, N, a, b);


  // plot:
  //rsPlotArray(x, N);
  //rsPlotArray(y, N);
  rsPlotArrays(N, y1);
  //rsPlotArrays(N, x, y1, y2);
  // todo: make the plot nicer - introduce a time-axis from 0 to 1, normalize the waveform, let
  // y-axis range from -1.1 to +1.1

  // Observations:
  // -a=2, b=0: good for saw, square becomes saw
  // -a=2, b=1: nice irregularities for saw
  // -a=1, b=1: sine becomes saw, square also interesting
  // -a=3, b=0: good for square
  // -fractalizing the fractalized wave again seems to emphasize the higher frequencies but doesn't
  //  seem to bring something qualitatively new (or does it? ...more experiments needed - maybe 
  //  compare results with using an exponent on the amplitudes to make them fall off more slowly)

  // todo: try, hwo it looks like on an xy-oscilloscope when we introduce a phase-shift between
  // left and right (maybe that phase shift can be time-varying)...this can be done in 
  // straighliner
}

void plotHistogram(const std::vector<double>& data, int numBins, double min, double max)
{
  // under construction - it's kinda ugly but at least, it works

  GNUPlotter plt;
  plt.addDataArrays((int)data.size(), &data[0]);
  plt.addCommand("reset");
  plt.addCommand("n=" + to_string(numBins));
  plt.addCommand("max=" + to_string(max));
  plt.addCommand("min=" + to_string(min));
  plt.addCommand("width=(max-min)/n");  // interval width
  plt.addCommand("hist(x,width)=width*floor(x/width)+width/2.0"); // function used to map a value to the intervals
  plt.addCommand("set xrange [min:max]");
  plt.addCommand("set yrange [0:]");
  plt.addCommand("set offset graph 0.05,0.05,0.05,0.0");
  plt.addCommand("set xtics min,(max-min)/5,max");
  plt.addCommand("set boxwidth width*0.9");
  plt.addCommand("set style fill solid 0.5");
  plt.addCommand("set tics out nomirror");
  //plt.addCommand("set xlabel \"x\"");
  //plt.addCommand("set ylabel \"Frequency\"");
  plt.addCommand("plot \"" + plt.getDataPath() 
    + "\" u (hist($1,width)):(1.0) smooth freq w boxes lc rgb\"black\" notitle");
  plt.invokeGNUPlot();

  // adapted from:
  // https://stackoverflow.com/questions/2471884/histogram-using-gnuplot
  // http://gnuplot-surprising.blogspot.com/2011/09/statistic-analysis-and-histogram.html
}
// move to plotting

// Computes the period of an Irwing-Hall PRNG of given order when the period of the underlying raw
// PRNG is prngPeriod
rsUint32 getIrwinHallPeriod(rsUint32 prngPeriod, rsUint32 order)
{
  rsUint32 state = 0;
  rsUint32 i = 0;
  // we move forward in steps of order, wrapping around at the raw prng period - when we encounter 
  // a zero again, we have hit the period of the Irwing-Hall generator
  while(true)
  {
    state = (state + order) % prngPeriod;
    i++;
    if(state == 0)
      break;
  }
  return i;
}

bool testIrwinHallPeriod(rsUint32 period)
{
  for(rsUint32 i = 1; i <= period; i++)
  {
    rsUint32 p1 = getIrwinHallPeriod(period, i);  // period of Irwing-Hall generator
    rsUint32 p2 = period / rsGcd(period, i);      // conjectured formula - seems to work
    if(p1 != p2)
      return false;
  }
  // Result with period = 20:
  // i: 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20
  // p: 20 10 20  5  4 10 20  5 20  2 20  5 20 10  4  5 20 10 20  1  the ihp
  // the full period is obtained when i and 20 have no common factors
  // maybe return a bool, if the test works out

  return true;
  // return true, if conjecture holds for given period
}

// tests all periods up to maxPeriod
bool testIrwinHallPeriods(rsUint32 maxPeriod)
{
  for(rsUint32 p = 1; p <= maxPeriod; p++)
  {
    if(!testIrwinHallPeriod(p))
      return false;
  }
  return true;
}


double irwinHall(double z, int n)
{
  int zi = (int) floor(z);
  double sum = 0;
  for(int k = 0; k <= zi; k++)
  {
    //sum += pow(-1, k) * rsBinomialCoefficient(n, k) * pow(z-k, n); // mass
    sum += pow(-1, k) * rsBinomialCoefficient(n, k) * n * pow(z-k, n-1); // density
  }
  return sum / rsFactorial(n);
}
// for z > n, the result is wrong
// https://www.youtube.com/watch?v=-2PA7SbWoJ0&t=17m50s (german)

void plotIrwinHall(int order)
{
  int N = 1001;
  double xMin = 0;
  double xMax = order;

  using Vec = std::vector<double>;
  Vec x = rsRangeLinear(xMin, xMax, N);
  Vec y(N);
  for(int n = 0; n < N; n++)
    y[n] = irwinHall(x[n], order);

  rsPlotVectorsXY(x, y);
}

void noise()
{
  // We generate noise with different distributions and plot the histograms...
  int numSamples = 100000;
  int numBins    = 50;
  int order      = 3;  // 5 looks good and is not a power of 2 -> sweet spot
  // i think, powers of two are not good for the order because the reduce the period length
  // todo: compute formulas for the period length as function of the order - the underlying PRNG
  // has a period of 2^32, when the order K is a power of two, the resulting period will be 
  // 2^32 / K (i think -> verify), but when K is not a power of two - what then? what if K==7,
  // for example. i think, maybe the period ength must be divided by the factor 2^k where k is 
  // the power of 2 which goes into K - for example, if K = 12, we must divide by 4 because 4
  // divides 12? maybe numbers which have no factor of two will leave the period as is, so
  // all odd numbers should be good? -> research needed (maybe make numeric experiments with a 
  // generator with a shorter period like 2^16)
  // But: if we run 4 or 8 prngs in parallel, each with its own state (each starting with a 
  // different seed), the period would still be 2^32. This may also be nicely parallelizable


  // some preliminray experiments:
  bool periodFormulaWorks = testIrwinHallPeriods(100);  // check conjectured formula
  plotIrwinHall(order);  // plot theoretical distribution

  rsNoiseGenerator2<double> prng;
  prng.setOrder(order);

  using Vec = std::vector<double>;

  int N = numSamples;
  Vec x(N);
  for(int n = 0; n < N; n++)
    x[n] = prng.getSample();

  //rsPlotVector(x);
  plotHistogram(x, numBins, -1.0, +1.0);

  std::string name = "NoiseIrwinHall" + to_string(order) + ".wav";
  rosic::writeToMonoWaveFile(name, &x[0], N, 44100); 

  // todo: make a function 
  // testNoiseGen(int numSamples, int order, bool plotHistogram, bool writeWaveFile)
}

void noiseTriModal()
{
  int numSamples = 100000;
  int numBins    = 50;
  int order      = 7;

  rsNoiseGeneratorTriModal<double> ng;
  ng.setOrder(order);
  ng.selectorLowpass.setCutoff(500);
  ng.selectorLowpass.setSampleRate(44100);
  ng.selectorLowpass.setMode(rsOnePoleFilter<double, double>::LOWPASS_IIT);

  using Vec = std::vector<double>;

  int N = numSamples;
  Vec x(N);
  for(int n = 0; n < N; n++)
    x[n] = ng.getSample();

  plotHistogram(x, numBins, -1.0, +1.0);

  std::string name = "NoiseTriModal" + to_string(order) + ".wav";
  rosic::writeToMonoWaveFile(name, &x[0], N, 44100); 
  // This noise has some crackly behavior

  // it seems like when we decrease the cutoff freq, we should also decrease the width of the
  // center section between thresh1 and thresh2 - otherwise the outer generators will be selected
  // to rarely
}

void noiseWaveShaped()
{
  // We create uniform white noise, apply a first order filter, predict the resulting amplitude
  // distribution and use the corresponding cumulative distribution to shape the amplitude 
  // distribution back to uniform. The eventual goal is to figure out a way how to shape the 
  // spectrum and amplitude distribution of some noise independently by using a combination of 
  // filtering and waveshaping....

  int numSamples = 50000;
  int numBins    = 61;
  double b0      = 1.0;    // filter coeff for direct signal
  double b1      = 1.0;    // filter coeff for delayed input signal
  double a1      = 0.8;    // filter coeff for delayed output signal
  int numConvs   = 8;      // determines accuracy of our approximation of the pdf


  // Create the raw input noise:
  rsNoiseGenerator<double> ng;
  using Vec = std::vector<double>;
  int N = numSamples;
  Vec x(N);
  for(int n = 0; n < N; n++)
    x[n] = ng.getSample();
  //plotHistogram(x, numBins, -1.0, +1.0);  // uniform in -1..+1

  // alternatively, create a sinewave:
  //createSineWave(&x[0], numSamples, 440, 1.0, 44100);
  //createWaveform(&x[0], numSamples, 1, 440., 44100., 0.0, true);


  // Apply a 2-point MA filter y[n] = b0*x[n] + b1*x[n-1]:
  double xOld = 0;
  double yOld = 0;
  for(int n = 0; n < N; n++)
  {
    double y = b0*x[n] + b1*xOld + a1*yOld;
    xOld = x[n];
    yOld = y;
    x[n] = y;
  }
  plotHistogram(x, numBins, -2.0, +2.0); // triangular in -2..+2, if b0 = b1 = 1, a1 = 0


  // Here, we predict the amplitude distribution based on the filter coeffs. When a1 == 0, i.e. the
  // filter is nonrecursive, we actually get an exact pdf and cdf. When the filter is recursive, we
  // obtain an approximation based on a truncated impulse response. The numConvs parameter 
  // determines, how many samples of the impulse response are taken into account. The impulse 
  // reponse of the filter is: h[0] = b0, h[1] = b1 + a1*b0, h[n] = a1*h[n-1] for n >= 2:
  double h0 = b0, hn = b1 + a1*h0;
  using PiecePoly = rsPiecewisePolynomial<double>; 
  PiecePoly p0  = PiecePoly::irwinHall(0, -rsAbs(h0), +rsAbs(h0));
  PiecePoly p1  = PiecePoly::irwinHall(0, -rsAbs(hn), +rsAbs(hn));
  PiecePoly pdf = p0.convolve(p1);
  if(a1 != 0) {
    for(int n = 1; n < numConvs; n++) {
      hn *= a1;
      PiecePoly pn = PiecePoly::irwinHall(0, -rsAbs(hn), +rsAbs(hn));
      pdf = pdf.convolve(pn); }}
  PiecePoly cdf = pdf.integral();
  plot(pdf);  // should match the histogram
  plot(cdf);  // should go from 0 to 1


  // Now that we know the cdf of the noise, we can apply this:
  // https://en.wikipedia.org/wiki/Inverse_transform_sampling
  // in order to shape the amplitude distribution into a uniform distribution between 0..1:
  auto f = [&](double x) -> double { return cdf.evaluate(x); };
  for(int n = 0; n < N; n++)
    x[n] = f(x[n]);
  plotHistogram(x, numBins, -0.5, 1.5);
  // yup - works - looks uniform in 0..1 
  // -now we should have noise that is not white but has a uniform amplitude distribution
  // -we can now further shape the amplitude distribution to taste

  // Map it to the range -1..+1
  for(int n = 0; n < N; n++)
    x[n] = 2*x[n] - 1;
  plotHistogram(x, numBins, -1.5, 1.5);
  rosic::writeToMonoWaveFile("NoiseWaveShaped.wav", &x[0], N, 44100); 

  // Observations:
  // -if b0 = b1 = 1, the filtered output has a triangluar pdf from -1..+1 
  // -if b0 = 1, b1 = 0.5 (or vice versa), we get a trapezoidal distribution from -1.5..+1.5
  // -i think the general rule is that we need to convolved two stretched/compressed uniform 
  //  distribtutions where the strech factors are the filter coeffs
  // -i guess, this generalizes to higher order filters - we just have to perform repeated 
  //  convolutions
  // -for IIR filters, we may have to obtain the impluse-response first to express it as an MA 
  //  filter and then do repeated convolutions using the so obtained non-recursive coeffs
  // -when numConvs is too low, we get an emphasis of extreme amplitude values - the distribution
  //  gets batman ears
  // -for a closer to 1, we need more numConvs to get a reasonable approximation of the cdf

  // ToDo:
  // -measure autocorrelation of result
  // -try waveshping first and then filtering
  // -try multiple stages of filter -> shape or shape -> filter
  // -try the process on other signals, like sinewaves, sawtooth waves, etc.
  // -alternatively, measure the amplitude distribution from an arbitrary signal
  // -try to improve the approximation of the amplitude distribution after the IIR filter by 
  //  convolving slightly wider 0th order irwin-hall distribtuions. the idea is to compensate for 
  //  the truncation after finitely many convolutions by using wider input distributions. assess 
  //  the quality of a low-order approximation by comparing with a high-order approximation, i.e. 
  //  one using a higher number of convolutions

  // See:
  // https://en.wikipedia.org/wiki/Probability_integral_transform
  // this is eactly what we need

  // also of interest is the inverse:
  // https://en.wikipedia.org/wiki/Inverse_transform_sampling

  // I think, in order to get a uniform distribution out of a filtered uniform noise, it may work
  // like this:
  // -filter the uniform noise and figure out the resulting amplitude distribution f_x, for 
  //  example, a 2-point MA filter turns a uniform distribution from -1..+1 into a triangular 
  //  distribution from -2..+2
  // -from f_x figure out the cumulative distribution F_x and apply it as waveshaper to the 
  //  filter's output - the new variable should be uniformly distributed between 0..1
  // -scale and shift back to the desired range -1..+1
  // I also think, we could first waveshape, then filter, if we use the inverse CDF instead. These
  // two operations: filtering and waveshaping, used together in an appropriate way should give 
  // rise to a framework to shape the spectral and amplitude distribution properties of noise to
  // taste. It may be interesting to explore how they could be chained to produce more interesting
  // noises.

  // Question: can we derive a formula or algorithm that tells us the resulting amplitude 
  // distribution for a given filter (given that input to the filter is uniformly distributed?). 
  // Ideally, if possible, from the filter coeffs but maybe we need to consider the freq-response. 
  // For MA filters of length k, we get the Irvin-Hall distribution of the same order (maybe 
  // divided by k) - but what about recursive filters? Maybe we can start from MA filters with 
  // weights and express the impulse response of a recursive filter as such. Let's consider the
  // filter y[n] = b0*x[n] + b1*x[n-1] ...i think, the resulting amplitude distribution is the 
  // uniform distribution scaled (in x-direction, i.e. stretched or compressed) by b0 convolved 
  // with the uniform distribution stretched by b1...is that right? Consider the filter:
  //   y[n] = b0*x[n] + b1*x[n-1]
  // and let's assume that x has a probability density function (pdf) p_x(x). What is the pdf of
  // the output p_y(y)? If we assume that the x[n] are independently drawn from p_x(x), i.e. 
  // uncorrelated, i think, it's p_y(y) = conv(p_x(x/b0)/b0, p_x(x/b1)/b1), i.e. the convolution of
  // appropriately stretched/compressed and scaled versions of p_x(x). When any of b-coeffs 
  // approach zero, the corresponding factor in the convolution product approaches the Dirac 
  // delta distribution (small b make the distribution narrower and scale it up). When using more
  // past x samples, i.e. a higher order non-recursive filter, we would have to do repeated 
  // convolutions...right? or wrong? after the first convolution, one of the inputs to the 2nd 
  // convolution is correlated, so the independence assumption does not hold anymore for one of the
  // factors (is it enough, if it holds for only one?)..hmm...let's assume, it works - then, a 
  // recursive filter could be expressed via its impulse response and then the formula for the 
  // non-recursive filter could be used (leading to an infinitely often iterated convolution)

  // To actually compute such a pdf, we need a way of convolving functions that are expressed as
  // piecewise polynomials. 
}

// maybe let it take a parameter for the length and produce various test signals with various 
// lengths to see, how the quality depends on the length:
void blit()  
{
  double f  = 1001;    // signal frequency
  double fs = 44100;   // sample rate
  double length = 0.1; // length in seconds
  double period = fs/f;

  //// settings that lead to artifacts with the linBlep:
  //double inc = 1. / 27;  // inc = f/fs
  //inc = 29. / (93*2);
  //period = 1/inc;
  //f  = inc*fs;


  int N = (int) (fs*length);

  int blitLength = 30;     // blit length in samples...16 seems to be a reasonable default value

  double amplitude = 0.5;
  //int N = 120;
  //double period = 10.25;

  rsTableLinBlep<double, double> linTableBlit;
  linTableBlit.setLength(blitLength);

  rsTableMinBlep<double, double> minTableBlit;
  minTableBlit.setLength(blitLength);


  // naive and anti-aliased signals:
  std::vector<double> x(N), ylt(N), ymt(N); 

  /*
  // preliminary test:
  sbl.setLength(3);
  rsSetZero(y);
  sbl.addImpulse(0.75, 1);       // delay is 0.75, fractional position of spike is 0.25
  y[0] = sbl.getSample(1);
  for(int n = 1; n < N; n++)
    y[n] = sbl.getSample(0);
  rsStemPlot(y);
  */


  int numSpikes = 0;
  int nextSpike = 0;
  double ts = 0.0;  // exact time of (next) spike
  double tf = 0.0;  // fractional part of ts
  rsSetZero(ylt);
  linTableBlit.reset();
  for(int n = 0; n < N; n++) {
    if(n == nextSpike) {
      // generate naive impulse train:
      x[n] = amplitude;

      // generate bandlimited impulse train:
      linTableBlit.prepareForImpulse(1-tf, amplitude);
      ylt[n] = linTableBlit.getSample(amplitude);

      minTableBlit.prepareForImpulse(1-tf, amplitude);
      ymt[n] = minTableBlit.getSample(amplitude);

      // housekeeping for next spike:
      numSpikes++;

      //ts += period;           // accumulate, but...
      ts  = numSpikes * period; // ...this is numerically better than accumulation

      nextSpike = (int) ceil(numSpikes*period);
      tf = ts - (nextSpike-1);  // maybe clip to interval [0,1]
    }
    else {
      x[n] = 0;
      ylt[n] = linTableBlit.getSample(0);
      ymt[n] = minTableBlit.getSample(0);
    }
  }


  // compensate delay for easier comparison:
  rsArrayTools::shift(&ylt[0], N, -linTableBlit.getDelay());
  rsArrayTools::shift(&ymt[0], N, -minTableBlit.getDelay());

  //rsPlotVector(x);
  rsPlotVectors(x, ylt, ymt);
  //rsPlotVectors(x, ylt);

  //rosic::writeToMonoWaveFile("BlitTestNon.wav", &x[0],   N, int(fs));
  //rosic::writeToMonoWaveFile("BlitTestLin.wav", &ylt[0], N, int(fs));
  //rosic::writeToMonoWaveFile("BlitTestMin.wav", &ymt[0], N, int(fs));


  // Observations:
  // -when trying settings that lead to spurious spikes with the linBlep, nothing special happens 
  //  here...aside from the correction signal coming out as all-zeros everytime a spike occurs
  //  ...could this be a hint? ..well..no - that hanppens only when inc = 1./n for some integer n

  // -the quality does not really seem to increase with increasing order - look at the spectra with
  //  linear freq axis - the higher the order, the more side-maxima occur in the noise-floor
  //  -maybe that's because of the simple rectabgular window? ...try better windows
  // -on the plus side, even with length=1, we see a significant improvement over naive synthesis
  // -although - even it doesn't look like an improvement in the spectra when ramping up the order,
  //  it sounds cleaner - compare lengths 1 and 4 - 4 sounds like the sweet spot
  //  ...however, 2 actually also does sound worse than 1...something weird is going on...
  // -maybe try higher values for setTablePrecision - maybe the table is just too imprecise
  // -try also to remove the normalization of the corrector signal to (scaled) unit sum
  //  replace: rsScale(tempBuffer, amplitude / rsSum(tempBuffer));
  //  by:      rsScale(tempBuffer, amplitude);
  //  -> done: without normalization, results are much worse -> normalization is good!
  //  ...but maybe we could try to normalize ths sum-of-squares (energy)? but that would scale 
  //  everything down a lot, when the kernel length increases...hmmm..

  // -when plotting the blep with the hanning window, it starts below 0 and ends above 1 - we may
  //  have to rescale it (i think, it may be due to inaccuracy of numeric intergration - try if it 
  //  gets worse with less table precision)

  // -i think, for sincLength = 1, the rectangular window is actually better than Hanning - so the 
  //  choice of the window should depend on the sincLength?
  // -L=2 with Hanning looks pretty good but 4 looks worse - there doesn't seem to be a simple 
  //  relationship between length/window/quality ...may that depend also on signal frequency?
  //  -i guess, the amplitude of an aliased frequency depends on whether it falls on the maximum or
  //   minimum of the window's sidelobe pattern
  // -i think, we should use windows with high sidelobe rolloff rate for this application - it will
  //  confine the aliasig components into the upper frequency range where it is less annoying

  // -for the MinBlit, the overall amplitude seems a bit low - maybe that's due to not normalizing
  //  the sum of the corrector values? ...but i really don't want to introduce a temp-buffer for
  //  that purpose ...maybe for prodcution code, scrap the blit capability from the 
  //  LinBlitBlepBlamp object, too - i don't currently have practical use for blits anyway - but we 
  //  should keep the prototype implementation around, just in case it will be needed someday
  //  ...blit is an oddball anyway due to the not-a-function property in the continuous domain
  // ..or maybe leave the blit in the production version but without the normalization ...but above
  // we found that the results are much worse without normalization...hmm...we'll see...


  //GNUPlotter plt;
}

void blep()
{

  double fs  = 44100;    // sample rate
  //double inc = 19.0/256;  // phase increment per sample
  //double inc = 7.0/512;  // phase increment per sample
  //double inc = 30. / (93*3);
  //double inc = 1. / 19;
  double inc = GOLDEN_RATIO / 100;
  int N      = 800;      // number of samples to produce
  int shape  = 3;        // 1: saw, 2: square, 3: triangle
  int prec   = 20;       // table precision
  int length = 30;       // blep length

  // try to figure out the periodicity of the rippled cycles - it has to do with how many times we
  // have to loop through through the cycle unitl we are back at the sample branch of the blep
  // could it be something like gcd(num % den, den)

  //double inc = 5.0/100;  // phase increment per sample
  // 5./100 gives an access violation with the triangle wave in rsTableLinBlep::readTable - i 
  // think, we should give the tables one extra guard sample that repeats the last actual sample 
  // - i think we may then also get rid of the unelegant shouldReturnEarly function
  // ..for rsTableMinBlep, we don't get such an access violation
  // ok - i added a guard-sample to the table in rsTableLinBlep - shouldn't the MinBlep also need 
  // one? -> add unit test - be sure to test limiting cases (delayFraction == 0 and 1)

  //double inc = 3./100;  // phase increment per sample


  // Create osc and blit/blep/blamp objects:
  rsBlepReadyOsc<double> osc;
  osc.setPhaseIncrement(inc);
  //osc.setAmplitude(0.5);

  // make the osc output phase consistent with the output of our additive renderer:
  //if(shape == 3)
  //  osc.setStartPosition(0.25);
  //else
  //  osc.setStartPosition(0.5);
  //osc.reset();

  // new:
  if(shape == 3)  osc.reset(0.25);
  else            osc.reset(0.5);



  rsTableLinBlep<double, double> linTableBlep;
  linTableBlep.setTablePrecision(prec);
  linTableBlep.setLength(length);  

  rsTableMinBlep<double, double> minTableBlep;
  minTableBlep.setTablePrecision(prec);
  minTableBlep.setLength(length);

  rsPolyBlep1<double, double> polyBlep1;
  rsPolyBlep2<double, double> polyBlep2;


  // naive and anti-aliased signal
  std::vector<double> x(N), ylt(N), ymt(N), yp1(N), yp2(N); 
  // lt: linear-phase/table, mt: minimum-phase/table, p1/2: poly-blep1/2

  double stepDelay, stepAmp;

  for(int n = 0; n < N; n++)
  {
    switch(shape)
    {
    case 1: x[n] = osc.getSampleSaw(&stepDelay, &stepAmp);      break;
    case 2: x[n] = osc.getSampleSquare(&stepDelay, &stepAmp);   break;
    case 3: x[n] = osc.getSampleTriangle(&stepDelay, &stepAmp); break;
    }

    if(stepAmp != 0.0) { // a step (or corner) did occur
      if(shape != 3) {   // saw or square - prepare for steps
        linTableBlep.prepareForStep(stepDelay, stepAmp);
        minTableBlep.prepareForStep(stepDelay, stepAmp);
        polyBlep1.prepareForStep(stepDelay, stepAmp);
        polyBlep2.prepareForStep(stepDelay, stepAmp);
      } else {           // triangle - prepare for corners
        linTableBlep.prepareForCorner(stepDelay, stepAmp);
        minTableBlep.prepareForCorner(stepDelay, stepAmp);
        //polyBlep1.prepareForCorner(   stepDelay, stepAmp); // not yet implemented
        polyBlep2.prepareForCorner(stepDelay, stepAmp);
      }
    }
    // maybe the blep objects could figure out the size of the step itself - maybe do something like
    // blep.prepareForStep(stepDelay, x[n] - sbl.previousInputSample() );
    // ...hmm bute the previousInputSample in the delayline already may have a correction applied
    // ...i think that would be wrong - we need the difference of the new sample and the 
    // uncorrected previous input sample


    ylt[n] = linTableBlep.getSample(x[n]);
    ymt[n] = minTableBlep.getSample(x[n]);

    //rsAssert(fabs(ylt[n]) <= 0.7); // for debug

    yp1[n] = polyBlep1.getSample(x[n]);
    yp2[n] = polyBlep2.getSample(x[n]);
  }

  // apply amplitude scaler of 0.5 to all signals:
  x   = 0.5 * x;
  ylt = 0.5 * ylt;
  ymt = 0.5 * ymt;
  yp1 = 0.5 * yp1;
  yp2 = 0.5 * yp2;
  // todo: implement *= operator for vector/scalar


  // Create a perfect bandlimited saw/square/triangle for reference:
  double f = fs*inc;
  std::vector<double> r(N);
  createWaveform(&r[0], N, shape, f, fs, 0.0, true);
  r = 0.5 * r;

  // delay compensation:
  rsArrayTools::shift(&ylt[0], N, -linTableBlep.getDelay()); // linBlep has about 5-times the delay of 
  rsArrayTools::shift(&ymt[0], N, -minTableBlep.getDelay()); // minBlep
  rsArrayTools::shift(&yp1[0], N, -polyBlep1.getDelay());
  rsArrayTools::shift(&yp2[0], N, -polyBlep2.getDelay());

  //rosic::writeToMonoWaveFile("BlepTestNoAA.wav",   &x[0],   N, int(fs));
  //rosic::writeToMonoWaveFile("BlepTestLinTbl.wav", &ylt[0], N, int(fs));
  //rosic::writeToMonoWaveFile("BlepTestMinTbl.wav", &ymt[0], N, int(fs));
  //rosic::writeToMonoWaveFile("BlepTestPoly1.wav",  &yp1[0], N, int(fs));
  //rosic::writeToMonoWaveFile("BlepTestPoly2.wav",  &yp2[0], N, int(fs));
  //rosic::writeToMonoWaveFile("BlepTestBL.wav",     &r[0],   N, int(fs));
  //rsPlotVector(x);
  //rsPlotVector(r-y);  // error-signal: reference minus blepped
  //rsPlotVectors(x, y);
  //rsPlotVectors(x, ylt);
  rsPlotVectors(x, ylt, ymt);
  rsPlotVectors(x, yp1, yp2); 
  //rsPlotVectors(x, ylt, r, r-ylt);
  //rsPlotSpectrum(y); // oh - no this doesn't work - it takes a spectrum as input..


  // Observations:
  // Bug: with inc = GOLDEN_RATIO / 100; we see weird spurious spikes with rsTableLinBlep
  //  inc = 29. / (93*2) - the spurious spikes occur at a distance of 93 samples (the 1st one at 
  //  94), 30./(93*2) and 28./(93*2) are also interesting - btw - this works also with blepLength=1
  //  other increments: 1./19, 1./22, 1./23, 1./24, 1./26, 1./27, 1./28
  // -increasing bufferSize by factor 2 doesn't help
  //  ...maybe try with a blit - with the blit, there don't seem to be problems
  // -with inc = 1./19, the downward steps are correct, the upward steps (except the initial one) 
  //  have an additional downward spike -> look at w´how upward and downward correctors differ
  // -with inc = 1./38, all steps (exept the initial) are wrong
  // -interestingly, we call prepareForStep at sample 0 and sample 1...that seems to be wrong, too
  //  but is probably unrelated
  // -the tempBuffer is clearly showing these additional spikes
  // -could it be related to reading out the table near the center (time = 0) where the *table* 
  //  itself has a jump?
  // -ok - yes - this seems indeed to be the culprit - adding
  //    if( i == halfLength*tablePrecision - 1 ) return 0;
  //  to readTable seems to fix it - but this not yet a proper solution - it's only applicable
  //  to the blep and even then it's questionable - maybe it would be better to not store the 
  //  residual but the blep itself - the problem is the linear interpolation around the 
  //  discontinuity in the table
  // -it jumps from +0.445 to -0.5 (with prec=20) ...would it help to use an odd number of samples
  //  and insert 0 at the discontinuity (creating a symmetric table)?
  // -adding blepTbl[ic] = 0; to updateTables gets rid of the spikes but derivative-discontinuities
  //  remain
  // -todo: try to store the blep itself instead of the residual and create residual in discrete
  //  time domain

  // replacing: rsScale(tempBuffer, TSig(0.5)*amplitude / rsMean(tempBuffer));
  // by:        rsScale(tempBuffer, amplitude );
  // seems to work better - so with the blep, normalizing seems not a good idea - this is 
  // surprising because with the blit, normalization gave better results
  // todo: try to pre-normalize the blep table (it doesn't range from 0 to 1) - done

  // -hmm...well, at least the time-domain signals show the typical expected ripple pattern - but 
  //  when plotting the spectra, the blep doesn't seem to help at all against aliasing
  //  verify the computation of stepDelay
  // -the saw has actually increased aliasing with the blep

  // -with inc = 7.0/256, the 4th, 11th, 18th etc. cycles are bad - generally 4 + 7*k where k is
  //  the cycle number when we use 1-osc.getStepDelay() - when we use osc.getStepDelay(), all 
  //  cycles are wrong
  // -with inc = 8.0/256, all cycles are wrong (with 1-..) - the applied corrector is always zero
  //  but should always produce the same ripples
  //  -it seems to be the case that at these samples shouldReturnEarly triggers - maybe we should 
  //   not return early? - delayFraction == 0 in these cases
  // -removing the first check in:
  //   return delayFraction < eps || delayFraction >= TTim(1) || sincLength == 0;
  //  fixes it
  // -ok, i changed stepDelay = T(1) - pos/inc; to stepDelay = pos/inc;  so we don't have to do the
  //  1-... here anymore - actually, that *should* be wrong but seems in fact to be right - why?

  // ...ok - it seems to work now - for the square at 7/256 (f ~= 1200 Hz @44.1), a length of 2 
  // seems to be already good enough -> try higher frequencies..

  // -varying TablePrecision 
  //  -at sincLength=16:
  //   - 5: large spikes in the error signal for rectangle
  //   -10: error below 0.01 everywhere - same at 20 ...maybe 16 is a good choice
  //  -at sincLength=32:
  //   -10: error at 0.002
  //   -20: error at 0.015
  //   -40: error at 0.015
  //    -> seems like 32 would be a good choice - should the tablesize be roughly in the same 
  //    ballpark as the sinc-length? ...it seems that only longer sincs can actually benefit from 
  //    more precise tables
  //  -maybe make a plot of the maximum error as function of sincLength and tablePrecision

  // -MinBlep has high overshoot - much more than LinBlep - but maybe that's to be expected due to 
  //  the steepness of the elliptic filter
  // -MinBlamp look weird for triangle waves - the delay does not fit, the corners look too rounded 
  //  and and the peak amplitude seems too high
  //  -spectrum also shows a lot aliasing - MinBlamp is apparently still very wrong!
  //  -using the numeric integral of the MinBlep instead of the filter's ramp-response doesn't help
  //   (results are very similar)
  //  ...ok - this seems to be fixed - the time-axis values were wrong (multiplied by 0.5) - time 
  //  signal still look strange (too high peak amplitude) but spectrum looks ok now
  // -todo: compare the MinBlep result with an appropriately oversampled version using the same
  //  elliptic filter - should give similar (or even same?) results
  // -todo: try a MinBlit

  // -what if we have multiple steps/corners within one sample?
  
  // maybe implement a hard-syncable blep oscillator for liberty

  // ..okay...that looks all pretty good already - next stop: polybleps(done)/blamps(todo) - then 
  // we have the full palette of what could ever be needed ...and then maybe indeed turn to 
  // liberty...or continue with the modal synth

}

//void blamp() {  }  // not yet implemented - some blamp tests are done int the blep function

void polyBlep()
{
  // Plots the polynomial blit, blep, blamp functions based on B-Spline polynomials as derived in 
  // my sage notebook. My derived coeffs are different from the ones shown in the papers
  // "Perceptually informed synthesis of bandlimited classical waveforms using integrated 
  //  polynomial interpolation" and "Rounding Corners with BLAMP" 
  // ...but the sage plots look exactly like the ones shown in the paper. Their "D" variable must
  // mean something else than the "x" in the sage notebook...or something.
  // I'm trying to re-derive their polynomials in the fractional delay "d" (which i got working 
  // well already) in order to generalize to even higher order B-Spline poly-blits/bleps/blamps and
  // maybe other polynomial blits.. too - also, they don't give formulas for the 1st order B-Spline 
  // blamp, which may be useful, too and will fall out of these derivations as by-product

  // B-Spline polynomials are defined by starting with a unit-height rectangular pulse from -1/2 to
  // +1/2 (and zero outside this interval) and repeatedly convolving that function by itself:
  // sage input:
  // B0 = piecewise([((-1/2,1/2),1)]) # 0th order B-Spline
  // B1 = B0.convolution(B0)          # 1st order B-Spline
  // B2 = B1.convolution(B0)          # 2nd order B-Spline
  // B3 = B2.convolution(B0)          # 3rd order B-Spline
  // B3  # output 3rd order B-Spline - this is our BLIT
  // sage output:
  // piecewise(
  // x|-->  1/6*x^3 + x^2 + 2*x + 4/3 on (-2, -1], 
  // x|--> -1/2*x^3 - x^2       + 2/3 on (-1,  0], 
  // x|-->  1/2*x^3 - x^2       + 2/3 on ( 0,  1], 
  // x|--> -1/6*x^3 + x^2 - 2*x + 4/3 on ( 1,  2]; x)

  // B-Spline BLIT coefficients:                // ..in section...
  double blitA[4] = { 4./3,  2.,  1.,  1./6 };  // A: -2..-1
  double blitB[4] = { 2./3,  0., -1., -1./2 };  // B: -1...0
  double blitC[4] = { 2./3,  0., -1.,  1./2 };  // C:  0...1
  double blitD[4] = { 4./3, -2.,  1., -1./6 };  // D:  1...2

  // The BLEP is then obtained by integrating the BLIT:
  // sage input:
  // blep = B3.integral()
  // blep
  // sage output:
  // piecewise(
  // x|-->  1/24*x^4 + 1/3*x^3 + x^2 + 4/3*x + 2/3 on (-2, -1], 
  // x|--> -1/8 *x^4 - 1/3*x^3       + 2/3*x + 1/2 on (-1,  0], 
  // x|-->  1/8 *x^4 - 1/3*x^3       + 2/3*x + 1/2 on ( 0,  1], 
  // x|--> -1/24*x^4 + 1/3*x^3 - x^2 + 4/3*x + 1/3 on ( 1,  2]; x)

  // B-Spline BLEP coefficients:
  double blepA[5] = { 2./3, 4./3,  1.,  1./3,  1./24 };  // A: -2..-1
  double blepB[5] = { 1./2, 2./3,  0., -1./3, -1./8  };  // B: -1...0
  double blepC[5] = { 1./2, 2./3,  0., -1./3,  1./8  };  // C:  0...1
  double blepD[5] = { 1./3, 4./3, -1.,  1./3, -1./24 };  // D:  1...2

  // The BLAMP is obtained by integrating the BLEP:
  // sage input:
  // blamp = blep.integral()
  // blamp
  // sage output:
  // piecewise(
  // x|-->  1/120*x^5 + 1/12*x^4 + 1/3*x^3 + 2/3*x^2 + 2/3*x + 4/15 on (-2, -1], 
  // x|--> -1/40 *x^5 - 1/12*x^4           + 1/3*x^2 + 1/2*x + 7/30 on (-1,  0], 
  // x|-->  1/40 *x^5 - 1/12*x^4           + 1/3*x^2 + 1/2*x + 7/30 on ( 0,  1], 
  // x|--> -1/120*x^5 + 1/12*x^4 - 1/3*x^3 + 2/3*x^2 + 1/3*x + 4/15 on ( 1,  2]; x)

  // B-Spline BLAMP coefficients:
  double blampA[6] = { 4./15, 2./3,  2./3,  1./3,  1./12,  1./120 };  // A: -2..-1
  double blampB[6] = { 7./30, 1./2,  1./3,  0.,   -1./12, -1./40  };  // B: -1...0
  double blampC[6] = { 7./30, 1./2,  1./3,  0.,   -1./12, +1./40  };  // C:  0...1
  double blampD[6] = { 4./15, 1./3,  2./3, -1./3,  1./12, -1./120 };  // D:  1...2

  // For the BLAMP residual, we could define the ideal ramp and subtract it from the blamp:
  // ramp = piecewise([((-oo,0),0), ((0,oo),x)])
  // res  = blamp - ramp
  // ..however, sage does not produce a nice single piecewise function but instead a difference of
  // two piecewise functions. But transforming BLEP and BLAMP into the corresponding residuals 
  // manually is easy: for the blep residual, we would just subtract 1 from the constant (i.e. x^0) 
  // coefficients in sections C,D from the blep coeffs (turning 1/2, 1/3 into -1/2, -2/3) and for 
  // the blamp residual, we would do the same with the x^1 coeffs. For the BLIT, it doesn't seem 
  // sensible to define a residual because that would have to contain the Dirac delta-distribution 
  // which is not really a function and awkward to deal with in the continuous time domain.

  // Evaluate the blit/blep/blamp piecewise polynomials between -2 and +2 for plotting:
  int N = 400;
  std::vector<double> x(N), yBlit(N), yBlep(N), yBlamp(N);
  RAPT::rsArrayTools::fillWithRangeLinear(&x[0], N, -2.0, 2.0);
  typedef RAPT::rsPolynomial<double> PL;
  int i = 0;
  while(x[i] <= -1.0) {
    yBlit[i]  = PL::evaluate(x[i], blitA,  3);
    yBlep[i]  = PL::evaluate(x[i], blepA,  4);
    yBlamp[i] = PL::evaluate(x[i], blampA, 5);
    i++;
  }
  while(x[i] <= 0.0) {
    yBlit[i]  = PL::evaluate(x[i], blitB,  3);
    yBlep[i]  = PL::evaluate(x[i], blepB,  4);
    yBlamp[i] = PL::evaluate(x[i], blampB, 5);
    i++;
  }
  while(x[i] <= 1.0) {
    yBlit[i]  = PL::evaluate(x[i], blitC,  3);
    yBlep[i]  = PL::evaluate(x[i], blepC,  4);
    yBlamp[i] = PL::evaluate(x[i], blampC, 5);
    i++;
  }
  while(i < N && x[i] <= 2.0) {
    yBlit[i]  = PL::evaluate(x[i], blitD,  3);
    yBlep[i]  = PL::evaluate(x[i], blepD,  4);
    yBlamp[i] = PL::evaluate(x[i], blampD, 5);
    i++;
  }

  // plot the functions:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0], &yBlit[0], &yBlep[0], &yBlamp[0]);
  plt.plot();
  // the plots look good - they resemble the ones in Fig.3 in the "Rounding corners..." paper
  // ..still, my polynomial coeffs are different from theirs - wtf is going on? is their "D" 
  // variable equal to -x? the "Perceptually informed..." paper says under Eq 14: "D is a real 
  // number that corresponds to the delay from the beginning (n=0) of the impulse response" 
  // ...but i think, that doesn't  work out, too...more research needed....
}

void superBlep()
{
  // Experimenting with the blep-based supersaw implementation...

  int N = 2000;
  double sampleRate = 48000;
  double freq = 200;

  int maxDensity = 8;

  // create and set up the osc object:


  std::vector<rsUint32> primes(maxDensity+1);  // +1 because of the difference-of-primes based algos
  rsFillPrimeTable(&primes[0], maxDensity+1);

  typedef rsRatioGenerator<double>::RatioKind RK;
  rsRatioGenerator<double> ratGen;
  ratGen.setPrimeTable(&primes);
  ratGen.setRatioKind(RK::rangeSplitOdd);
  //ratGen.setParameter1(1.0/GOLDEN_RATIO);
  ratGen.setParameter1(0.7);

  typedef rsBlepOscArray<double, rsBlepReadyOscBase<double>, rsPolyBlep1<double, double>> SBO;
  //typedef rsBlepOscArray<double, rsBlepReadyOscBase<double>, rsPolyBlep2<double, double>> SBO;
  //typedef rsBlepOscArray<double, rsBlepReadyOscBase<double>, rsTableMinBlep<double, double>> SBO;
  SBO osc(&ratGen);
  osc.setReferenceIncrement(freq / sampleRate);
  osc.setMaxDensity(maxDensity);
  osc.setDetune(0.07);
  osc.setDensity(3);
  osc.setStereoSpread(0.0);


  // produce output signal:
  std::vector<double> x(N);
  //for(int n = 0; n < N; n++)
  //  x[n] = osc.getSample();

  // trying to figure out stereo spread loudness bug - alternatingly call getSample and 
  // getSampleFrameStereo - with stereoSpread=0, the result should be the same - the even numbered
  // samples use the mono function, the odd numbered samples use the stereo function - if something
  // is wrong, we'll see the amplitude alternating between even and odd samples
  rsArrayTools::fillWithZeros(&x[0], N);
  for(int n = 0; n < N; n += 2) {
    double dummy = 0;
    x[n] = osc.getSample();
    osc.getSampleFrameStereo(&x[n+1], &dummy);
  }




  rsArrayTools::normalize(&x[0], N);


  // plot:
  rsPlotVector(x);
  //rosic::writeToMonoWaveFile("SuperSaw.wav", &x[0], (int)x.size(), (int)sampleRate);



  // Observations:
  // -we should try to choose the ratios in such a way that the minimum and maximum of the envelope
  //  are not too far apart (of the sustained section - there will be transient envelope peak at 
  //  the start, if the oscs start in phase)
  // -using 1.f where f is the fractional part of the metallic ratios works better than using the
  //  metallic ratios themselves

  // -using the naive version (that aliases as hell) and slapping a highpass on it may actually be
  //  useful - the aliasing creates a somewhat pleasant hiss that makes the whole sound brighter
  //  (3rd order butterworth sounds good at f=1kHz, d=0.1) - but maybe we can create the hiss
  //  differently?
  //  -maybe it sounds brighter because the blep is introducing some lowpass charachter
  //   in the passband? -> try LinBlep/MinBlep -> ok - yes - that seems to be part of the reason
  //  -the highpassed naive one sounds still brighter/different
  //  -a high-shelving filter helps somwhat
  //  -maybe PolyBlep1 is better than PolyBlep2 


  // todo: distribute the start-phases - maybe make the start-phase spread a user parameter
  // -call the number of oscs "density" - maybe the density can be made continuous by somehow 
  //  fading in/out the new frequencies? (we would also have to somehow change the 
  //  freq-distribution continuously - maybe fade in the new freq in the center while shifting
  //  the other freqs toward the sides)
  // -apply amplitude bell curve according to frequency (make inner freqs louder)

  //
  // -i tried the naive supersaw at 2kHz at 44.1 and 48 kHz sample-rate - it seems to make the 
  //  aliasing artifacts even worse when the samplerate is a multiple of the fundamental - if it's
  //  unrelated, the aliasing tends to become more noisy which is better (but still very bad)

  // -whatever anti-aliasing is used (with 2kHz saw), when playing it via headphones, there seems 
  //  to be some sort of "hum" - but i don't see anything in the spectrum - what's going on?
  //  ...oh- it becomes visible when using a smaller FFT size (like 512) - apparently, it averages
  //  out over long enough time or something
  // -could i be somehow hearing the envelope? ...if that's the case, how about applying a 
  //  leveller that somehow cancels the envelope? ...that might be a good idea anyway
  // -there definitely is some sort of audible "flutter" or something - the character can be 
  //  modified by dynamics compression

  // Conclusion:
  // -use PolyBlep1
  // -use optional 3rd order butterworth highpass tuned to fundamental (or maybe other order - look
  //  at the original
  // -don't bother trying to tune the sample-rate to the saw-freq

  // hmm - the golden ratio is only hardest to approximate by continued fractions in the sense that 
  // you need a lot of CF coefficients - but, of course, you don't need higher denominators - but 
  // i think, the size of the denominator is what is really relevant for us here with respect to
  // assessing the irrationality
  // ...what we really want is to maximize the time unitl which the phases (approxiamtely) line
  // up again ...maybe between any pair of frequencies...

  // some suggestions from aciddose:
  // https://www.kvraudio.com/forum/viewtopic.php?p=5444677#p5444677
  // detune = detune + depth * (sqrt(prime) - sqrt(prime))
  // detune = depth * (sqrt(prime) - sqrt(prime))
  // detune = depth * (sqrt(i) - sqrt(i - 1))
  // ...differences of square-roots of primes or integers. square roots of integers are either
  // themselves integers or irrational. there are no non-integer rational square-roots of integers

  // the fractional parts of metallic ratios sound a bit smoother than the primeSqrts

  // essentially, with all thes algorithms
  // maybe in order to have the freq distributed roughly evenly, we should start with an arithmetic
  // progression and then add some irrational offsets?
  // or: recursively split the interval between fl, fu at a "midpoint" that is not exactly in the 
  // middle but at a user defined ratio - maybe apply the ratio and 1-ratio alternatingly to avoid
  // skeing the whole distribution ...iirc, that was what my old algos were also doing

  // todo: tune the freq-ratios by ear: implement a superosc in ToolChain, set it to 3 oscs - the 
  // outer two have fixed freqs, tune the inner one until it sounds best, then introduce a 4th saw,
  // and tune that until it sounds best and so on - bring in one saw at a time and leave the 
  // already tuned ones alone when tuning the new one - this will us later allow continuous density
  // by just ramping up the additional saw gradually
  // ..but: what sort of start-phase configuration whould we use for tuning - all at zero or
  // maximally spread out?

  // -the range-splitting results seem to be somewhat worse that the metallic ratio fractional 
  //  parts (mrf) - looking at the spectrum, it seems that mrf has more partials concentrated 
  //  around the center 
  //  -> we should have an option to somehow contract the ratios toward the center - maybe start
  //  with the distribution from range splitting but then apply a monotonic mapping function - 
  //  maybe something like the bipolar rational map ...maybe the unipolar version can be used
  //  to introduce skew - maybe have pre-skew -> contract/spread -> post-skew
}

void superSawDensitySweep() // rename to SuperSawFeedbackFM
{
  // We try to implement a continuos sweep of the density. For simplicity, we use an array of 
  // rsBlepSawOsc. Production code should probably optimize this (share the blep object, etc.).

  int maxDensity = 7;
  int N = 1000000;          // number of samples

  double fMin = 95;
  double fMax = 105;
  double sampleRate = 44100;
  double p = 0.5;  // shape of distribution of frequencies

  double modAmount = 0.0;

  std::vector<rsBlepSawOsc<double>> oscs(maxDensity);

  // filter for modulating the frequencies
  rosic::TwoPoleFilter flt;
  flt.setBandwidth(0.25); // maybe that should be set in Hz, not octaves
  flt.setFrequency(0.25*(fMin+fMax));
  flt.setMode(flt.BANDPASS);
  flt.setSampleRate(sampleRate);

  // test - altenrtively use a subsonic lowpass:
  flt.setMode(flt.LOWPASS12);
  flt.setFrequency(7.0);


  using Vec = std::vector<double>;
  using AT  = rsArrayTools;

  Vec freqs(maxDensity), amps(maxDensity);  // frequencies and amplitudes

  Vec y(N);

  for(int i = 0; i < maxDensity; i++)
    oscs[i].reset(i / double(maxDensity));

  double mod = 0;  // modulation signal

  for(int n = 0; n < N; n++)
  {
    int D = maxDensity;  // instantaneous density - preliminary

    AT::fillWithRange(&freqs[0], D, fMin, fMax, p);
    AT::fillWithRange(&amps[0],  D, 1.0,  1.0,  1.0); // preliminary - all 1

    y[n] = 0;
    for(int i = 0; i < D; i++)
    {
      double fi = freqs[i];

      if(rsIsEven(i)) fi += modAmount*mod*fi;
      else            fi -= modAmount*mod*fi;

      oscs[i].setPhaseIncrement(fi/sampleRate);

      // test: invert every other saw:
      if(rsIsEven(i)) y[n] += oscs[i].getSample();
      else            y[n] -= oscs[i].getSample();

      //y[n] += oscs[i].getSample();  // without inversion
    }

    //mod = flt.getSample(y[n]);   // for bandpass
    mod = flt.getSample(y[n] * y[n]);  // for lowpass
  }

  AT::normalize(&y[0], N);
  //rosic::writeToMonoWaveFile("SuperSawDensitySweep.wav", &y[0], N, sampleRate);
  rosic::writeToMonoWaveFile("SuperSaw.wav", &y[0], N, (int)sampleRate);


  std::cout << "Done";
  //rsPlotVector(y);
}
// idea to randomize the frequencies:
// -pass output through bandpass, tuned to center frequency
//  ..or an octave below is also very interesting
// -use that output as LFO for the oscillators - alternatingly up and down
// -alternatively to a bandpass: use a lowpass with subsonic cutoff freq on the squared output
// -maybe as an alternative spacing strategy: strat from equidistant frequencies and offset them
//  alternatingly by some user adjustable amount

void superSawStereo()
{
  // Testing the stereo spread algorithm by creating a stereo signal using getSampleFrameStereo 
  // mixing it to mono and comparing it to what would be generated using getSample (the results 
  // should be the same)

  int N = 2000;

  //rosic::rsOscArrayPolyBlep1 osc;
  TestOscArrayPolyBlep osc;
  osc.setMaxDensity(8);    // so the arrays don't get too long in the debugger
  osc.setSampleRate(44100);
  osc.setFrequency(300*GOLDEN_RATIO);
  osc.setDensity(4);
  osc.setStereoSpread(1);
  osc.setDetunePercent(5.0);
  osc.setInitialPhaseCoherence(1);


  // produce output signals:
  std::vector<double> xMono(N), xLeft(N), xRight(N), xMid(N), xSide(N);

  // produce mono signal:
  osc.reset();
  for(int n = 0; n < N; n++)
    xMono[n] = osc.getSample();

  // produce stereo signals and ontain mid/side signals from them:
  osc.reset();
  rsArrayTools::fillWithZeros(&xLeft[0],  N);
  rsArrayTools::fillWithZeros(&xRight[0], N);
  for(int n = 0; n < N; n++)
    osc.getSampleFrameStereo(&xLeft[n], &xRight[n]);
  xMid  = 0.5 * (xLeft + xRight);
  xSide = 0.5 * (xLeft - xRight);
  // times 0.5, because it works as follows: when stereoSpread is set to zero, the left and right 
  // outputs are equal to output that the mono version produces - so when we add the two, we get a 
  // factor of two which we should compensate for here

  //rsPlotVectors(xMono, xMid, xMono-xMid);  // xMid should be equal to xMono
  //rsPlotVectors(xLeft, xRight);

  // we want the stereo spread work as follows - let stereoSpread = 1, for increasing densities, the 
  // left and right amplitude arrays should look like (underscore _ stands for: doesn't matter):
  // D = 1: L = 1,_,_,_,_,_,_,_   R = 1,_,_,_,_,_,_,_
  // D = 2: L = 2,0,_,_,_,_,_,_   R = 0,2,_,_,_,_,_,_
  // D = 3: L = 2,1,0,_,_,_,_,_   R = 0,1,2,_,_,_,_,_
  // D = 4: L = 2,0,2,0,_,_,_,_   R = 0,2,0,2,_,_,_,_
  // D = 5: L = 2,0,1,2,0,_,_,_   R = 0,2,1,0,2,_,_,_
  // D = 6: L = 2,0,2,0,2,0,_,_   R = 0,2,0,2,0,2,_,_
  // D = 7: L = 2,0,2,1,0,2,0,_   R = 0,2,0,1,2,0,2,_
  // D = 8: L = 2,0,2,0,2,0,2,0   R = 0,2,0,2,0,2,0,2
  // and for stereoSpread = -1, the roles of L and R should be reversed - these numbers are only 
  // valid for a linear pan - but currently an (approximate) contant power law is used

  osc.initAmpArrays(); osc.setDensity(1);  // ok
  osc.initAmpArrays(); osc.setDensity(2);  // ok
  osc.initAmpArrays(); osc.setDensity(3);  // ok
  osc.initAmpArrays(); osc.setDensity(4);  // ok
  osc.initAmpArrays(); osc.setDensity(5);  // ok
  osc.initAmpArrays(); osc.setDensity(6);  // ok
  osc.initAmpArrays(); osc.setDensity(7);  // ok
  osc.initAmpArrays(); osc.setDensity(8);  // ok
  //int dummy = 0;
}

// computes algorithm parameters for the two-pice oscillator from user parameters
// pulse-width: -1..+1, sawVsSquare: -1..+1
void getTwoPieceAlgoParams(double pulseWidth, double sawVsSquare,
  double& h, double& a1, double& b1, double& a2, double& b2)
{
  typedef RAPT::rsBlepReadyOscBase<double> Osc;

  double d = 0.5; // preliminary - todo: compute from pulseWidth
  h = 0.5; // preliminary - todo: compute from pulse-width


  double a1u, b1u, a2u, b2u, a1d, b1d, a2d, b2d;
  a1u = a2u = +2.0;      // 1st and 2nd segment slopes for saw-up
  a1d = a2d = -2.0;      // 1st and 2nd segment slopes for saw-down
  b1u = b2u = -1.0 - d;  // 1st and 2nd segment offsets for saw-up
  b1d = b2d = +1.0 - d;  // 1st and 2nd segment offsets for saw-down
  // shifted saw-up:    x(t) =  2*(t-d/2) - 1 =  2*t - d - 1
  // shifted saw-down:  x(t) = -2*(t+d/2) + 1 = -2*t - d + 1


  double s = 0.5 * (sawVsSquare + 1); // convert shape parameter from -1..+1 to 0..1
  //s = 0; // test
  double wu = 1-s, wd = s;            // weights for up/down coeffs
  a1 = wu*a1u + wd*a1d;
  a2 = wu*a2u + wd*a2d;
  b1 = wu*b1u + wd*b1d;
  b2 = wu*b2u + wd*b2d;
}


void getTwoPieceAlgoParams2(double h, double slope, double slopeDiff, 
  double& a1, double& b1, double& a2, double& b2)
{
  double s = slope;
  double d = slopeDiff;

  a1 = (s+d) / h;
  a2 = (s-d) / h; 
  // how about using 1-h as denominator in the 2nd equation as alternative? ..and maybe an overall 
  // factor of 2? ...that wouldn't change any ofe the computations below, btw. - or a normalization
  // by peak amplitude rather than energy? maybe the formulas would come out simpler? maybe provide
  // different versions of the user-to-algo parameter mapping (also including the mix-of-two-saws 
  // approach) and make them static member functions of rsBlepReadyOsc(Base), so client code can 
  // select, which control scheme it wants to use


  // for a1 = +-2, we want b1 to be zero - maybe that can be used as additional constraint instead 
  // of the energy?


  double a12 = a1*a1, a22 = a2*a2;
  double h2 = h*h, h3 = h2*h, h4 = h2*h2, h5 = h3*h2;

  // energy normalized to 1/2:
  double r = sqrt(3*(a12-a22)*h5 - 18*a22*h3 - 3*(a12-4*a22)*h4 + 6*(2*a22-3)*h2 - 3*(a22-6)*h);
  // inside the sqrt is a 5th order polynomial in h - maybe evaluate via Horner scheme

  // 1st solution:
  b1 = -1./6 * (3*a1*h2        + r) / h;
  b2 = -1./6 * (3*a2*h2 - 3*a2 + r) / (h-1);

  //// 2nd solution:
  //b1 = -1./6 * (3*a1*h2        - r) / h;
  //b2 = -1./6 * (3*a2*h2 - 3*a2 - r) / (h-1);

  //// new: energy normalized to 1/4:
  //double r = sqrt(3*(a12-a22)*h5 - 18*a22*h3 - 3*(a12-4*a22)*h4 + 3*(4*a22-3)*h2 - 3*(a22-3)*h);
  //// inside the sqrt is a 5th order polynomial in h - maybe evaluate via Horner scheme
  //b1 = (-1./6) * (3*a1*h2        + r) / h;
  //b2 = (-1./6) * (3*a2*h2 - 3*a2 + r) / (h-1);

  //// energy normalized to 1/12:
  //r = sqrt(3*(a12-a22)*h5 - 18*a22*h3 - 3*(a12-4*a22)*h4 + 3*(4*a22-1)*h2 - 3*(a22-1)*h);
  //b1 = (-1./6) * (3*a1*h2        + r) / h ;
  //b2 = (-1./6) * (3*a2*h2 - 3*a2 + r) / (h-1);
  // nnnaaahhh!! it's all wrong! why?


  // the sawtooth is still wrong - maybe 1/4 is not the correct target value? integrate 
  // f(x) = 2*x from 0 to 1/2, take times two:
  // 2 * integral(x^2, x, 0, 1/2)  ->  1/12 
}
/*
given a1,a2 and applying the constraints of no DC and energy = 1/2 (maybe we should use 1/4?) gives
two solutions for b1,b2 via sage:

b1 == -1/6*(3*a1*h^2        + sqrt(3*(a1^2 - a2^2)*h^5 - 18*a2^2*h^3 - 3*(a1^2 - 4*a2^2)*h^4 + 6*(2*a2^2 - 3)*h^2 - 3*(a2^2 - 6)*h))/h, 
b2 == -1/6*(3*a2*h^2 - 3*a2 + sqrt(3*(a1^2 - a2^2)*h^5 - 18*a2^2*h^3 - 3*(a1^2 - 4*a2^2)*h^4 + 6*(2*a2^2 - 3)*h^2 - 3*(a2^2 - 6)*h))/(h-1) 

b1 == -1/6*(3*a1*h^2        - sqrt(3*(a1^2 - a2^2)*h^5 - 18*a2^2*h^3 - 3*(a1^2 - 4*a2^2)*h^4 + 6*(2*a2^2 - 3)*h^2 - 3*(a2^2 - 6)*h))/h, 
b2 == -1/6*(3*a2*h^2 - 3*a2 - sqrt(3*(a1^2 - a2^2)*h^5 - 18*a2^2*h^3 - 3*(a1^2 - 4*a2^2)*h^4 + 6*(2*a2^2 - 3)*h^2 - 3*(a2^2 - 6)*h))/(h-1)

// 

b1 == -1/6*(3*a1*h^2        + sqrt(3*(a1^2 - a2^2)*h^5 - 18*a2^2*h^3 - 3*(a1^2 - 4*a2^2)*h^4 + 3*(4*a2^2 - 1)*h^2 - 3*(a2^2 - 1)*h))/h 
b2 == -1/6*(3*a2*h^2 - 3*a2 + sqrt(3*(a1^2 - a2^2)*h^5 - 18*a2^2*h^3 - 3*(a1^2 - 4*a2^2)*h^4 + 3*(4*a2^2 - 1)*h^2 - 3*(a2^2 - 1)*h))/(h-1)

b1 == -1/6*(3*a1*h^2        - sqrt(3*(a1^2 - a2^2)*h^5 - 18*a2^2*h^3 - 3*(a1^2 - 4*a2^2)*h^4 + 3*(4*a2^2 - 1)*h^2 - 3*(a2^2 - 1)*h))/h
b2 == -1/6*(3*a2*h^2 - 3*a2 - sqrt(3*(a1^2 - a2^2)*h^5 - 18*a2^2*h^3 - 3*(a1^2 - 4*a2^2)*h^4 + 3*(4*a2^2 - 1)*h^2 - 3*(a2^2 - 1)*h))/(h-1)


the input was:

var("x a1 b1 a2 b2 h E")
y1(x) = a1*x + b1  # 1st segment
y2(x) = a2*x + b2  # 2nd segment
A = integral( y1(x),    x, 0, h) + integral( y2(x),    x, h, 1) # signed area
E = integral((y1(x))^2, x, 0, h) + integral((y2(x))^2, x, h, 1) # energy
# gives:
# A = 1/2*a1*h^2 - 1/2*a2*h^2 + b1*h - b2*h + 1/2*a2 + b2
# E = 1/3*a1^2*h^3 - 1/3*a2^2*h^3 + a1*b1*h^2 - a2*b2*h^2 
#     + b1^2*h - b2^2*h + 1/3*a2^2 + a2*b2 + b2^2
# A should be zero (no DC) and E should be 1/2 (energy of square-wave with
# amplitude 1/2), so we get two equations and solve them for b1,b2:
eq1 = 0   == 1/2*a1*h^2 - 1/2*a2*h^2 + b1*h - b2*h + 1/2*a2 + b2
eq2 = 1/2 == 1/3*a1^2*h^3 - 1/3*a2^2*h^3 + a1*b1*h^2 - a2*b2*h^2 + b1^2*h - b2^2*h + 1/3*a2^2 + a2*b2 + b2^2
solve([eq1,eq2],[b1,b2])

can we get rid of the sqrt (and maybe get a generally simpler formula) by normalizing the amplitude 
differently - for example with respect to peak amplitude? or area of absolute values?

maybe we could say that the b1,b2 coeffs should be computed from a1,a2 via a linear (or affine) 
transforn and determine the coeffs of the transform form some constraints, like: if a1 = 2, then b1 
should come out as 0, etc. - with enough such constraints, we could fix the coeffs of the trafo

...or maybe we should use a1,b1 as inputs and compute a2,b2? ...soo many possibilities to devise a 
control scheme - maybe we should implement them all and let client code pick one - or even allow 
the user to select one in a combobox

*/


void twoPieceOsc()
{
  // Experiments with an oscillator waveform consisting of two linear segments

  // x1(t) = a1 * t + b1    for t in 0...h
  // x2(t) = a2 * t + b2    for t in h...1
  //
  // this allows for sawtooth (up and down), triangle and square waves (h=0.5):
  // saw-up:    x1(t) = x2(t) =  2*t - 1
  // saw-down:  x1(t) = x2(t) = -2*t + 1
  // square:    x1(t) = -1, x2(t) = +1
  // triangle:  x1(t) = 4*t - 1, x2(t) = -4*t + 3
  //
  // a pulse-wave can be obtained by using a saw-up and saw-down and adding them up. one of them 
  // has to be shifted by half a cycle to obtain a square - other shifts will give other 
  // pulse-widths. problem: if they are exactly in phase, they will sum to zero, good: they will 
  // naturally produce a DC-free pulse-wave (because none of them has any DC). we could provide two 
  // user paremeters: one for the shift (determining the pulse-width) and one for the mix of the 
  // two saws (-1: saw-down, 0: square/pulse, +1: saw up) - problem: at -1 and +1 the "pulse-width"
  // paremeter would do nothing
  // potential problem when modulating the midpoint h: we may miss a blep - solution: keep the old 
  // h and update to the new value in getSample - with the new and the old value, we can check 
  // against both and if there's a transition for any of the two, we can determine the exact 
  // instant of the step (maybe by some line equation or something)

  // We may use a mix of phase-shifted saw-up and saw-down waves to obtain all sorts of pulse
  // waves. The equation is
  // y1(t) =  2 * wrap(t - d/2) - 1    saw-up
  // y2(t) = -2 * wrap(t + d/2) + 1    saw-down
  // y(t)  = (1-b)*y1 + b*y2           mix-of-saws
  // where "wrap" means wrap-into-interval 0..1, see also here:
  // https://www.desmos.com/calculator/wbot9ibar5
  // Where we have two user parameters d and b. d determines the phase-offset between the two saws 
  // and is applied symettrically, halfway to each of the saws. b is the mixing cofficient. At b=0.5
  // There are two transition points (step discontinuities) at d/2 and 1-d/2.

  // But maybe it's better to let the user set h and the (appropriately scaled) slope at time zero.

  int numCycles       = 21;
  int samplesPerCycle = 100;
  int numSamples      = numCycles * samplesPerCycle;
  //int N = 1000;
  //double inc = 0.01;
  double inc = 1.0 /  samplesPerCycle;

  RAPT::rsBlepReadyOsc<double> osc;
  osc.setPhaseIncrement(inc);

  double a1, b1, a2, b2, h;
  h = 0.5;
  //a1 = a2 = +2; b1 = b2 = -1;      // saw up
  a1 = a2 = -2; b1 = b2 = +1;      // saw down
  //a1 = 0; b1 = -1; a2 = 0; b2 = 1; // square (low first) - but here it starts high (start-phase wrong?)

  std::vector<double> xc(numSamples);  // (c)ontinuous parameter change
  double dummy;             // for blep/blamp info - not yet used

  //double k = sqrt(2.);
  double k = 2.0;

  double h1 =  0.5, h2 =  0.5;  // transition time at start and end
  double s1 =  k,   s2 = -k;    // slope sum at start and end
  double d1 =  0,   d2 =  0;    // slope difference at start and end

  // continuous parameter sweep:
  for(int n = 0; n < numSamples; n++) {
    double p = n / (numSamples-1.0);
    double h = rsLinToLin(p, 0.0, 1.0, h1, h2);
    double s = rsLinToLin(p, 0.0, 1.0, s1, s2);
    double d = rsLinToLin(p, 0.0, 1.0, d1, d2);
    getTwoPieceAlgoParams2(h, s, d, a1, b1, a2, b2);
    xc[n] = osc.getSampleTwoPiece(&dummy, &dummy, &dummy, h, a1, b1, a2, b2);
  }
  // todo: don't change the parameters continuously but per cycle - have a total of 21 cycles for 
  // steps of 0.1 - or maybe have two plots - one with continuous change, one with per-cycle change
  // algo may produce NaN when d has large absolute value (for h= 0.5 and s = 0)

  // parameters jump at the cycle boundaries:
  std::vector<double> xd(numSamples);  // (d)iscrete parameter change
  for(int cycle = 0; cycle < numCycles; cycle++)
  {
    double p = cycle / (numCycles-1.0);
    double h = rsLinToLin(p, 0.0, 1.0, h1, h2);
    double s = rsLinToLin(p, 0.0, 1.0, s1, s2);
    double d = rsLinToLin(p, 0.0, 1.0, d1, d2);
    getTwoPieceAlgoParams2(h, s, d, a1, b1, a2, b2);
    int n0 = cycle * samplesPerCycle;
    for(int n = 0; n < samplesPerCycle; n++)
      xd[n0 + n] = osc.getSampleTwoPiece(&dummy, &dummy, &dummy, h, a1, b1, a2, b2);
  }

  // saw ist still wrong (shows additional steps) - presumably because the energy should be 
  // normalized to 1/4, not 1/2 - overall amplitude of square should be 1/2
  // -> nope - that doesn't work either

  //// user parameters at start and end:
  //double pw1 = -1, pw2 = +1;    // pulse-width
  //double sh1 =  0, sh2 =  0;    // shape: -1: saw-down, 0: square/pulse, +1: saw-up

  //for(int n = 0; n < N; n++)
  //{
  //  double p = n / (N-1.0);
  //  double pw = rsLinToLin(p, 0.0, 1.0, pw1, pw2);
  //  double sh = rsLinToLin(p, 0.0, 1.0, sh1, sh2);
  //  getTwoPieceAlgoParams(pw, sh, h, a1, b1, a2, b2);
  //  x[n] = osc.getSampleTwoPiece(&dummy, &dummy, &dummy, h, a1, b1, a2, b2);
  //}


  // todo: create a square/pulse wave by using a mix of saw-up/saw-down (for square, they must have
  // a phase-shift of 90°)

  //rsPlotVector(xc);
  rsPlotVector(xd);
  //rsPlotVectors(xd, xc);
}


// maybe move to test-signal generation
template<class TOsc>
std::vector<double> syncSweep(TOsc& osc, double masterFreq, double slaveStartFreq, 
  double slaveEndFreq, double lengthInSeconds, double sampleRate, bool withBlep)
{
  osc.setMasterIncrement(masterFreq     / sampleRate);
  osc.setSlaveIncrement( slaveStartFreq / sampleRate);
  osc.reset();
  int numSamples = (int) round(lengthInSeconds * sampleRate);
  std::vector<double> y(numSamples);
  for(int n = 0; n < numSamples; n++) {
    double freqMix = n/(numSamples-1.0);  // sweep freq of slave osc...
    double slaveFreq = (1-freqMix) * slaveStartFreq + freqMix * slaveEndFreq;
    osc.setSlaveIncrement(slaveFreq / sampleRate);
    if(withBlep)  y[n] = osc.getSample();
    else          y[n] = osc.getSampleNaive();
  }
  return y;
}
std::vector<double> downSample(const std::vector<double>& x, int decimationFactor,
  bool antiAlias = false)
{
  int Nx = (int) x.size();
  int d  = decimationFactor;
  int Ny = Nx / d;
  std::vector<double> y(Ny);


  if(antiAlias) {
    RAPT::rsEllipticSubBandFilter<double, double> flt;
    flt.setSubDivision(d);
    for(int n = 0; n < Ny; n++) {
      y[n] = flt.getSample(x[n*d]);
      for(int k = 1; k < d; k++)
        flt.getSample(x[n*d+k]); // dummy outputs which are discarded
      // maybe have other anti-aliasing options (for example, zero-phase, no-overshoot, etc)
    }
  }
  else
    for(int n = 0; n < Ny; n++)
      y[n] = x[n*d];

  return y;
}

void syncSweep()
{
  // We create a hardsync sweep once naively, once via blep and once via oversampling to compare 
  // them

  //int N = 5000;
  double sampleRate     = 44100;
  double masterFreq     = 100 * GOLDEN_RATIO;  // master osc freq
  double slaveStartFreq = 20 * masterFreq;     // start freq of slave osc 
  double slaveEndFreq   = 1  * masterFreq;     // end   freq of slave osc 2
  double amplitude      = 0.25;                // amplitude
  double length         = 5.0;                 // length in seconds
  int    overSampling   = 16;                  // oversampling factor

  typedef rsSyncPhasor<double, rsPolyBlep1<double, double>> SyncOsc;
  //typedef rsSyncPhasor<double, rsPolyBlep2<double, double>> SyncOsc;
  //typedef rsSyncPhasor<double, rsTableMinBlep<double, double>> SyncOsc;
  SyncOsc osc;

  typedef std::vector<double> Vec; 
  Vec xNaive = syncSweep(osc, masterFreq, slaveStartFreq, slaveEndFreq, length, sampleRate, false);
  Vec xBlep  = syncSweep(osc, masterFreq, slaveStartFreq, slaveEndFreq, length, sampleRate, true);
  Vec xOver  = downSample(syncSweep(osc, masterFreq, slaveStartFreq, slaveEndFreq, length, 
    overSampling*sampleRate, false), overSampling, true);

  xNaive = xNaive - 0.5;
  xBlep  = xBlep  - 0.5;
  xOver  = xOver  - 0.5;

  rosic::writeToMonoWaveFile("SyncSweepNaive.wav", &xNaive[0], (int)xNaive.size(), (int)sampleRate);
  rosic::writeToMonoWaveFile("SyncSweepBlep.wav",  &xBlep[0],  (int)xBlep.size(),  (int)sampleRate);
  rosic::writeToMonoWaveFile("SyncSweepOver.wav",  &xOver[0],  (int)xOver.size(),  (int)sampleRate);
}

void syncPhasor()
{
  int N = 5000;
  double fs       = 44100;
  double f1       = 100 * GOLDEN_RATIO;  // master osc freq
  double f2_start = 20 * f1;             // start freq of slave osc 
  double f2_end   = 2  * f1;             // end   freq of slave osc 2
  double a        = 0.25;                 // amplitude

  //typedef rsSyncPhasor<double, rsPolyBlep1<double, double>> SP;
  typedef rsSyncPhasor<double, rsPolyBlep2<double, double>> SP;
  //typedef rsSyncPhasor<double, rsTableMinBlep<double, double>> SP;
  SP sp;
  sp.setMasterIncrement(f1      /fs);
  sp.setSlaveIncrement( f2_start/fs);
  sp.reset();


  // produce naive and anti-aliased signal:
  std::vector<double> xNaive(N), xBlep(N);
  for(int n = 0; n < N; n++) 
  {
    // sweep freq of slave osc:
    double f2_mix = n/(N-1.0);
    double f2     = (1-f2_mix) * f2_start + f2_mix * f2_end;
    sp.setSlaveIncrement(f2/fs);

    // produce naive and anti-aliased signals:
    xNaive[n] = sp.getSampleNaive();
    xBlep[n]  = sp.applyBlep(xNaive[n]);

    // convert range and apply amplitude scaler:
    xNaive[n] = 2 * xNaive[n] - 1;
    xBlep[n]  = 2 * xBlep[n]  - 1;
    xNaive[n] *= a;
    xBlep[n]  *= a;
  }

  rsArrayTools::shift(&xBlep[0], N, -sp.getBlep().getDelay());

  rsPlotVectors(xNaive, xBlep);
  rosic::writeToMonoWaveFile(  "SyncPhasorNaive.wav", &xNaive[0],            N, (int) fs);
  rosic::writeToMonoWaveFile(  "SyncPhasorBlep.wav",  &xBlep[0],             N, (int) fs);
  rosic::writeToStereoWaveFile("SyncPhasorBoth.wav",  &xNaive[0], &xBlep[0], N, (int) fs);

  // master-then-slave wraparound produces artifacts

  // todo: produce the same sweep naively at a much higher sample-rate and then downsample and 
  // compare with blep version - make a function 
  // createSyncOutput(masterFreq, slaveStartFreq, slaveEndFreq, length, sampleRate, bool withBlep)
  // -produce difference signals between naive and polyblep versions - Elan says, the PolyBlep1 
  //  works better for that - but why should that be the case? ...it shouldn't - i think. 
  //  -> figure out
}

template<class TOsc>
std::vector<double> createSyncOutputNaive(TOsc& osc, int numSamples)
{
  std::vector<double> y(numSamples);
  osc.reset();
  for(int n = 0; n < numSamples; n++)
    y[n] = osc.getSampleNaive();
  return y;
}
template<class TOsc>
std::vector<double> createSyncOutputBlep(TOsc& osc, int numSamples)
{
  std::vector<double> y(numSamples);
  osc.reset();
  for(int n = 0; n < numSamples; n++)
    y[n] = osc.getSample();
  return y;
}
template<class TOsc>
void plotSyncOutput(TOsc& osc, int numSamples)
{
  std::vector<double> yNaive = createSyncOutputNaive(osc, numSamples);
  std::vector<double> yBlep  = createSyncOutputBlep( osc, numSamples);
  RAPT::rsArrayTools::shift(&yBlep[0], numSamples, -osc.getBlep().getDelay());
  rsPlotVectors(yNaive, yBlep);
}
template<class TOsc>
void plotSyncOutput(TOsc& osc, int numSamples, double masterInc, double slaveInc)
{
  osc.setMasterIncrement(masterInc);
  osc.setSlaveIncrement( slaveInc);
  plotSyncOutput(osc, numSamples);
}
template<class TOsc>
void writeSyncOutput(TOsc& osc, int numSamples, double masterInc, double slaveInc, int fs)
{
  int N = numSamples;

  osc.setMasterIncrement(masterInc);
  osc.setSlaveIncrement( slaveInc);
  std::vector<double> yNaive = createSyncOutputNaive(osc, N);
  std::vector<double> yBlep  = createSyncOutputBlep( osc, N);
  RAPT::rsArrayTools::shift(&yBlep[0], N, -osc.getBlep().getDelay());

  yNaive = yNaive - 0.5;
  yBlep  = yBlep  - 0.5;

  rosic::writeToMonoWaveFile(  "SyncNaive.wav", &yNaive[0],            N, (int) fs);
  rosic::writeToMonoWaveFile(  "SyncBlep.wav",  &yBlep[0],             N, (int) fs);
  rosic::writeToStereoWaveFile("SyncBoth.wav",  &yNaive[0], &yBlep[0], N, (int) fs);
}
void syncPhasor2()
{
    // we use PolyBlep1 because it has a blep corrector whose desired state can be most easily 
  // predicted
  typedef rsSyncPhasor<double, rsPolyBlep1<double, double>> SP;
  //typedef rsSyncPhasor<double, rsPolyBlep2<double, double>> SP;
  SP sp;


  plotSyncOutput(sp, 3200, 31./1024, 128./1024);
  writeSyncOutput(sp, 320000, 31./1024, 128./1024, 44100);
  // alternately enters the "slave-then-master" branch 1x then "slave-only" branch 4x
  // hmm - actually here even the naive signal looks strange


  // masterInc 32/1024 results in a master freq of 1378.125 Hz @ 44.1

  //plotSyncOutput( sp, 500,   32./1024, 129./1024);
  //writeSyncOutput(sp, 50000, 32./1024, 129./1024, 44100);
  // alternately enters the "slave-then-master" branch 1x then "slave-only" branch 3x
  // that looks pretty good

  // slaveFreq = 4*masterFreq: master and slave resets always occur exactly simultaneously - we 
  // expect the master reset to be inconsequential, i.e. the code for master-reset is called and 
  // the blep is prepared, but always with zero step-amplitude:
  //plotSyncOutput( sp, 100,   32./1024, 128./1024);
  //writeSyncOutput(sp, 10000, 32./1024, 128./1024, 44100);
  // alternately enters the "slave-then-master" branch 1x then "slave-only" branch 3x

  // we expect the blep-signal to be equal to the naive signal (but one sample delayed) because the
  // slave resets occur exactly at integer sample instants and master resets are supposed to be 
  // inconsequential anyway

  // maybe we should manually walk through the samples that are supposed to be generated and define
  // the desired output signals and then check, if they are actually produced (don't check internal
  // states of the osc and its embedded blep directly - just look at the output signal)

  // hmm...that doesn't look like a pure delay
  // ahh - not we should expect a delay and some sort of lowpass character - that seems to fit with
  // what we see
  // the slave-wraparounds write 0.5 nonzero numbers into corrector of the blep (i think, the delayed
  // sample remains untouched). the master wraparounds are indeed inconsequential

  // maybe we should first unit-test ts rsPolyBlep1 object - why would it write 0.5 into the corrector 
  // when the wraparound occurs exactly at the sample-instant?
  // rsPolyBlep1::prepareForStep does this:
  //   corrector += a * (-d2/2 + d - 1./2);
  // this creates the .5 in the corrector - is this really correct? ..but the informal anti-aliasing 
  // tests with it look good - they are probably correct and what we see is the lowpass character
  // that affects the passband




  //plotSyncOutput( sp, 500,   32./1024, 127./1024);
  //writeSyncOutput(sp, 50000, 32./1024, 127./1024, 44100);
  // alternately enters the "master-only" branch 1x and "slave-only" branch 3x
  // that looks wrong! i think, it exposes the artifact that i'm battling with
  // ..or..well - the black (naive) looks wrong but the blue (blepped) looks actually nice


  // we need a case that sometimes enters the "master-then-slave" branch - i think, we need a 
  // master increment that causes master-wrap-arounds at non-integers
  // ...maybe 31./1024, 33./1024, 31./1023 or so

  // compare with a signal generated by oversampling a naively synthesized sync signal - maybe the
  // artifact is inherent to sync and not a bug?

  // for the implementation - check this:
  // -we need figure outthe normalized position of the slave at the sync event before reset/wrap 
  //  for computing the step amplitude (for slave-resets, this is always 1 because that's the 
  //  slave's wraparound point, but for master-resets, it may be anywhere)
}

void syncOsc()
{
  int N = 500000;
  double fs       = 44100;
  double f1       = 100 * GOLDEN_RATIO;  // master osc freq
  double f2_start = 20 * f1;             // start freq of slave osc 
  double f2_end   = 2  * f1;             // end   freq of slave osc 2
  double a        = 0.5;                 // amplitude


  //f1 = fs/100;
  //f2_start = f2_end = 4.1 * f1;

  //double pm = 100.25; // master period
  //double ps = 30.;     // slave period
  //f1 = fs/pm;
  //f2_start = f2_end = fs/ps;


  typedef rsSyncOsc<double, rsPolyBlep2<double, double>> SO;
  //typedef rsSyncOsc<double, rsTableMinBlep<double, double>> SO;
  SO osc;
  osc.setMasterIncrement(f1/fs);
  osc.reset();

  // produce naive and anti-aliased signal:
  std::vector<double> xNaive(N), xBlep(N);
  for(int n = 0; n < N; n++) 
  {
    // sweep freq of slave osc:
    double f2_mix = n/(N-1.0);
    double f2     = (1-f2_mix) * f2_start + f2_mix * f2_end;
    osc.setSlaveIncrement(f2/fs);

    // produce naive and anti-aliased signals:
    xNaive[n] = osc.getSampleNaive();
    xBlep[n]  = osc.applyBlep(xNaive[n]);
    xNaive[n] *= a;
    xBlep[n]  *= a;
  }


  rsArrayTools::shift(&xBlep[0], N, -osc.blep.getDelay());

  //rsPlotVectors(xNaive, xBlep);
  rosic::writeToMonoWaveFile(  "SyncNaive.wav", &xNaive[0],            N, (int) fs);
  rosic::writeToMonoWaveFile(  "SyncBlep.wav",  &xBlep[0],             N, (int) fs);
  rosic::writeToStereoWaveFile("SyncBoth.wav",  &xNaive[0], &xBlep[0], N, (int) fs);

  // when using slave.resetPhase(newPhase + slave.getPhaseIncrement()); the phase seems to be in
  // advance after a sync trigger, without the +inc, we get wrong step-amplitudes

  // hmm. i think, +inc is actually correct

  // i need to devise a test scenario with clearly defined master and slave peridicities and
  // walk through it manually, to figure out what is actually desired
  // maybe use pMaster = 100 samples ( = 2^2 * 5^2), pSlave = 21 samples ( = 3*7)
  // ...waveform seems to be strictly periodic with master period - shouldn't each cycle look
  // slightly different with respect ot blep phase? ah - no - only if, the master-cycle length is
  // non-integer (try 100.25), 100.25/30.0 looks actually really good

  // with a slow sweep, it looks like 

}

void dualBlepOsc()
{
  int N = 800;
  double fs = 44100;
  //double f1 = 200;
  double f1 = 100 * GOLDEN_RATIO;    // freq of osc 1
  //double f2 = f1 * GOLDEN_RATIO;
  double f2_start = 20 * f1;         // start freq of osc 2
  double f2_end   = 2  * f1;         // end   freq of osc 2
  double a        = 1.0;             // amplitude


  // test:
  //f2_start = f2_end = 5000;

  f1 = fs/100;
  f2_start = f2_end = 4.99 * f1;





  //f1 = 100;
  //f2 = 420;
  //f2 = 4.25 * f1;

  // maybe make a sweep from 5 kHz to 500 Hz for osc2





  // create and set up dual osc:
  //typedef rsDualBlepOsc<double, rsPolyBlep2<double, double>> DBO;
  typedef rsDualBlepOsc<double, rsTableMinBlep<double, double>> DBO;
  DBO dualOsc;
  dualOsc.setPhaseIncrement1(f1/fs);
  //dualOsc.setPhaseIncrement2(f2/fs);

  // try sync with a slave increment of 0.1 and master increment of 0.019..0.021:
  //dualOsc.setPhaseIncrement1(0.015);  // osc1 is master
  //dualOsc.setPhaseIncrement2(0.1);   // osc2 is slave
  dualOsc.setSync12(true);

  // generate both output channels and their sum:
  dualOsc.reset();
  std::vector<double> x1n(N), x2n(N), x1(N), x2(N), x(N);
  for(int n = 0; n < N; n++) 
  {

    // sweep freq of osc2:
    double f2_mix = n/(N-1.0);
    double f2     = (1-f2_mix) * f2_start + f2_mix * f2_end;
    dualOsc.setPhaseIncrement2(f2/fs);


    dualOsc.getSamplePairNaive(&x1n[n], &x2n[n]);
    dualOsc.applyBleps(&x1n[n], &x2n[n], &x1[n], &x2[n]);

    x1n[n] *= a;
    x2n[n] *= a;
    x1[n]  *= a;
    x2[n]  *= a;
    x[n] = x1[n] + x2[n];
  }

  // plot and/or write wavefile:
  //rsPlotVectors(x1, x2, x);
  //rsPlotVectors(x1, x2);
  //rsPlotVector(x2);
  rsPlotVectors(x2n, x2);
  //rosic::writeToStereoWaveFile("DualBlepOsc.wav",  &x1[0], &x2[0], N, (int) fs);
  //rosic::writeToMonoWaveFile(  "DualBlepOsc2.wav", &x2[0],         N, (int) fs);

  // seems like the osc2 responds to the sync triggers with 1 sample delay - something is wrong

  // somehow it seems, the step-amplitude is wrong for the sync-steps

  // maybe plot naive and anti-aliased versions side by side
  // -ok - done - this is strange: the sync-steps somehwo do not hsow the same blep-delay latency 
  //  as the regular swatooth steps - also they overshoot more than they should
  // -with f2 = 4.9 * f1, it's very visible...or even 4.99 * f1
  //  f1 = fs/100 works well for the plots
  // -ok - using new-old instead of old-new for the step-amp fixes the delay - this seems to be the
  //  correct thing to do also theoretically - but now there are additional spikes - is this
  //  because the blep is added/prepared twice?
  // -with f2 = 5.1 * f1, it looks fine
  // -with f2 = 3.99 * f1, the same thing occurs, but now we have spikes at every 4th cycle instead
  //  of every 5th
  // -the amplitude of the additional spikes is 2*a, where a is the waveform amplitude - so they 
  //  spike up to a + 2*a = 3*a and it seems it is really only a single sample that is wrong
  //  -interestingly, with the polyblep, the spike is only at 2a, not at 3a
  // -maybe try using a separate blep object for the sync...but i really see nor reason why we 
  //  couldn't use the same blep for both sorts of steps
  // -moving *x2 = osc2.getSampleSaw(); below the sync-handling for osc2 removes the spike but 
  //  introduces and additional delay of 1 sample
  // -so, it seems first, osc2 detects a step and then due to sync another step is detected
  //  ...maybe we should reset the stepAmplitude member in osc2 when we do a step due to to sync?
  // -adding  oldPhase += osc2.getPhaseIncrement();  osc2.wrapPhase(oldPhase); seeems to fix the
  //  really bad spike - but then the regular delay is missing for the sync-steps

  // -maybe it's better to not treat the 2 oscs on the same footing - use a simple phasor as
  //  master and let only the slave be a full osc - then we only need a single blep object

  int dummy = 0;
}

void bouncillator()
{
  rsBouncillator<float> rb;
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

  rsBouncillator<float> bnc;
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
  rsRayBouncer<float> rb;
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

  for(int n = 0; n < N; n++)
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

void snowFlake()
{
  // We test the rendering and anti-aliasing of snowflake - the anti aliasing shall be done by 
  // means of bleps (for resets) and blamps (for turns)

  int N    = 44100;  // number of samples
  int fs   = 44100;  // sample rate
  double f = 1000;   // signal frequency
  double resetRatio = 0.251;  // use 0 to turn resetting off

  //resetRatio  = 0.0;
  f *= sqrt(2.0);  // test
  //f = 1102.5;

  //fs = 32; f = 1.1;
  //fs = 32; f = 1.01;
  //fs = 32; f = 1.001;  // fails - adds upward spikes
  //fs = 32; f = 1.0001;  // fails - adds upward spikes
  //fs = 32; f = 1;     // works with edge-case check
  //fs = 32; f = 0.999;     // fails - adds downward spikes
  //fs = 32; f = 0.99;
  //fs = 32; f = 0.9;

  using Vec = std::vector<double>;

  rosic::Snowflake sf;
  sf.setSampleRate(fs);
  sf.setFrequency(f);
  sf.setAmplitude(0.5);
  //sf.setPhaseOffset(0.0625 * 360);
  //sf.setPhaseOffset(0.125 * 360);
  //sf.setPhaseOffset(0.6 * 360);
  sf.setAxiom("F+F+F+F+"); sf.setTurnAngle(90);    // square
  //sf.setAxiom("F+F+F+F+F+"); sf.setTurnAngle(72);    // pentagon
  //sf.setAxiom("F+F+F+F+F+F+"); sf.setTurnAngle(60);    // hexagon
  sf.setUseTable(false);
  sf.setAntiAlias(false);
  sf.setNumIterations(0);
  sf.setResetRatio(0, resetRatio);

  // produce signal without anti-aliasing:
  Vec xL1(N), xR1(N);
  for(int n = 0; n < N; n++)
    sf.getSampleFrameStereo(&xL1[n], &xR1[n]);

  // produce signal with anti-aliasing:
  Vec xL2(N), xR2(N);
  sf.setAntiAlias(true);
  sf.reset();
  for(int n = 0; n < N; n++)
    sf.getSampleFrameStereo(&xL2[n], &xR2[n]);

  // plot:
  int numPlotPoints = 1000;
  Vec dL1 = rsDifference(xL1);  // difference
  Vec dL2 = rsDifference(xL2);
  Vec cL1 = rsDifference(dL1);  // 2nd difference (curvature)
  Vec cL2 = rsDifference(dL2);
  //rsPlotArrays(numPlotPoints, &xL1[0], &xL2[0]);     // left without and with anti-alias
  //rsPlotArrays(numPlotPoints, &dL1[0], &dL2[0]);   // left difference without and with anti-alias
  //rsPlotArrays(numPlotPoints, &cL1[0], &cL2[0]);
  //rsPlotArrays(numPlotPoints, &xR1[0], &xR2[0]);
  //rsPlotArrays(numPlotPoints, &xL1[0], &xR1[0]);   // left/right non anti-aliased
  rsPlotArrays(numPlotPoints, &xL2[0], &xR2[0]);   // left/right anti-aliased


  // write files:
  rosic::writeToStereoWaveFile("SnowflakeTest.wav",   &xL1[0], &xR1[0], N, fs);
  rosic::writeToStereoWaveFile("SnowflakeTestAA.wav", &xL2[0], &xR2[0], N, fs);


  // Observations:
  // -Tests passed with 1 resetter with resetRatios of: 0.49, 0.51, 1.5, 1.9, 2.1
  // -Tests failed with 1 resetter with resetRatios of: 0.25, 0.50, 1.0, 2.0
  //  -it seems only the right channel is wrong (that's strange!)
  //  -it has probably to do with the reset occuring at a moment exactly when a cycle finishes
  //  -with 0.25, we see 3 correct cycles alternate with 1 messed up cycle
  //  -we actually do not fall into the "if(iPos != lineIndex)" branch in cases where a reset 
  //   occured before for this sample
  //  -apparently, stepX is always 0 and stepY is always 1 when the reset occurs - but that is 
  //   wrong! both steps should be 0, right?
  //  -using "F+F+F+F" as seed withut the final + does not seem to help
  //   with f = 1102.5 there 8 initial correct cycles before it goes wrong
  //  -when using an inc = 1/32 (via fs=32, f=1), the glitches occur at k*128 samples
  // -a resetRatio of 0.252 has occasional spikes - it seems like the closer the ratio is to a 
  //  fraction of the cycle-length, the denser the spurious spikes get, so maybe we should assume
  //  they are always there, regardless of the ratio, just more or less frequently
  // -probably they always happen when two events (reset and line-switch?) occur within one sample?
  //  ...but the line-switch does not happen anymore, after a reset occured...maybe it should? but 
  //  no, that doesn't seem to make sense - a reset actually switches to a new line segment, too
  //  but maybe a line-switch may happen a sample *before* the reset and then we have two 
  //  correction bleps...or at least 2 blamps
  // -oh: it does happen to x, too when we use, for example setPhaseOffset(0.125 * 360)
  // -try to figure out what conditions are different when the spikes occur (fPos, targetPhase, 
  //  etc.)

  // -i think, TurtleSource::resetPhase must take into account the starPhase set by 
  //  setPhaseOffset(0.6 * 360); otherwise we trigger an assert right in first getSample call 
  //  because we get an invalid blepTi - we should make sure that the (iPos != lineIndex) branch is
  //  not entered in the first call...but that could give rise to more porblems when we jump over
  //  several line segments in one sample - in thins case, the blepTime = fPos / (numLines*inc) may 
  //  not work

  // -the inc depends on the resetRatio - why? probably to tune the funcdamental frequency 
  //  independently from the reset ratio
  // -the first resetter is by default not off but instead set to a reset-ratio of 1, turning 
  //  resetting off actually produces a clean signal even without anti-aliasing

  // Notes:
  // -let N be the number of line-segments, then: pos = [0..1), linePos = [0..N), 
  //  iPos = floor(linePos) = [0..N-1] (is N-1 really included - why?), 
  //  fPos = linePos - iPos = [0..1), fPos is used to interpolate between start and end of current 
  //  line segment  
  // -to be free from alisiasng when using blamps for the line segment switches, we need 
  //  inc <= 1/N (i think)
  // -if a blamp occurs due to going from line i to i+1 during nomral(!) playback, the formula
  //  blepTime = fPos/(N*inc) works (why?), but in the very first sample or when we skip segments 
  //  due to high inc (violating inc <= 1/N), this formula seems inappropriate...we may get 
  //  blepTime > 1 ...maybe we should fmod the value by 1? does that make sense? or do we need an 
  //  entirely different approach?
  // -actually, the line-segment change could also create a step-discontinuity, but in normal 
  //  operation (where pos increases by inc (<= 1/N) each sample), this only is the case when we
  //  wrap around, i.e. i = N-1, i+1 = wrapAround(N) = 0 and the turtle's path is not closed
  // -maybe it would be more convenient to have a turtle-source that allows for random access 
  //  values: double Turtle::getSampleAt(double phase) ...it could internally compute the increment
  //  by inc = phase - oldPhase, if phase and oldPhase are both in [0..1), inc is in (-1..+1) for 
  //  each sample, oldPhase would be member variable
  //  -this may facilitate phase-modulation more easily than the current implementation
  //  -we could build a driver around it that creates a swatooth shape for the phase


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
  RAPT::rsArrayTools::fillWithRangeExponential(&f[0], N, fL, fU);

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

void shepardTone()
{
  using Real = double;

  Real sampleRate = 44100.f;
  Real length     = 3.0;     // length of 1 cycle in seconds
  int  numCycles  = 4;       // number of cycles to generate


  Real speed = 12.0 / length;

  rsShepardToneGenerator<Real> stg;
  stg.setSampleRate(sampleRate);
  stg.setBellFloor(0.00001);
  stg.setSpeed(speed);
  stg.updateCoeffs();

  using AT  = RAPT::rsArrayTools;
  using Vec = std::vector<Real>;

  int N = (int) round(length * sampleRate * numCycles);
  Vec y(N);
  for(int n = 0; n < N; n++)
    y[n] = stg.getSample();

  AT::normalize(&y[0], N, Real(0.5));
  rosic::writeToMonoWaveFile("ShepardTone.wav", &y[0], N, (int)sampleRate);


  /*
  // Generate and plot the frequency response of the pseudo "filter"
  int numBins = 1000;
  Vec freqs(numBins), amps(numBins);
  AT::fillWithRangeExponential(&freqs[0], numBins, 20.0, 20000.0);
  for(int i = 0; i < numBins; i++)
    amps[i] = stg.getGainForFrequency(freqs[i]);
  GNUPlotter plt;
  plt.addDataArrays(numBins, &freqs[0], &amps[0]);
  plt.setLogScale("x");
  plt.plot();
  */



  // ToDo: 
  // -there are still discontinuities at the wrap-arounds of the time variable. i think, they 
  //  happen because the sine that is switched off is not yet totally zero and the one that is 
  //  switched on does not start at zero?
  //  -figure out by generating only bottom and top sine and/or 
  //  -using a lower floor does not seem to help
  //  -commenting out the if(p < centerPitch - halfWidth) ... stuff also does not help
  //  -maybe inconsistent use of > and >= in the various wraparounds? nope: using >= consistently
  //   also doesn't fix is - but is probably a good idea anyway
  //  -it seems like the phase jumps suddenly - maybe it's not a good idea to multiply the phase
  //   itself by the factor?
  // -what about using using an actual filter instead of the pseduo-filter?

  int dummy = 0;
}