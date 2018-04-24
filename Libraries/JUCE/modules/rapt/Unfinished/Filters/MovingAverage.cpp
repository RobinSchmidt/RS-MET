using namespace RSLib;

rsMovingAverage::rsMovingAverage()
{
  delayLine.setDelayInSamples(100);

  sampleRate = 44100.0;
  setLengthInSeconds(0.01);
  setDeviation(0.001);
  reset();
}

void rsMovingAverage::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  setLengthInSeconds(length);
}

void rsMovingAverage::setLengthInSeconds(double newLength)
{
  int N = rsRoundToInt(newLength * sampleRate);
  length = N / sampleRate;
  setLengthInSamples(N);
}

void rsMovingAverage::setLengthInSamples(int newLength)
{
  delayLine.setDelayInSamples(newLength);
  updateCoefficients();
}

void rsMovingAverage::setDeviation(double newDeviation)
{
  d = newDeviation;
  updateCoefficients();
}

void rsMovingAverage::reset()
{
  delayLine.reset();
  y1 = 0.0;
}

void rsMovingAverage::updateCoefficients()
{
  int    N = delayLine.getDelayInSamples();
  double k = (d+1)/N;     // d is the relative deviation
  double c = 1-k;

  // solve k*x^N - x + c via Newton iteration (x == a1):
  double x = 0.0;
  double f, fp;
  double dx = 1.0;
  double xN;
  while( fabs(dx) > 0.0 )
  {
    xN  = pow(x, N-1);
    fp  = N*k*xN - 1;
    xN *= x;
    f   = k*xN - x + c;
    dx  = f/fp;
    x  -= dx;
  }
  a1 = x;
  bN = pow(a1, N); 

  // Remark:
  // I'm not sure anymore, how exactly i defined the deviation "d" which led to the above equation,
  // but i think, defining it as d = 1 - a1^N = 1 - bN leads to a simpler calculation without
  // need for Newton iteration, namely: bN = 1-d, a1 = bN^(1/N). Maybe, i should switch to that. 
  // Defining the error d this way, seems to make a lot of sense anyway - it's the normalized 
  // difference between first and last value in the impulse response ("normalized" in the sense 
  // that it's divided by the first value, which is 1, before applying the overall gain).

  // calculate gain as reciprocal of the sum of the impulse-response values:
  double sn = 1.0;
  double an = a1;
  for(int n = 1; n < N; n++)
  {
    sn += an;
    an *= a1;
  }
  g = 1 / sn;
  
  // below is the analytical solution for the gain calculation as comment - when the leakage is 
  // close to zero, it tends to get less accurate than the loop-based calculation above (the sum 
  // of the impulse-response is not as close to unity as with the loop). This is because 
  // differences between very close numbers are involved which leads to precision loss.
  // ...can it be reformulated to avoid the precision loss?
  // if( 1-bN < DBL_EPSILON )  // would lead to division by very small number
  //   g = 1.0 / delayLine.getDelayInSamples();
  // else 
  //   g = (a1-1) / (bN-1);

}
