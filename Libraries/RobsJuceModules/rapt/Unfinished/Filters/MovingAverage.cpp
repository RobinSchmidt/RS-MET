template<class TSig, class TPar>
rsMovingAverage<TSig, TPar>::rsMovingAverage()
{
  delayLine.setDelayInSamples(100);

  sampleRate = 44100.0;
  setLengthInSeconds(0.01);
  setDeviation(0.001);
  reset();
}

template<class TSig, class TPar>
void rsMovingAverage<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  sampleRate = newSampleRate;
  setLengthInSeconds(length);
}

template<class TSig, class TPar>
void rsMovingAverage<TSig, TPar>::setLengthInSeconds(TPar newLength)
{
  int N = rsRoundToInt(newLength * sampleRate);
  length = N / sampleRate;
  setLengthInSamples(N);
}

template<class TSig, class TPar>
void rsMovingAverage<TSig, TPar>::setLengthInSamples(int newLength)
{
  delayLine.setDelayInSamples(newLength);
  updateCoefficients();
}

template<class TSig, class TPar>
void rsMovingAverage<TSig, TPar>::setDeviation(TPar newDeviation)
{
  d = newDeviation;
  updateCoefficients();
}

template<class TSig, class TPar>
void rsMovingAverage<TSig, TPar>::reset()
{
  delayLine.reset();
  y1 = 0.0;
}

template<class TSig, class TPar>
void rsMovingAverage<TSig, TPar>::updateCoefficients()
{
  int  N = delayLine.getDelayInSamples();
  TPar k = (d+1)/N;     // d is the relative deviation
  TPar c = 1-k;

  // solve k*x^N - x + c via Newton iteration (x == a1):
  TPar x = 0.0;
  TPar f, fp;
  TPar dx = 1.0;
  TPar xN;
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
  TPar sn = 1.0;
  TPar an = a1;
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
