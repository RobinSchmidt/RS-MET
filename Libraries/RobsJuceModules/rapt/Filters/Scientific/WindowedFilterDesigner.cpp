

template<class T>
void rsWindowedFilterDesigner::hilbert(T* h, int numTaps, rsWindowFunction::WindowType type)
{
  // Create the window:
  RAPT::rsWindowFunction::createWindow(h, numTaps, type, false);

  // Multiply in the Hilbert-filter weights:
  int m = numTaps/2;                          // Middle tap
  if(rsIsOdd(numTaps))
  {
    for(int k = m % 2; k < numTaps; k += 2)   // k starts at 0 if m is even and at 1 if m is odd
      h[k] = T(0);
    for(int k = 1; k <= m; k += 2)
    {
      T hk = T(2) / T(k*PI);
      h[m+k] *= +hk;
      h[m-k] *= -hk;
    }
  }
  else
  {
    for(int k = 0; k < m; k++)
    {
      T t  = T(k) + T(0.5);
      T hk = T(1) / (t*PI);
      h[m+k]   *= +hk;
      h[m-k-1] *= -hk;
    }
  }

  // See:
  // https://www.kvraudio.com/forum/viewtopic.php?t=608320
  // https://en.wikipedia.org/wiki/Hilbert_transform#Discrete_Hilbert_transform
  // https://www.dsprelated.com/freebooks/sasp/Hilbert_Transform_Design_Example.html
  // https://www.intechopen.com/chapters/39362
  //
  // ToDo:
  // -Compare the results of this routine with those of some reference implementations from octave 
  //  or numpy/scipy
}


template<class T>
void rsWindowedFilterDesigner::hilbertSmoothed(T* h, int numTaps, 
  rsWindowFunction::WindowType type, bool evenNominalLength)
{
  rsAssert(rsIsOdd(numTaps));

  using WFD = rsWindowedFilterDesigner;
  int M = numTaps;
  if(evenNominalLength)
  {
    WFD::hilbert(h, M-1, type);
    h[M-1] = T(0);
    rsArrayTools::movingAverage2ptForward(h, M, h);
  }
  else
  {
    WFD::hilbert(&h[1], M-2, type);
    h[0]   = 0;
    h[M-1] = 0;
    rsArrayTools::weightedAverage3pt(h, M, h, T(0.25), T(0.5), T(0.25));
  }

  // ToDo: 
  // -In the case of even nominal length, the main purposes of the MA smoothing is to bring the 
  //  delay from a half-integer to a full integer. Maybe that can be done in better ways by 
  //  interpolating the impulse response with polynomial interpolators. The 2-sample MA is 
  //  basically a linear interpolator that reads out the raw impulse response at half-integer 
  //  positions. We can do better than linear! Doing so may give more desirable magnitude 
  //  responses, i.e. less lowpassing. Although, that lowpassing might not be such a bad thing 
  //  in the context of instantaneous envelope detection.
}
