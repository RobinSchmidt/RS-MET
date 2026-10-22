
/*
template<class TSig, class TPar>
void rsHilbertFilter<TSig, TPar>::computeCoeffs(TPar* h, int numTaps, 
  RAPT::rsWindowFunction::WindowType type)
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

  // - Compare the results of this routine with those of some reference implementations from octave
  //   or numpy/scipy
  //
  // - Explain why this is commented out. I think, it may be obsolete due to the implementation of
  //   updateCoeffs() in the .h file? See also: rsWindowedFilterDesigner. It has a member function 
  //   hilbert() which seems to implement this algorithm here. Verify if the algorithm there is 
  //   really exactly the same as the one here. If so, delete this code here (but maybe move over 
  //   the comments from here before).
}
*/

//=================================================================================================
/*

Ideas:

- Implement a stereoizer based on such a Hilbert filter in the following way: (1) Convert L/R 
  (left/right) inputs to M/S (mid/side). (2) Apply the Hilbert filter to the mid signal. (3) Mix 
  the so produced signal into the side signal (weighted by an "Amount" parameter). (4) Convert back
  to L/R. Variations: Maybe we should have "Amount" slider for the Hilbert output as well as for 
  the original side signal. The length of the Hilbert filter could be another user parameter. 
  Shorter ones could actually make sense because they would leave the bass mono which is often 
  desired in mixing.
  

*/