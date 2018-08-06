template<class T>
void rsTriSawOscillator<T>::updateTriSawCoeffs()
{
  // range variables - make settable members later:
  //T min = 0;
  //T max = 1;
  // hmm...but actually we need thme to be 0 and 1 because the waveshaping expects its input in
  // this range - if the range should really be different we may scale/shift as last step...

  //// coeffs:
  //a0 = min;
  //a1 = (max-a0) / h;
  //b0 = (max-h*min) / (1-h);
  //b1 = min-b0;

  // coeffs:
  a0 = 0;  // can actually be optimized away
  a1 = T(1) / h;
  b0 = T(1) / (1-h);
  b1 = -b0;

  // todo: we need to catch cases when h = 0 or 1-h = 0
  // if(h < eps)
  //   a1 = 0;
  // if(1-h < eps)
  //   b1 = 0;
  // or something
}