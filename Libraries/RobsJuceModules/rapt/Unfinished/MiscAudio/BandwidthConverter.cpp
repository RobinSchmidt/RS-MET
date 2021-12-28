template<class T>
T rsBandwidthConverter::bandedgesToCenterFrequency(T fl, T fu)
{
  return sqrt(fl*fu);
}

template<class T>
T rsBandwidthConverter::bandedgesToAbsoluteBandwidth(T fl, T fu)
{
  return fu - fl;
}

template<class T>
T rsBandwidthConverter::relativeBandwidthToBandedgeFactor(T br)
{
  return 0.5*br + sqrt(0.25*br*br + 1);
}

template<class T>
void rsBandwidthConverter::absoluteBandwidthToBandedges(T bw, T fc, T *fl, 
  T *fu)
{
  T br = bandwidthAbsoluteToRelative(bw, fc);
  T k  = relativeBandwidthToBandedgeFactor(br);
  *fl = fc/k;
  *fu = fc*k;
}

template<class T>
T rsBandwidthConverter::bandwidthAbsoluteToRelative(T bw, T fc)
{
  return bw / fc;
}

template<class T>
T rsBandwidthConverter::absoluteBandwidthToQ(T bw, T fc)
{
  return 1.0 / bandwidthAbsoluteToRelative(bw, fc);
}

template<class T>
T rsBandwidthConverter::bandedgesToBandwidthInOctaves(T fl, T fu)
{
  return rsLog2(fu/fl);
}

template<class T>
T rsBandwidthConverter::absoluteBandwidthToOctaves(T bw, T fc)
{
  T fl, fu;
  absoluteBandwidthToBandedges(bw, fc, &fl, &fu);
  return bandedgesToBandwidthInOctaves(fl, fu);
}

template<class T>
T rsBandwidthConverter::multipassScalerButterworth(int M, int N, T g)
{
  return pow(pow(g, -2.0/M)-1, -0.5/N);

  // The formula was derived by considering the analog magnitude-squared response of a N-th order
  // Butterworth lowpass applied M times, given by: g^2(w) = (1 / (1 + (w^2/w0^2)^N))^M. To find 
  // the frequency w, where the magnitude-squared response has some particular value p, we solve 
  // for w. Chosing some particluar value for g (such as g = sqrt(0.5) for the -3.01 dB 
  // "half-power" point), we can compute the scaler s, at which we see this value of p. Scaling our 
  // cutoff frequency w0 by 1/s, we shift the half-power point to w0 in the multipass case. 

  // \todo to apply this to bilinear-transform based digital filters, we should replace w by 
  // wd = tan(...) and solve....
}

// for energy normalization, use the total energy formula (that i have obtained via sage)
// E = pi*gamma(M - 1/2/N)/(N*gamma(M)*gamma(-1/2/N + 1)*sin(1/2*pi/N)) // 1/2 -> 0.5 
// E = pi*gamma(M - 0.5/N)/(N*gamma(M)*gamma(-0.5/N + 1)*sin(0.5*pi/N)) // k = 0.5/N
// k = 0.5/N
// E = pi*gamma(M-k) / (N*gamma(M)*gamma(1-k)*sin(k*pi)) 
// maybe compare the energy normalization to the formula above
// this formula has been implemented in
// rsPrototypeDesigner<T>::butterworthEnergy

template<class T>
T rsBandwidthConverter::lowpassResoGainToQ(T a)
{
  a *= a;
  T b = sqrt(a*a - a);
  T P = T(0.5) * (a + b);  // P := Q^2
  return sqrt(P);

  // The formula was derived by considering the magnitude response of a 2nd order analog lowpass 
  // with transfer function H(s) = 1 / (1 + s/Q + s^2). The magnitude-squared response of such a 
  // filter is given by: M(w) = Q^2 / (Q^2 (w^2 - 1)^2 + w^2). This function has maxima of height
  // (4 Q^4)/(4 Q^2 - 1) at w = -1/2 sqrt(4 - 2/Q^2). The expression for the height was solved 
  // for Q. We get a quadratic equation for P := Q^2. The correct solution was picked empirically.
  // It was also observed empirically that the formula also works for highpass filters.
}

/*

More formulas that need some numeric checks:

fu = fl * 2^bo
fc = fl * sqrt(2^bo) = fl * 2^(bo/2)
k  = 2^(bo/2)
Q  = 2^(bo/2) / (2^bo - 1)

*/