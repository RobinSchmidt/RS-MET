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
