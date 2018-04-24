using namespace RSLib;

double rsBandwidthConverter::bandedgesToCenterFrequency(double fl, double fu)
{
  return sqrt(fl*fu);
}

double rsBandwidthConverter::bandedgesToAbsoluteBandwidth(double fl, double fu)
{
  return fu - fl;
}

double rsBandwidthConverter::relativeBandwidthToBandedgeFactor(double br)
{
  return 0.5*br + sqrt(0.25*br*br + 1);
}

void rsBandwidthConverter::absoluteBandwidthToBandedges(double bw, double fc, double *fl, 
  double *fu)
{
  double br = bandwidthAbsoluteToRelative(bw, fc);
  double k  = relativeBandwidthToBandedgeFactor(br);
  *fl = fc/k;
  *fu = fc*k;
}

double rsBandwidthConverter::bandwidthAbsoluteToRelative(double bw, double fc)
{
  return bw / fc;
}

double rsBandwidthConverter::absoluteBandwidthToQ(double bw, double fc)
{
  return 1.0 / bandwidthAbsoluteToRelative(bw, fc);
}

double rsBandwidthConverter::bandedgesToBandwidthInOctaves(double fl, double fu)
{
  return rsLog2(fu/fl);
}

double rsBandwidthConverter::absoluteBandwidthToOctaves(double bw, double fc)
{
  double fl, fu;
  absoluteBandwidthToBandedges(bw, fc, &fl, &fu);
  return bandedgesToBandwidthInOctaves(fl, fu);
}

double rsBandwidthConverter::multipassScalerButterworth(int M, int N, double g)
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
