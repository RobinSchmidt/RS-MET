
void rsModalFilterFloatSSE2::setParametersTwoEnvs(
  double w, double A, double p, 
  double att1, double att2, double attB,
  double dec1, double dec2, double decB)
{
  //double A = E;  // preliminary
  // todo: compute energy, optionally divide A by sqrt(energy), optimize computations (some 
  // values are computed 4 times (cw,sw,cp,sp) in the 4 calls below - but do this in production 
  // code):

  //auto calc = &RAPT::rsDampedSineFilterCoeffs<float, float>;
  //calc(w, (1-attB)*A, att1, p, &b0[0], &b1[0], &a1[0], &a2[0]);
  //calc(w,    attB *A, att2, p, &b0[1], &b1[1], &a1[1], &a2[1]);
  //calc(w, (1-decB)*A, dec1, p, &b0[2], &b1[2], &a1[2], &a2[2]);
  //calc(w,    decB *A, dec2, p, &b0[3], &b1[3], &a1[3], &a2[3]);

  RAPT::rsDampedSineFilterCoeffs(w, (1-attB)*A, att1, p, &b0[0], &b1[0], &a1[0], &a2[0]);
  RAPT::rsDampedSineFilterCoeffs(w,    attB *A, att2, p, &b0[1], &b1[1], &a1[1], &a2[1]);
  RAPT::rsDampedSineFilterCoeffs(w, (1-decB)*A, dec1, p, &b0[2], &b1[2], &a1[2], &a2[2]);
  RAPT::rsDampedSineFilterCoeffs(w,    decB *A, dec2, p, &b0[3], &b1[3], &a1[3], &a2[3]);

  // Filters 0 and 1 need a minus sign to make the final getSum work (todo: explain why):
  b0[0] = -b0[0], b1[0] = -b1[0];
  b0[1] = -b0[1], b1[1] = -b1[1];
}

void rsModalFilterFloatSSE2::setParameters(double w, double A, double p, 
  double att, double dec, double dw, double dp, double b, double attScl, double decScl)
{
  //auto calc = &RAPT::rsDampedSineFilterCoeffs<float, float>;

  // 1st attack/decay filter pair:
  double c = 1-b;
  RAPT::rsDampedSineFilterCoeffs(w-0.5*dw, c*A, att/attScl, p-0.5*dp, &b0[0], &b1[0], &a1[0], &a2[0]);
  RAPT::rsDampedSineFilterCoeffs(w-0.5*dw, c*A, dec/decScl, p-0.5*dp, &b0[1], &b1[1], &a1[1], &a2[1]);

  // 2nd attack/decay filter pair:
  RAPT::rsDampedSineFilterCoeffs(w+0.5*dw, b*A, att*attScl, p+0.5*dp, &b0[2], &b1[2], &a1[2], &a2[2]);
  RAPT::rsDampedSineFilterCoeffs(w+0.5*dw, b*A, dec*decScl, p+0.5*dp, &b0[3], &b1[3], &a1[3], &a2[3]);

  // filters with indices 0 and 2 are for the attacks and need to be subtracted in final sum:
  b0[0] = -b0[0], b1[0] = -b1[0];
  b0[2] = -b0[2], b1[2] = -b1[2];

  // maybe dw, dp should also be scaled by b(lend) and not just by 0.5
}

