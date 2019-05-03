rsEllipticQuarterBandFilter::rsEllipticQuarterBandFilter()
{
  reset();  
}

void rsEllipticQuarterBandFilter::reset()
{
  for(int i = 0; i < 12; i++)
    w[i] = 0.0;
}

