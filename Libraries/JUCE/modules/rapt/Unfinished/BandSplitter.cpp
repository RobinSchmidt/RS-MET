template<class TSig, class TPar>
void rsTwoBandSplitter<TSig, TPar>::setOmega(TPar newOmega)
{
  TPar c = tan(TPar(0.5) * newOmega);
  c = (c-1) / (c+1);
  b0 = b1 = TPar(0.5) * (1+c);
  a1 = c;  // use a1 directly in computation instead of temporary c
  // Formulas from DAFX, page 40.
}