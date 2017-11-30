template<class TSig, class TPar>
void rsTwoBandSplitter<TSig, TPar>::setOmega(TPar newOmega)
{
  w = newOmega;
  TPar c = tan(TPar(0.5) * w);
  c = (c-1) / (c+1);
  b0 = b1 = TPar(0.5) * (1+c);
  a1 = c;  // use a1 directly in computation instead of temporary c
  // Formulas from DAFX, page 40.
}


template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::setSampleRate(TPar newSampleRate)
{

}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::setSplitFrequencies(const std::vector<TPar>& newFrequencies)
{

}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::setNumBands(int newNumBands)
{

}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::setSplitFrequency(int bandIndex, TPar newFrequency)
{

}

template<class TSig, class TPar>
void rsMultiBandSplitter<TSig, TPar>::updateSplitters()
{

}