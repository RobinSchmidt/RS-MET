template<class T>
void rsHeatEquation1D<T>::setMaxCycleLength(int newLength) 
{ 
  rodArray1.resize(newLength);
  rodArray2.resize(newLength);
  rsArray::fillWithZeros(&rodArray1[0], newLength);
  rsArray::fillWithZeros(&rodArray2[0], newLength);
  reset(); // possibly re-adjust pointers
}

template<class T>
void rsHeatEquation1D<T>::setHeatDistribution(T* d, int N)
{
  rsAssert( N <= (int) rodArray1.size() );
  rodLength = N;

  for(int i = 0; i < N; i++)
    rodIn[i] = d[i];
}

template<class T>
void rsHeatEquation1D<T>::setRandomHeatDistribution(int seed, int N)
{
  rsAssert( N <= (int) rodArray1.size() );
  rodLength = N;

  // todo: allow different probability densities for the prng



  rsNoiseGenerator<T> prng;
  prng.setSeed(seed);
  for(int i = 0; i < N; i++)
    rodIn[i] = prng.getSample();
}


template<class T>
void rsHeatEquation1D<T>::setTwoValueDistribution(T highFraction, int N)
{
  rsAssert( N <= (int) rodArray1.size() );
  rodLength = N;

  int nh = (int) round(highFraction*N);
  for(int i = 0; i < nh; i++) rodIn[i] = +1.0;
  for(int i = nh; i < N; i++) rodIn[i] = -1.0;
}


template<class T>
void rsHeatEquation1D<T>::normalizeHeatDistribution(T targetMean, T targetVariance)
{
  // set mean to desired target value (maybe factor out):
  //int N = (int) rodArray1.size();
  int N = (int) rodLength;
  typedef rsArray AR;
  T mean = AR::mean(rodIn, N);
  AR::add(rodIn, -mean, rodIn, N); 
  // ...hmm..this actually just set the mean to zero...which is the most reasonable target value
  // nayway (in case of sound generation)...but what if the user wants some other mean..we'll see

  // set variance to desired target variance:
  // ....
}