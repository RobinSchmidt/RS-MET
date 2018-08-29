template <class T>
rsFilterSpecificationZPK<T>::rsFilterSpecificationZPK(const std::vector<std::complex<T>>& poles,
  const std::vector<std::complex<T>>& zeros, T gain, T sampleRate)
{
  this->poles = poles;
  this->zeros = zeros;
  this->gain  = gain;
  this->sampleRate = sampleRate;
}

//=================================================================================================

template <class T>
rsFilterSpecificationBA<T>::rsFilterSpecificationBA(
  const std::vector<T>& num, const std::vector<T>& den, T sampleRate)
{
  this->sampleRate = sampleRate;
  b = num;
  a = den;
}
