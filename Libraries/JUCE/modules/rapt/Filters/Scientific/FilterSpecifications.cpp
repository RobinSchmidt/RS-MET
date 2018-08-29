template <class T>
rsFilterSpecificationZPK<T>::rsFilterSpecificationZPK(const std::vector<std::complex<T>>& poles,
  const std::vector<std::complex<T>>& zeros, T gain, T sampleRate)
{
  this->p = poles;
  this->z = zeros;
  this->k = gain;
  this->sampleRate = sampleRate;
}

template <class T>
std::complex<T> digitalTransferFunctionZPK(const std::complex<T>* zeros, size_t numZeros, 
  const std::complex<T>* poles, size_t numPoles, std::complex<T> k, std::complex<T> z)
{
  Complex zr  = 1.0/z;       // z^-1
  Complex num = 1, den = 1;  // numerator and denominator
  for(size_t i = 0; i < numZeros; i++) num *= (T(1) - zeros[i] * zr);
  for(size_t i = 0; i < numPoles; i++) den *= (T(1) - poles[i] * zr);
  return k * num/den;
  // formula: https://ccrma.stanford.edu/~jos/filters/Factored_Form.html
} // maybe move to rsFilterAnalyzer

template <class T>
std::complex<T> analogTransferFunctionZPK(const std::complex<T>* zeros, size_t numZeros, 
  const std::complex<T>* poles, size_t numPoles, std::complex<T> k, std::complex<T> s)
{
  Complex num = 1, den = 1;  // numerator and denominator
  for(size_t i = 0; i < numZeros; i++) num *= (s - zeros[i]);
  for(size_t i = 0; i < numPoles; i++) den *= (s - poles[i]);
  return k * num/den;
} // maybe move to rsFilterAnalyzer

template <class T>
std::complex<T> rsFilterSpecificationZPK<T>::transferFunctionAt(std::complex<T> s_or_z)
{
  if(sampleRate != RS_INF(double))
    return digitalTransferFunctionZPK(&z[0], z.size(), &p[0], p.size(), k, s_or_z);
  else
    return analogTransferFunctionZPK( &z[0], z.size(), &p[0], p.size(), k, s_or_z);
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
