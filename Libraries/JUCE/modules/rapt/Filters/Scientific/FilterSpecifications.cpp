template <class T>
rsFilterSpecificationZPK<T>::rsFilterSpecificationZPK(const std::vector<std::complex<T>>& poles,
  const std::vector<std::complex<T>>& zeros, T gain, T sampleRate)
{
  this->poles = poles;
  this->zeros = zeros;
  this->gain  = gain;
  this->sampleRate = sampleRate;
}

template <class T>
std::complex<T> rsFilterSpecificationZPK<T>::transferFunctionAt(std::complex<T> s_or_z)
{
  std::complex<T> num = 1, den = 1;    // numerator and denominator
  size_t i;
  if(sampleRate != RS_INF(double)) {   // digital
    std::complex<T> zr  = 1.0/s_or_z;  // z^-1
    for(i = 0; i < zeros.size(); i++) num *= (T(1) - zeros[i] * zr);
    for(i = 0; i < poles.size(); i++) den *= (T(1) - poles[i] * zr);
    // formula: https://ccrma.stanford.edu/~jos/filters/Factored_Form.html
    // maybe factor out into:
    // digitalTransferFunctionZPK(const Complex* zeros, int numZeros, 
    // const Complex* poles, int numPoles, Complex k, Complex z);
  }
  else {                               // analog (not yet tested)
    std::complex<T> s = s_or_z;
    for(i = 0; i < zeros.size(); i++) num *= (s - zeros[i]);
    for(i = 0; i < poles.size(); i++) den *= (s - poles[i]);
  }
  return k * num/den;
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
