template <class T>
std::complex<T> digitalTransferFunctionZPK(const std::complex<T>* zeros, size_t numZeros, 
  const std::complex<T>* poles, size_t numPoles, std::complex<T> k, std::complex<T> z)
{
  std::complex<T> zr  = T(1)/z;      // z^-1, reciprocal of z
  std::complex<T> num = 1, den = 1;  // numerator and denominator
  for(size_t i = 0; i < numZeros; i++) num *= (T(1) - zeros[i] * zr);
  for(size_t i = 0; i < numPoles; i++) den *= (T(1) - poles[i] * zr);
  return k * num/den;
  // formula: https://ccrma.stanford.edu/~jos/filters/Factored_Form.html
} // maybe move to rsFilterAnalyzer

template <class T>
std::complex<T> analogTransferFunctionZPK(const std::complex<T>* zeros, size_t numZeros, 
  const std::complex<T>* poles, size_t numPoles, std::complex<T> k, std::complex<T> s)
{
  std::complex<T> num = 1, den = 1;  // numerator and denominator
  for(size_t i = 0; i < numZeros; i++) num *= (s - zeros[i]);
  for(size_t i = 0; i < numPoles; i++) den *= (s - poles[i]);
  return k * num/den;
} // maybe move to rsFilterAnalyzer

template <class T>
std::complex<T> digitalTransferFunctionBA(const std::complex<T>* b, size_t Nb, 
  const std::complex<T>* a, size_t Na, std::complex<T> z)
{
  std::complex<T> num = 0, den = 0;
  for(size_t i = 0; i < Nb; i++) num += b[i] * pow(z, -T(i)); // can be optimized
  for(size_t i = 0; i < Na; i++) den += a[i] * pow(z, -T(i));
  return num/den;
} // maybe move to rsFilterAnalyzer

template <class T>
std::complex<T> analogTransferFunctionBA(const std::complex<T>* b, size_t Nb, 
  const std::complex<T>* a, size_t Na, std::complex<T> s)
{
  std::complex<T> num = 0, den = 0;
  for(size_t i = 0; i < Nb; i++) num += b[i] * pow(s, T(i)); // can be optimized
  for(size_t i = 0; i < Na; i++) den += a[i] * pow(s, T(i));
  return num/den;
} // maybe move to rsFilterAnalyzer

template <class T>
std::complex<T> dcGainNormalizer(const std::complex<T>* zeros, size_t numZeros,
  const std::complex<T>* poles, size_t numPoles)
{
  std::complex<T> H1 = digitalTransferFunctionZPK(zeros, numZeros, poles, numPoles, 
    std::complex<T>(1.0, 0.0), std::complex<T>(1.0, 0.0)); // H(z) at z=1
  return T(1) / H1;
}

//=================================================================================================

template <class T>
rsFilterSpecificationZPK<T>::rsFilterSpecificationZPK(
  const std::vector<std::complex<T>>& zeros, 
  const std::vector<std::complex<T>>& poles,
  std::complex<T> gain, T sampleRate)
{
  this->z = zeros;
  this->p = poles;
  this->k = gain;
  this->sampleRate = sampleRate;
  sortPolesAndZeros();
}

template <class T>
std::complex<T> rsFilterSpecificationZPK<T>::transferFunctionAt(std::complex<T> s_or_z) const
{
  if(isDigital())
    return digitalTransferFunctionZPK(&z[0], z.size(), &p[0], p.size(), k, s_or_z);
  else
    return analogTransferFunctionZPK( &z[0], z.size(), &p[0], p.size(), k, s_or_z);
}

template <class T>
rsFilterSpecificationBA<T> rsFilterSpecificationZPK<T>::toBA() const
{
  rsFilterSpecificationBA<T> ba;
  ba.sampleRate = sampleRate;
  ba.a = rsPolynomial<T>::rootsToCoeffs(p);
  ba.b = rsPolynomial<T>::rootsToCoeffs(z);
  if(isDigital()){
    rsReverse(ba.a);
    rsReverse(ba.b);
  }
  for(size_t i = 0; i < ba.b.size(); i++)
    ba.b[i] *= k;
  ba.normalizeDenominator(); // maybe make the normalization optional (on by default)
  return ba;
}

template <class T>
void rsFilterSpecificationZPK<T>::sortPolesAndZeros()
{
  if(p.size() > 0) rsHeapSort(&p[0], (int)p.size(), rsComplexLessByImRe<T>);
  if(z.size() > 0) rsHeapSort(&z[0], (int)z.size(), rsComplexLessByImRe<T>);
}

template <class T>
bool rsFilterSpecificationZPK<T>::equals(const rsFilterSpecificationZPK& zpk, T tol)
{
  bool r = sampleRate == zpk.sampleRate; // sample-rate uses exact comparison - it's not supposed to be subject to roundoff error
  r &= abs(k-zpk.k) <= tol;
  r &= rsAreVectorsEqual(z, zpk.z, tol);
  r &= rsAreVectorsEqual(p, zpk.p, tol);
  return r;
}

//=================================================================================================

template <class T>
rsFilterSpecificationBA<T>::rsFilterSpecificationBA(
  const std::vector<std::complex<T>>& num, const std::vector<std::complex<T>>& den, T sampleRate)
{
  this->sampleRate = sampleRate;
  b = num;
  a = den;
}

template <class T>
std::complex<T> rsFilterSpecificationBA<T>::transferFunctionAt(std::complex<T> s_or_z) const
{
  if(isDigital())
    return digitalTransferFunctionBA(&b[0], b.size(), &a[0], a.size(), s_or_z);
  else
    return analogTransferFunctionBA( &b[0], b.size(), &a[0], a.size(), s_or_z);
}

template <class T>
rsFilterSpecificationZPK<T> rsFilterSpecificationBA<T>::toZPK() const
{
  rsFilterSpecificationZPK<T> zpk;
  zpk.sampleRate = sampleRate;
  zpk.p.resize(a.size()-1);
  zpk.z.resize(b.size()-1);
  if(isDigital()) {
    std::vector<std::complex<T>> tmp = a;
    rsReverse(tmp);
    if(zpk.p.size() > 0) rsPolynomial<T>::roots(&tmp[0], (int)tmp.size()-1, &zpk.p[0]);
    tmp = b; 
    rsReverse(tmp);
    if(zpk.z.size() > 0) rsPolynomial<T>::roots(&tmp[0], (int)tmp.size()-1, &zpk.z[0]);
    zpk.k = b[0] / a[0];
  }
  else {
    if(zpk.p.size() > 0) rsPolynomial<T>::roots(&a[0], (int)a.size()-1, &zpk.p[0]);
    if(zpk.z.size() > 0) rsPolynomial<T>::roots(&b[0], (int)b.size()-1, &zpk.z[0]);
    zpk.k = b[b.size()-1] / a[a.size()-1];
  }
  zpk.sortPolesAndZeros();
  return zpk;
}

template <class T>
void rsFilterSpecificationBA<T>::normalizeDenominator()
{
  std::complex<T> s;
  if(isDigital()) s = T(1) / a[0];
  else            s = T(1) / a[a.size()-1];

  for(size_t i = 0; i < a.size(); i++) a[i] *= s;
  for(size_t i = 0; i < b.size(); i++) b[i] *= s;
}
