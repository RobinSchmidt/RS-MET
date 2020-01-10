template<class T>
T polyEval(std::vector<T>& p, T x)
{
  int k = (int)p.size()-1;  // last valid index
  if(k < 0)
    return 0;
  T y = p[k];
  while(k > 0) {
    k -= 1;
    y = y*x + p[k]; }
  return y;
}

template<class T>
void polyTrunc(std::vector<T>& p, T tol)
{
  int i = (int)p.size();
  while(i > 1) {
  //while(i > 0) {  // old - this is a bug: a polynomial should have at least 1 coeff
    if(fabs(p[i-1]) > tol)
      break;
    i -= 1; }
  p.resize(i);
}

template<class T>
T makeMonic(std::vector<T>& p)
{
  T lc = rsLast(p);
  for(size_t i = 0; i < p.size(); i++)
    p[i] /= lc;
  return lc;
}

template<class T>
std::vector<T> polyAdd(
  const std::vector<T>& p, const std::vector<T>& q, 
  T tol, T wp, T wq)
{
  int np = (int) p.size();
  int nq = (int) q.size();
  std::vector<T> r;
  if(np >= nq) {
    r = wp*p;
    for(int i = 0; i < nq; i++) {
      r[i] += wq*q[i];
      polyTrunc(r, tol); }}
  else {
    r = wq*q;
    for(int i = 0; i < np; i++) {
      r[i] += wp*p[i];
      polyTrunc(r, tol); }}
  return r;
}

template<class T>
std::vector<T> polySub(const std::vector<T>& p, const std::vector<T>& q,
  T tol)
{
  return polyAdd(p, q, 1, -1, tol);
}

template<class T>
std::vector<T> polyMul(const std::vector<T>& x, const std::vector<T>& h,
  T tol)
{
  int L = (int)x.size() + (int)h.size() - 1;  // length of result
  std::vector<T> y(L);
  for(int n = 0; n < L; n++) {
    y[n] = 0;  
    for(int k = std::max(0, n-(int)x.size()+1); k < std::min(n+1, (int)h.size()); k++)
      y[n] += h[k] * x[n-k]; }
  polyTrunc(y, tol);
  return y;
}

template<class T>
void polyDivMod(std::vector<T> p, std::vector<T> d, 
  std::vector<T>& q, std::vector<T>& r, T tol)
{ 
  q.resize(p.size());
  r = p;                // init remainder with copy of product
  rsFill(q, 0.0);       // init quotient to all zeros
  int k = (int)p.size() - (int)d.size();
  while(k >= 0) {
    q[k] = r[(int)d.size()-1+k] / d[(int)d.size()-1];
    int j = (int)d.size()+k-2; 
    while(j >= k) {
      r[j] -= q[k] * d[j-k];
      j -= 1; }
    k -= 1; }
  for(int i = (int) d.size()-1; i < (int) p.size(); i++)
    r[i] = 0;
  polyTrunc(q, tol);
  polyTrunc(r, tol);
}

template<class T>
std::vector<T> polyDiv(std::vector<T> p, std::vector<T> d, T tol)
{
  std::vector<T> q, r;
  polyDivMod(p, d, q, r, tol);
  return q;
}

template<class T>
std::vector<T> polyMod(std::vector<T> p, std::vector<T> d, T tol)
{
  std::vector<T> q, r;
  polyDivMod(p, d, q, r, tol);
  return r;
}

template<class T>
bool isAllZeros(const std::vector<T>& v, T tol)
{
  for(size_t i = 0; i < v.size(); i++)
    if(fabs(v[i]) > tol)
      return false;
  return true;
}

template<class T>
std::vector<T> polyGCD(
  const std::vector<T>& p, const std::vector<T>& q, T tol, bool monic)
{
  std::vector<T> a = p, b = q, t;
  while(!isAllZeros(b, tol)) {
    t = b;
    b = polyMod(a, b, tol);
    a = t; }
  if(monic)
    makeMonic(a);
  return a;
}

template<class T>
std::vector<T> polyNest(const std::vector<T>& a, const std::vector<T>& b)
{
  int aN = (int)a.size()-1;               // degree of a
  int bN = (int)b.size()-1;               // degree of b
  int cN = aN*bN;                         // degree of result c
  std::vector<T> an(cN+1), c(cN+1); 
  rsFill(an, 0.0); an[0] = 1;             // powers of a, i.e. a^n - initially a^0 = [1 0 0...]
  rsFill(c,  0.0);                        // accumulator for result
  int K = 1;
  for(int n = 1; n <= bN; n++) {
    an = polyMul(an, a, 0.0);
    K += aN;
    for(int k = 0; k < K; k++)  
      c[k] += b[n] * an[k]; }
  return c;
}

template<class T>
void ratReduce(
  const std::vector<T>& pIn, const std::vector<T>& qIn,
  std::vector<T>& pOut, std::vector<T>& qOut, T tol)
{
  std::vector<T> gcd = polyGCD(pIn, qIn, tol);
  pOut = polyDiv(pIn, gcd, tol);
  qOut = polyDiv(qIn, gcd, tol);
}

template<class T>
void ratMul(
  const std::vector<T>& p, const std::vector<T>& q,
  const std::vector<T>& r, const std::vector<T>& s,
  std::vector<T>& u, std::vector<T>& v, T tol, bool reduced)
{
  u = polyMul(p, r, tol);
  v = polyMul(q, s, tol);
  if(reduced)
    ratReduce(u, v, u, v, tol);
}

template<class T>
void ratDiv(
  const std::vector<T>& p, const std::vector<T>& q,
  const std::vector<T>& r, const std::vector<T>& s,
  std::vector<T>& u, std::vector<T>& v, T tol, bool reduced)
{
  ratMul(p, q, s, r, u, v, tol, reduced); // r and s are swapped
}

template<class T>
void ratAdd(
  const std::vector<T>& n1, const std::vector<T>& d1,
  const std::vector<T>& n2, const std::vector<T>& d2,
  std::vector<T>& nr, std::vector<T>& dr, 
  T tol, T w1, T w2)
{
  std::vector<T> gcd, f1, f2, s1, s2;
  gcd = polyGCD(d1, d2, tol);
  f1 = polyDiv(d2, gcd, tol);
  f2 = polyDiv(d1, gcd, tol);
  dr = polyMul(f1, d1, tol);
  s1 = polyMul(f1, n1, tol);         // 1st summand in numerator of result
  s2 = polyMul(f2, n2, tol);         // 2nd summand
  nr = polyAdd(s1, s2, tol, w1, w2); // numerator of result
}

template<class T>
void ratPolyNest(
  const std::vector<T>& ni, const std::vector<T>& di,
  const std::vector<T>& po,
  std::vector<T>& nr, std::vector<T>& dr, T tol)
{
  std::vector<T> nt;
  nr = { po[0] };   // numerator of result
  dr = { 1.0 };     // denominator of result
  nt = ni;          // temporary numerator (for convolutive accumulation)
  for(int k = 1; k < po.size(); k++) {
    dr = polyMul(dr, di, tol);  
    nr = polyMul(nr, di, tol);
    nr = polyAdd(nr, nt, tol, 1.0, po[k]);
    nt = polyMul(nt, ni, tol); }
}

template<class T>
void ratNest(
  const std::vector<T>& nI, const std::vector<T>& dI,
  const std::vector<T>& nO, const std::vector<T>& dO,
  std::vector<T>& nR, std::vector<T>& dR, T tol)
{
  std::vector<T> nU, dU, nL, dL;
  ratPolyNest(nI, dI, nO, nU,	dU, tol);  // compute upper num and den
  ratPolyNest(nI, dI, dO, nL,	dL, tol);  // compute lower num and den
  ratDiv(nU, dU, nL, dL, nR, dR, tol);
}
// It's not optimal to call ratPolyNest two times - inside this function, there are values that are
// calculated just the same in both calls, namely the successive powers of ni - but this is not 
// meant to be optimized, high performance code. Maybe in production code, this optimization should
// be done.