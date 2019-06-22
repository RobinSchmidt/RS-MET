double polyEval(std::vector<double>& p, double x)
{
  int k = (int)p.size()-1;  // last valid index
  if(k < 0)
    return 0;
  double y = p[k];
  while(k > 0) {
    k -= 1;
    y = y*x + p[k]; }
  return y;
}

void polyTrunc(std::vector<double>& p, double tol)
{
  int i = (int)p.size()-1;
  while(i > 0) {
    if(fabs(p[i]) <= tol)
      break;
    i -= 1; }
  p.resize(i);
}

double makeMonic(std::vector<double>& p)
{
  double lc = rsLast(p);
  for(size_t i = 0; i < p.size(); i++)
    p[i] /= lc;
  return lc;
}

std::vector<double> polyAdd(
  const std::vector<double>& p, const std::vector<double>& q, 
  double tol, double wp, double wq)
{
  int np = (int) p.size();
  int nq = (int) q.size();
  std::vector<double> r;
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

std::vector<double> polySub(const std::vector<double>& p, const std::vector<double>& q,
  double tol)
{
  return polyAdd(p, q, 1, -1, tol);
}

std::vector<double> polyMul(const std::vector<double>& x, const std::vector<double>& h,
  double tol)
{
  int L = (int)x.size() + (int)h.size() - 1;  // length of result
  std::vector<double> y(L);
  for(int n = 0; n < L; n++) {
    y[n] = 0;  
    for(int k = max(0, n-(int)x.size()+1); k < min(n+1, (int)h.size()); k++)
      y[n] += h[k] * x[n-k]; }
  polyTrunc(y, tol);
  return y;
}

void polyDivMod(std::vector<double> p, std::vector<double> d, 
  std::vector<double>& q, std::vector<double>& r, double tol)
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

std::vector<double> polyDiv(std::vector<double> p, std::vector<double> d, double tol)
{
  std::vector<double> q, r;
  polyDivMod(p, d, q, r, tol);
  return q;
}

std::vector<double> polyMod(std::vector<double> p, std::vector<double> d, double tol)
{
  std::vector<double> q, r;
  polyDivMod(p, d, q, r, tol);
  return r;
}

bool isAllZeros(const std::vector<double>& v, double tol)
{
  for(size_t i = 0; i < v.size(); i++)
    if(fabs(v[i]) > tol)
      return false;
  return true;
}

std::vector<double> polyGCD(
  const std::vector<double>& p, const std::vector<double>& q, double tol, bool monic)
{
  std::vector<double> a = p, b = q, t;
  while(!isAllZeros(b, tol)) {
    t = b;
    b = polyMod(a, b, tol);
    a = t; }
  if(monic)
    makeMonic(a);
  return a;
}

std::vector<double> polyNest(const std::vector<double>& a, const std::vector<double>& b)
{
  int aN = (int)a.size()-1;               // degree of a
  int bN = (int)b.size()-1;               // degree of b
  int cN = aN*bN;                         // degree of result c
  std::vector<double> an(cN+1), c(cN+1); 
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

void ratReduce(
  const std::vector<double>& pIn, const std::vector<double>& qIn,
  std::vector<double>& pOut, std::vector<double>& qOut, double tol)
{
  std::vector<double> gcd = polyGCD(pIn, qIn, tol);
  pOut = polyDiv(pIn, gcd, tol);
  qOut = polyDiv(qIn, gcd, tol);
}

void ratMul(
  const std::vector<double>& p, const std::vector<double>& q,
  const std::vector<double>& r, const std::vector<double>& s,
  std::vector<double>& u, std::vector<double>& v, double tol, bool reduced)
{
  u = polyMul(p, r, tol);
  v = polyMul(q, s, tol);
  if(reduced)
    ratReduce(u, v, u, v, tol);
}

void ratDiv(
  const std::vector<double>& p, const std::vector<double>& q,
  const std::vector<double>& r, const std::vector<double>& s,
  std::vector<double>& u, std::vector<double>& v, double tol, bool reduced)
{
  ratMul(p, q, s, r, u, v, tol, reduced); // r and s are swapped
}

void ratAdd(
  const std::vector<double>& n1, const std::vector<double>& d1,
  const std::vector<double>& n2, const std::vector<double>& d2,
  std::vector<double>& nr, std::vector<double>& dr, 
  double tol, double w1, double w2)
{
  std::vector<double> gcd, f1, f2, s1, s2;
  gcd = polyGCD(d1, d2, tol);
  f1 = polyDiv(d2, gcd, tol);
  f2 = polyDiv(d1, gcd, tol);
  dr = polyMul(f1, d1, tol);
  s1 = polyMul(f1, n1, tol);         // 1st summand in numerator of result
  s2 = polyMul(f2, n2, tol);         // 2nd summand
  nr = polyAdd(s1, s2, w1, w2, tol); // numerator of result
}

void ratPolyNest(
  const std::vector<double>& ni, const std::vector<double>& di,
  const std::vector<double>& po,
  std::vector<double>& nr, std::vector<double>& dr, double tol)
{
  std::vector<double> nt;
  nr = { po[0] };   // numerator of result
  dr = { 1.0 };     // denominator of result
  nt = ni;          // temporary numerator (for convolutive accumulation)
  for(int k = 1; k < po.size(); k++) {
    dr = polyMul(dr, di, tol);  
    nr = polyMul(nr, di, tol);
    nr = polyAdd(nr, nt, tol, 1.0, po[k]);
    nt = polyMul(nt, ni, tol); }
}

void ratNest(
  const std::vector<double>& nI, const std::vector<double>& dI,
  const std::vector<double>& nO, const std::vector<double>& dO,
  std::vector<double>& nR, std::vector<double>& dR, double tol)
{
  std::vector<double> nU, dU, nL, dL;
  ratPolyNest(nI, dI, nO, nU,	dU, tol);  // compute upper num and den
  ratPolyNest(nI, dI, dO, nL,	dL, tol);  // compute lower num and den
  ratDiv(nU, dU, nL, dL, nR, dR, tol);
}
// It's not optimal to call ratPolyNest two times - inside this function, there are values that are
// calculated just the same in both calls, namely the successive powers of ni - but this is not 
// meant to be optimized, high performance code. Maybe in production code, this optimization should
// be done.