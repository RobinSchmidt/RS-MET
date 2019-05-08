
template<class T>
void rsRatioGenerator<T>::fillRatioTable(T* r, int N)
{
  typedef RatioKind RK;
  switch(kind)
  {
  case RK::metallic:      for(int i = 0; i < N; i++) r[i] = metallic(T(i));   break;
  case RK::primeSqrt:     for(int i = 0; i < N; i++) r[i] = primeSqrt(i);     break;
  case RK::primeSqrtDiff: for(int i = 0; i < N; i++) r[i] = primeSqrtDiff(i); break;
  }
}


/*
Ideas:
what about ratios whose continued fraction expansion (CFE) starts with a single n and then has only
1s thereafter - they should be also "very irrational" in the sense of needing a long CFE - i think, 
more so than the other metallic ratios since what matters for the convergence of the CFE are those 
coeffs that come late in the sequence
http://www.peacefromharmony.org/docs/7-27_Stakhov_Math_of_Harmony_EN.pdf
here's a generalization witha 2nd parameter
"Further Generalization of Golden Mean in Relation to Euler's 'Divine' Equation":
https://arxiv.org/ftp/math/papers/0611/0611095.pdf



*/