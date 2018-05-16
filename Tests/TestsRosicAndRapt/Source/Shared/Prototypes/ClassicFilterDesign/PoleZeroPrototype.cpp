template<class T>
void PoleZeroPrototype<T>::butterworth(size_t N, T* k, Complex* p, Complex* z)
{
  size_t L = (N+1)/2;

  // intermediate variables:
  Complex j(0.0, 1.0);               // imaginary unit
  T Gp = sqrt(T(0.5));               // use -3.01 dB point as cutoff frequency for Butterworths
  T ep = sqrt(T(1)/(Gp*Gp)-T(1));    // (1),Eq.2
  T ep_pow = pow(ep, T(-1)/T(N));

  // poles/zeros/gain:
  for(int i = 0; i < L; i++)
  {
    Complex u_i  = (T) (2*i-1) / (T) N;     // Eq.69
    p[i] = ep_pow*j*exp(j*u_i*T(PI*0.5));   // Eq.70
    z[i] = RS_INF(T);                       // zeros are at infinity
  }
  *k = T(1); // is this correct?
}

template<class T>
void PoleZeroPrototype<T>::butterworth(size_t N, T G0, T G, T* k, Complex* p, Complex* z)
{
  if(G0 == T(0) && G == T(1)) { butterworth(N, k, p, z); return; } // delegate to simpler case
  size_t L = (N+1)/2;

  // intermediate variables:
  T GB   = sqrt(G0*G);                           // (2),Eq.52
  T ep   = sqrt( (G*G-GB*GB) / (GB*GB-G0*G0) );  // (2),Eq.12
  T g0   = pow(G0, T(1) / T(N));                 // (2),Eq.94
  T g    = pow(G,  T(1) / T(N));                 // (2),Eq.94
  T wb   = 1.0;                                  // unit cutoff prototype
  T beta = wb * pow(ep, T(-1) / T(N));           // (2),Eq.94

  // poles/zeros/gain:
  T phi, s, c;
  for(int i = 0; i < L; i++)
  {
    phi = T(2*i-1)*T(PI) / T(2*N);           // (2),Eq.95
    s   = sin(phi);                          // (2),Eq.95
    c   = cos(phi);                          // (2),Eq.95
    z[i].real(-s*g*beta/g0);                 // (2),Eq.93
    z[i].imag( c*g*beta/g0);                 // (2),Eq.93
    p[i].real(-s*beta);                      // (2),Eq.93
    p[i].imag( c*beta);                      // (2),Eq.93
  }
  *k = T(1); // is this correct?
}

//-------------------------------------------------------------------------------------------------
/*
template<class T>
void PoleZeroPrototype<T>::elliptic(size_t N, T* k, Complex* p, Complex* z, T ripple, T rejection)
{
  Complex j(0.0, 1.0);                          // imaginary unit
  Complex zeta_i;
  size_t L = (N+1)/2;
  T  Gp  = pow(T(10), -Ap/T(20));               // Eq. 1
  T  Gs  = pow(T(10), -As/T(20));               // Eq. 1
  T  ep  = sqrt(T(1)/(Gp*Gp) - T(1));           // Eq. 2
  T  es  = sqrt(T(1)/(Gs*Gs) - T(1));           // Eq. 2
  T  k1  = ep/es;                               // Eq. 3
  T  k   = ellipdeg(N, k1);                     // solve degree equation for k
  T  v_0 =  (-j*rsAsnC(j/ep, k1)/(T)N).real();  // from ellipap.m



  // calculate the position of the real pole (if present):
  if( r == 1 )
  {
    p[L+r-1] = j * rsSnC(j*v_0, k);                                // from ellipap.m
    z[L+r-1] = RS_INF(T);
  }
  // calculate the complex conjugate poles and zeros:
  for(int i = 0; i < L; i++)
  {
    Complex u_i = (T) (2*(i+1)-1) / T(N);                        // Eq. 69
    zeta_i = rsCdC(u_i, k);                                      // from ellipap.m
    z[i]   = j / (k*zeta_i);                                     // Eq. 62
    p[i]   = j*rsCdC((u_i-j*v_0), k);
  }
}
*/
