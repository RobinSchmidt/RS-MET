template<class T>
void rsHeatEquation1D<T>::setMaxCycleLength(int newLength) 
{ 
  rodArray1.resize(newLength);
  rodArray2.resize(newLength);
  rsArrayTools::fillWithZeros(&rodArray1[0], newLength);
  rsArrayTools::fillWithZeros(&rodArray2[0], newLength);
  reset(); // possibly re-adjust pointers
}

template<class T>
void rsHeatEquation1D<T>::setHeatDistribution(T* d, int N)
{
  rsAssert( N <= (int) rodArray1.size() );
  rodLength = N;

  for(int i = 0; i < N; i++)
    rodIn[i] = rodOut[i] = d[i];
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
  typedef rsArrayTools AR;
  T mean = AR::mean(rodIn, N);
  AR::add(rodIn, -mean, rodIn, N); 
  // ...hmm..this actually just set the mean to zero...which is the most reasonable target value
  // nayway (in case of sound generation)...but what if the user wants some other mean..we'll see

  // set variance to desired target variance:
  // ....
}

//=================================================================================================

template<class T>
T rsWaveEquation1D<T>::getCourantNumber(T timeStep) const
{
  int N = getNumGridPoints()-1; // see (1), section 5.2.8
  T c   = waveSpeed;            // what unit?
  T L   = T(1);                 // we use the unit interval [0,1], so the spatial length is 1
  T gam = c / L;                // "gamma" - (1), Eq. 6.5
  T k   = timeStep;             // temporal sampling interval
  T h   = L / N;                // spatial sampling interval
  return gam*k / h;             // "lambda", the Courant number
}
// see:
// https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition

template<class T>
T rsWaveEquation1D<T>::getOmegaForWaveNumber(T waveNumber, T timeStep) const
{
  T beta   = waveNumber;
  T k      = timeStep;
  T h      = getGridSpacing();
  T lambda = getCourantNumber(timeStep);
  return (T(2)/k) * asin(lambda * sin(beta*h/2)); // (1), Eq. 6.43
  // does this depend on the scheme? ...probably - we may need different formulas for different 
  // schemes
}

template<class T>
void rsWaveEquation1D<T>::setInitialConditions(T* newPositions, T* newVelocities,
  int length, T timeStep)
{
  rsAssert(length == getNumGridPoints(), "array length should match number of grid points");
  T k = timeStep;
  for(int l = 0; l < length; l++) {
    u[l]  = newPositions[l];
    u1[l] = u[l] + k*newVelocities[l];  // (1), Eq 6.36 ...why + not -?
  }
}

template<class T>
void rsWaveEquation1D<T>::updateState(T timeStep)
{
  computeInteriorPoints(timeStep);
  //computeInteriorPointsSimple();      // simplified formula for Courant number == 1
  computeBoundaryPoints(timeStep);
  updateStateArrays();
}

template<class T>
void rsWaveEquation1D<T>::computeInteriorPoints(T timeStep)
{
  // We implement the scheme in (1), Eq. 6.34: u_tt = g^2 * u_xx where u_tt and u_xx are central 
  // difference approximations to the second temporal and spatial derivative respectively. g is a 
  // constant (gamma).

  // these intermediate variables are mostly for clarity and consistency with the mathematical 
  // notation in (1) - production code may get away without them ...or maybe the compiler 
  // optimizes them away anyway:
  int N = getNumGridPoints()-1; // see (1), section 5.2.8
  T c   = waveSpeed;     // what unit?
  T L   = T(1);          // we use the unit interval [0,1], so the spatial length is 1
  T gam = c / L;         // "gamma" - (1), Eq. 6.5
  T k   = timeStep;      // temporal sampling interval
  T h   = L / N;         // spatial sampling interval
  T lam = gam*k / h;     // "lambda", the Courant number  use getCourantNumber
  rsAssert(lam <= T(1), "scheme unstable with these settings!"); // (1), Eq. 6.40
  // actually lamda == 1 is most desirable - in this special case, the numerical solution becomes 
  // exact
  // maybe factor out into getLamda() or getCourantNumber(timeStep)

  // compute updated solution at interior points (factor out to allow switching between schemes):
  T l2 = lam*lam;  // lambda-squared
  for(int l = 1; l <= N-1; l++)                           // (1), section 5.2.8
    tmp[l] = 2*(1-l2)*u[l] + l2*(u[l-1]+u[l+1]) - u1[l];  // (1), Eq 6.35
}
// todo: allow to switch between different schemes
// it's a bit surprising (in a good way), that a central difference in time leads to an explicit 
// rather than implicit recursion
// -> figure out, why -> derive the recursion from the operators
// other schemes are perhaps only of academic interest because this scheme is actually the best
// in various respects (see (1))

template<class T>
void rsWaveEquation1D<T>::computeInteriorPointsSimple()
{
  for(size_t l = 1; l < u.size()-1; l++) 
    tmp[l] = u[l+1] + u[l-1] - u1[l];    // (1), Eq. 6.54
}

template<class T>
void rsWaveEquation1D<T>::computeBoundaryPoints(T timeStep)
{
  int N  = getNumGridPoints()-1; // see (1), section 5.2.8
  tmp[0] = tmp[N] = T(0);        // endpoints fixed at zero - "Dirichlet" conditions
  // todo: allow to let client code choose from various boundary conditions (Dirichlet, Neumann, 
  // mixed, etc.) ..maybe also allow for circular/wrap-around conditions - that's non-physical but
  // may be interesting - it would sort of model a circular string - maybe like strung around a 
  // cyclinder

  // on page 138, (1) says something about "at l=0, scheme 6.35 may be modified to..." - try that
}

template<class T>
void rsWaveEquation1D<T>::updateStateArrays()
{
  int N = getNumGridPoints()-1;             // see (1), section 5.2.8
  RAPT::rsArrayTools::copy(&u[0],   &u1[0], N);  // u goes into u1
  RAPT::rsArrayTools::copy(&tmp[0], &u[0],  N);  // tmp goes into u
}

// todo: implement 6.38, 6.45, 146: bottom, 149: u_l^{n+1} =..., 6.59, 
// 151: implicit scheme, 6.62, 6.66: recursion

// see:
// https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition



//=================================================================================================

template<class T>
T rsRectangularMembrane<T>::getCourantNumber() const
{
  int N    = u.getNumRows()-1; // later have separate Nx, Ny
  T h      = T(1)/N;           // later have separate hx, hy
  T k      = timeStep;
  T gamma  = waveSpeed;        // gamma = waveSpeed/length and length is normalized to unity
  T lambda = k*gamma/h;        // same as in 1D case, (1), pg 310 
  return lambda;
}


template<class T>
void rsRectangularMembrane<T>::updateState()
{
  computeInteriorPoints();
  computeBoundaryPoints();
  updateStateMatrices();
}

template<class T>
void rsRectangularMembrane<T>::computeInteriorPoints()
{
  int N    = u.getNumRows()-1;    // later have separate Nx, Ny
  T lambda = getCourantNumber();

  // stability limit: lambda = 1/sqrt(2) - maybe add an rsAssert? but maybe assert that 
  // l2=lambda^2 <= 0.5 because 0.5 is exactly representable

  // compute new displacements via (1), Eq. 11.10:
  T l2 = lambda*lambda;
  for(int l = 1; l <= N-1; l++)
    for(int m = 1; m <= N-1; m++)
      tmp(l,m) = l2*(u(l+1,m) + u(l-1,m) + u(l,m+1) + u(l,m-1)) + 2*(1-2*l2)*u(l,m) - u1(l,m);

  // how are we supposed to incorporate the aspect-ratio "epsilon"?
}
// -maybe rename to computeInteriorPointsExplicit - indicating that this uses the explicit scheme, 
//  then also implement implicit scheme(s)
// -maybe also implement different explicit schemes that use different approximations to the 
//  Laplacian operator - this scheme here uses a 5-point stencil using direct neighbours, there's
//  also a scheme that uses diagonal neighbours and one that uses a mix of both - compare the 
//  schemes in terms of numerical dispersion and (an)isotropy of wave-propagation


template<class T>
void rsRectangularMembrane<T>::computeInteriorPointsSimple()
{
  // implements the simplified scheme for special case lambda = 1/sqrt(2) - Eq. 11.12 
  int N = u.getNumRows()-1;
  for(int l = 1; l <= N-1; l++)
    for(int m = 1; m <= N-1; m++)
      tmp(l,m) = 0.5*(u(l+1,m) + u(l-1,m) + u(l,m+1) + u(l,m-1)) - u1(l,m);
}
// not yet tested

// todo: implement schemes 11.16, 11.18 (implicit)
// make a class for circular membranes and implement schemes 11.25, 11.27

template<class T>
void rsRectangularMembrane<T>::computeBoundaryPoints()
{
  // fix ends to zero (todo: let user switch between different boundary conditions): 
  int M = tmp.getNumRows();
  int N = tmp.getNumColumns();
  int m, n;
  for(m = 0; m < M; m++) tmp(m,   0  ) = T(0);
  for(m = 0; m < M; m++) tmp(m,   N-1) = T(0);
  for(n = 0; n < N; n++) tmp(0,   n  ) = T(0);
  for(n = 0; n < N; n++) tmp(M-1, n  ) = T(0);
  // ...well...actually, if the points are initialized to zero, wo don't really need to do anything
  // here as they are not modified in computeInteriorPoints - but for other boundary conditions, 
  // the situation may be different

  // if we would use circular/wrap boundary conditions for x,y, we would get a toroidal topology,
  // i.e. the membrane would model the surface of a torus - how abut modeling the surface of a 
  // sphere - maybe use wrap-around in spherical coordinates? hmm...that seems appropriate for the 
  // meridians but not the elevation (or whatever it is called)...maybe the poles have to be 
  // treated in a special way


  // test: place an "obastacle" into the membrane:
  tmp(20, 30) = 0;
  tmp(20, 31) = 0;
  tmp(20, 32) = 0;
  tmp(20, 33) = 0;
  tmp(21, 31) = 0;
  tmp(22, 31) = 0;
}
// todo: maybe introduce "obstacles" - set the displacement to zero at some arbitrary points inside
// the interior region - maybe this can be realized using a "mask" - a matrix with multipliers that 
// scale the displacements at every grid point - values of 1 do nothing, values of 0 fix the 
// displacements at these points and values in between do something in between - a user could even
// "draw" such a mask

template<class T>
void rsRectangularMembrane<T>::updateStateMatrices()
{
  u1 = u;
  u  = tmp;
}
// can be optimized - use class rsMatrix


//=================================================================================================
/*
template<class T>
rsRectangularRoom<T>::rsRectangularRoom(int Nx, int Ny, int Nz)
{
  // maybe factor out into function setGridResolution:
}
*/

template<class T>
void rsRectangularRoom<T>::setGridDimensions(int _Nx, int _Ny, int _Nz)
{
  Nx = _Nx;
  Ny = _Ny;
  Nz = _Nz;
  std::vector<int> shape({ Nx, Ny, Nz });
  u.setShape(shape);
  u_t.setShape(shape);
  u_tt.setShape(shape);

  // maybe we should have function template that can be called like:
  // u.setShape(Nx,Ny,Nz); ...needs to be a variadic template
}


template<class T>
void rsRectangularRoom<T>::updateState()
{
  computeLaplacian3D(u, u_tt);




  T k = timeStep;

  // update the "velocities" by adding a fraction of the Laplacian (which is proportional to
  // the "acceleration"):
  int N = u.getSize();
  rsArrayTools::addWithWeight(u_t.getDataPointer(), N, u_tt.getDataPointer(), k); // is k the right scaler?

  // maybe we should set the velocities to zero at the boundary...currently they already are 
  // because the lapalcian is only computed at inetrior points

  // update the pressures by adding a fraction of the "velocities":
  rsArrayTools::addWithWeight(u.getDataPointer(), N, u_t.getDataPointer(), k); // is k the right scaler?

  int dummy = 0;
}

template<class T>
void rsRectangularRoom<T>::reset()
{
  u.setToZero();
  u_t.setToZero();
  u_tt.setToZero();
}

template<class T>
void rsRectangularRoom<T>::computeLaplacian3D(const rsMultiArray<T>& u, rsMultiArray<T>& L)
{
  // In 1D, a 2nd order accurate approximation to the 2nd spatial derivative of a function u(x,t) 
  // is given by:
  //   u_xx ~= (u(x-h) - 2*u(x) + u(x+h)) / h^2
  // where h is the grid spacing. See (1), Eq. 5.12a. Here in 3D, we use the same formula for all
  // 3 dimensions to get approximations to the 2nd partial derivatives along the 3 spatial 
  // directions and add them all up to get an approximation of the Laplacian.

  // todo: there are other ways to approximate the Laplacian - this here is the simplest - maybe 
  // implement various variants...

  // spatial sampling intervals:
  T hx = Lx / (Nx-1);
  T hy = Ly / (Ny-1);
  T hz = Lz / (Nz-1);

  // scaling coeffs used in the finite difference approximations:
  T cx = 1 / (hx*hx);
  T cy = 1 / (hy*hy);
  T cz = 1 / (hz*hz);

  // compute Laplacian for the interior points of the cuboid:
  for(int i = 1; i < Nx-1; i++) {
    for(int j = 1; j < Ny-1; j++) {
      for(int k = 1; k < Nz-1; k++) {
        T u2     = 2*u(i,j,k);
        T u_xx   = cx * (u(i-1,j,k) - u2 + u(i+1,j,k));
        T u_yy   = cy * (u(i,j-1,k) - u2 + u(i,j+1,k));
        T u_zz   = cz * (u(i,j,k-1) - u2 + u(i,j,k+1));
        L(i,j,k) = u_xx + u_yy + u_zz; }}}

  // compute Laplacian for the boundary planes (faces of the cuboid):
  // ...this will need 6 double-loops

  // compute Laplacian for the edges of cuboid:
  // ...this requires 12 single-loops

  // compute Laplacian in the corners
  // ...requires 8 assignments



  // todo: compute Laplacian for boundary points....but how? ...maybe using a one-sided 
  // approximation? or maybe using linear extrapolation from the interior points? ...maybe first
  // implement a 1D and 2D version of the laplacian
  // maybe for fixed boundary conditions, we don't need it...hmm..but what actually are appropriate
  // boundary conditions for the pressure at the room walls? intuitively, the sound-velocity must 
  // go zero at the boundary...does that mean, the pressure gradient must go to zero?
}
// see (1), 5.2 for 1D and 10.2 for 2D difference operators
// todo: 
// -implement bi-laplacian (aka biharmonic)...but is this actually a thing in 3D? in 1D and 2D
//  it's used for modeling stiffness...but in 3D? hmm...maybe for modeling stiffness in solids?
// -implement Laplacian for cylindrical and spherical coordinates
// -maybe rename this function to reflect that we a doing a 7 point approximation in cartesian 
//  coordinates - there are so many other possibilities...




//=================================================================================================

/*
 Background:



 questions:
 -the so called dispersion-relations seem to always relate frequency and wavelength (right?) - and 
  frequency and wavelength in turn are related via wavelength * frequency = wavespeed, see here:
  https://en.wikipedia.org/wiki/Wavelength#Sinusoidal_waves
  ...so should we think of the dispersion relation as a function that assigns a wave-propagation 
  speed to each frequency? and for non-dispersive wave-propagation that wave-speed is constant 
  (seen as function of frequency)


Schemes:
-a "scheme" for solving a PDE for some function u(x,t) = F(u,u_x,u_t,u_xx,u_tt,...) numerically 
 involves using sampled values of the function u(x,t-k),u(x,x-2k),u(x-h,t-k),...
-intuitively, when we have a PDE like u_tt = c*u_xx, we may express it as:
 u(x,t) = integral_0^s v(x,s) ds, v(x,t) = integral_0^t a(x,s) ds for acceleration and velocity
 and our PDE (in this case, the wave-equation) may say: a(x,t) = c*u_xx(x,t) which we may calculate 
 directly unsing finite differences in space, then we may update the velocity v by adding k*a, then 
 update the displacement u by adding k*v - this would amount to forward Euler integration in time
-all of the intermediate variables a(x,t), v(x,t) can also be numerically approximated by finite 
 differences in time, for example v(x,t) ~= u(x,t) - u(x,t-k), a(x,t) ~= v(x,t) - v(x,t-k), so in 
 the end, we may express everything in terms of delayed values of u(x,t) itself without resorting
 to intermediates like v and a. This expression in terms of u(x,t) and its delayed and shifted 
 values is a "scheme"
-todo: figure out what the scheme for this simple, intuitive (any maybe naive) 
 twice-forward-Euler-in-time strategy would look like and compare it to other schemes that are in 
 common use - how can these other schemes be "translated back" into the language of velocities and 
 accelerations?
-in some schemes, the values to be computed (typically, u(x,t+k)) along the grid 
 ...x-2h,x-h,x,x+h,... may not be isolated on the left hand side and it may not be possible to do so 
 -in this case, the scheme is implicit
 -in implicit schemes, we may have to solve a system of equations for all values u(x+m*h,t+k) at 
  once where index m runs over the spatial samples along the grid
  -it may make sense to further subdivide implicit schemes depending on whether the system of 
   equations to be solved is linear or nonlinear
  -in the linear case, we typically deal with sparse (often band-diagonal) matrices
 -implicit schemes are sometimes easier to analyze in terms of numerical stability and they also to
  tend to *be* more stable
 -hmm...i actually think, it's more like the grid values at the next time-step are all 
  interdependent and it's not possible to solve for a single grid-value explicitly in terms of 
  known values
-i think, the difference between using a scheme with delayed values and explicit representations of 
 deriviatives is somewhat similar to consolidating two first order (integrator) filters into a 
 single second order (biquad) filter - the biquad may produce the same output but one loses 
 physical interpretability of its internal state variables

-maybe instead of using past values of the displacement u itself with a temporal stencil, 
 explicitly compute and use things like u_xx, u_t, u_tt, etc. (for u_xxxx also use the notation 
 u_4x) that is more readily interpretable and mor intuitive for research purposes - only at the 
 very end, when everything has settled, convert to a scheme using u1, u2 for delayed displacement 
 ..or maybe use notation u_d, u_dd, u_3d, u_4d, etc. - subscript d meaning delayed
-with this representation, it is also intuitively obvious, how to include driving forces - they are
 just added to the acceleration (with a proportionality factor) - coupling forces between multiple
 systems can also be represented in a straightforward manner
-re-derive the schemes from (1) using that notation with finite difference approximations like:
 u_t ~= (u - u_d) / k for the forward-difference approximation, u_t ~= (u_f - u_d) / 2k (where u_f 
 is a future value) for the central difference, etc.
-maybe we can derive stability criteria for the representation via u, u_t, u_tt - maybe consider 
 the total energy = sum(u^2)/2 + sum(u_t^2)/2 = E_pot + E_kin - for a stable scheme, the derivative of 
 this energy should be strictly nonpositive - for a conservative scheme, the derivative should be 
 zero
-to (very crudely) simulate a rattling element on a string or membrane, hardclip the displacements
 at the position(s) of the rattling element - to refine it, maybe use smoother functions or model 
 it as forces (against the displacement) that acts only when the displacement is above some 
 threshold


sort of pre-print of (1):
https://ccrma.stanford.edu/~bilbao/nssold/booktoplast
it has some extra chapters that are not in the book - about finite element and spectral methods

https://www2.ph.ed.ac.uk/~sbilbao/matlabpage.html
https://www2.ph.ed.ac.uk/~sbilbao/nsstop.html

http://dafx09.como.polimi.it/proceedings/papers/paper_68.pdf
https://www2.ph.ed.ac.uk/~sbilbao/CFApaper.pdf

https://hplgit.github.io/fdm-book/doc/pub/book/sphinx/index.html
https://hplgit.github.io/fdm-book/doc/pub/book/sphinx/._book008.html
http://hplgit.github.io/num-methods-for-PDEs/doc/pub/wave/sphinx/._main_wave005.html



Notes from Numerical Sound Synthesis:
Chapter 5:
-a general PDE in 1D for a function u(x,t) can be expressed as F(x,t,u,u_x,u_t,u_xx,u_tt,...) = 0
-if F does not explicitly depend on x, it is shift-invariant. shift-invariance models uniform 
 material properties
-if F does not explicitly depend on t, it is time-invariant. time-invariance models static 
 operating (playing) conditions
-in the book, that function F is expressed as an operator P = P(d/dt,d/dx,d^2/dt^2,d^2/dx^2,...) 
 acting on u(x,t): P u = 0 (Eq. 5.1)...but i think, that notation is a bit obscure...
-linear time-invariant and linear shift-invariant systems are abbreviated as LTI and LSI
-in addition to the PDE, boundary conditions and initial conditions must be supplied. The PDE 
 itself determines an infinite set of possible solutions, the initial- and boundary conditions pick
 one solution from that set
-PDEs in musical acoustics usually: 
 -are 2nd order in time (i.e. u_t and u_tt appears but no higher time-derivatives)
 -involve even order spatial derivatives such as u_xx, u_xxxx which reflects direction independence
 -are of the hyperbolic type (what does that mean?)
-the wave-equation in 1D is: u_tt = g^2 * u_xx, where g^2 (in the book: gamma-squared) is the 
 square of the wave-velocity (or its reciprocal? -> figure out!)
 -the wave equation is linear, iff g does not depend on u
 -if g depends on x, the system is not LSI
 -if g depends on t, the system is not LTI
-von Neumann analysis (an extension of the z-transform to distributed systems) is only applicable 
 to LTI/LSI systems
-LTI but non-LSI systems still allow analysis in terms of modes and frequencies but there's no 
 global wave-velocity
-for analyzing PDEs, one may apply the Laplace transform in the time domain:
   û(x,s)  = integral u(x,t) * exp(-s*t) dt     Eq. 5.2
 or the Fourier transform in the spatial domain:
   u'(b,t) = integral u(x,t) * exp(-j*b*x) dx   Eq. 5.3, b is used for beta
 or both, where the integrals both run from -inf to inf. For problems over a finite spatial domain, 
 Fourier series  may be used. The derivative operators transform to multiplications to 
 corresponding powers of s or j*b, for example: u_tt L-> s^2 * û, u_xx F-> (j*b)^2 * u', where 
 L-> means Laplace-transforms-to, likewise for F->
-one may also plug the ansatz: u(x,t) = exp(s*t + j*b*x) into the PDE and analyze, under which 
 conditions such a special function (representing a wave of frequency s and and wavenumber b)
 is a solution of the PDE. ...doing this will lead to constraints for how frequencies and 
 wavenumbers must be related? or what? -> try it!
-If the PDE is LSI, both Laplace and Fourier transform may be used:
   P u = 0  L,F->  ^P' û' = 0  (the "hat" should actually appear above P like it does for u)
 where ^P' = ^P'(s, j*b) is a bivariate polynomial in s and j*b and often called the "symbol" of 
 the PDE. The condition that ^P' û' = 0 means that a nontrivial solution must satisfy ^P' = 0 (why? 
 because otherwise û' would have to be identically zero?), so we are looking for zeros of the 
 bivariate polynomial ^P'(s, j*b)
 -solutions/zeros will be of the form s_p = s_p(j*b) and we will get M such solutions, where M is 
  the highest temporal derivative - whih is typically 2, so the solutions are often written as 
  s+, s-. These two solution correspond to waves travelling left and right
 -these solutions (solution curves?) are called dispersion relations
 -in lossless systems, the solutions are typically of the form: s_p(j*b) = j*w_p(j*b), so one may 
  write w_p = w_p(b), i.e. the temporal frequency w can be written as a function of the spatial 
  frequency b. In this case, phase velocity v_p and group velocity v_g are defined as v_p = w/b and
  v_g = dw/db
-Energy analysis is based on spatial inner products, i.e. integrals over a product of two functions
 f and g of the form <f,g> = integral f*g dx where the integral runs over some finite or infinite 
 domain. The integration over space removes spatial dependence - the energies become functions of 
 time alone. If f and g are vector valued functions (as is the case for non-planar 
 string-vibration), the product f*g inside the integral must be  replaced by a scalar product aka 
 dot-product. [todo: find examples for kinetic and potential energies expressed this way]


todo: von Neumann analysis, 
  numerical schemes/energy/dispersion/boundary-conditions




Notes on the 3D wave equation

in Ingenieurakustik, it's derived in an interesting way:
  rho * v_t = -p_x           Eq. 1.1
  rho * v_x = -p_t / c^2     Eq. 1.8
so, it's expressed as a system of two first order PDEs (where does *that* derive from?). rho is the 
density (i think), p is the pressure and v is the velocity. From this system of two first order 
PDEs, we get a second order PDE for p by differentiating the 1st equation with respect to x and the 
2nd with respect to t:
  p_xx = p_tt / c^2          Eq. 1.9
We may also get a 2nd order PDE for v by differentiating (1.1) with respect to t and (1.8) with 
respect to x:
  v_xx = v_tt / c^2          Eq. 1.10
So far, this applies to the 1D wave-equation - the generalization to 3D is:
  rho   v_t    = -grad(p)       Eq. 1.14
  rho * div(v) = -p_t / c^2     Eq. 1.15
where v is now a vector field, div(v) is its divergence (= vx_x + vy_y + vz_z where vx_x means the 
x-derivative of the x-component of the velocity, etc.) and grad(p) is the gradient (p_x, p_y, p_z) 
of the scalar pressure field p. They also introduce a velocity potential phi (i call it V) which is 
a scalar field with the property grad(V) = v - the vector v is the gradient of the scalar V. The 
system of two 1st order PDEs can again be written as a single 2nd order PDE either in the pressure 
p or the velocity potential V:
  Lap(p) = p_tt / c^2    Eq. 1.17
  Lap(V) = V_tt / c^2    Eq. 1.16
where Lap(...) means the Laplacian operator: Lap(p) = p_xx + p_yy + p_zz


-In 1 time step, an impulse can travel at most 1 spatial step. That implies that if space is sampled 
 so densely that a physical wave would travel multiple space intervals in one time step, i can't work
 anymore. And that's the inuitive reason why instabilities occur - we have to increase the temporal
 sample-rate as well, wehn we want to increase the spatial sampling density.

-A PDE solver could be parametrized by a function (or better: functor) that can compute partial 
 derivatives.

-Maybe implement the finite volume method with rectangular cells - the function u(x,t) that 
 represents the conserved quantity could be seen as a fluid density.

-To approximate the Laplacian numerically, we could take a weighted average of the straight 
 approximant and the diagonal approximant with weights inversely proprtional to distance, i.e.
 1 and 1/sqrt(2).

-We could implement an explicit scheme by keeping track of the physical variables v(x,t), a(x,t) 
 for velocity and acceleration. Such a scheme could perhaps be artificially stabilized by computing
 the total energy and then scaling velocity and displacement so as to keep that constant..or we 
 could actually apply an arbitrary external envelope to the total energy. ...but this energy 
 normalization will not suppress oscillations at the Nyquits freq. If such oscillation occur in 
 space and time, it would probably make most sense to suppress them with a Nyquist blocker in the 
 time domain (i've already experimented with applying such a filter in the spatial domain with the
 Schrödinger equation - but time domain seems to be the better approach - it would replace Euler 
 steps with trapezoidal steps, i think)



  
by


Ideas:




-imagine two strings that cross each other at a common point, maybe one goes along the x-direction 
 and the other along the y-direction, so the physical configuration could look like this:

     |
     |
  ---+------------
     |

-now consider displacing the upper string in the z-direction (out of the screen) and letting it go, 
 let's assume, the upper string is one in the x-direction
-until the strings collide, it behaves the same as if the other string wasn't there
-let's call u the z-displacement of the upper string and v the z-displacement of the lower string
-the condition that u is above v means that u >= v at the contact point, if u < v, there should be
 suddenly a very strong force pushing u up and v down - that's the model of the collision - perhaps
 some function with a threshold and a very steep buildup beyond the threshold...maybe it should also
 act in some neighbourhood of the contact point (due to stiffness)...maybe some sort of bell-shape
-in a refined model, we could also consider displacement in the x and y directions - a traverse 
 displacement in the y direction of the horizontal string would couple to longitudinal motion in the
 vertical string and vice versa, when the strings are in contact
-we could consider different angles at which the strings cross - but for modeling only the z 
 displacement, that angle should be irrelevant
-we would have two 1D PDEs coupled by a nonlinear mechanism, so it's a nice starting point for 
 studying the relevant modeling techniques
-we could actually make the coupling linear as well, by letting the force just be proportional and 
 opposite to the displacement-difference at the contact point - this would correspond to a 
 bi-directional Hooke-spring with equilibirium when it's squeezed into zero space and the strings 
 would be allowed to pass through each other - certainly a very non-physical situation but maybe
 interesting anyway
-to make it more realistic, the strings could also couple to a soundboard at their ends
-i would expect to see interesting musical effects when the lengths (or tensions) of the string are
 tuned to harmonic ratios

*/