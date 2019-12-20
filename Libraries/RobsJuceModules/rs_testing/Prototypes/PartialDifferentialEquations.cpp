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
  typedef rsArray AR;
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
  // mixed, etc.) ..maybe also allow for circular/wrap-around conditions - that's no-physical but
  // may be interesting
}

template<class T>
void rsWaveEquation1D<T>::updateStateArrays()
{
  int N = getNumGridPoints()-1;             // see (1), section 5.2.8
  RAPT::rsArray::copy(&u[0],   &u1[0], N);  // u goes into u1
  RAPT::rsArray::copy(&tmp[0], &u[0],  N);  // tmp goes into u
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

  /*
  // compute Laplacian for interior points:
  for(int i = 1; i < Nx-1; i++) {
    for(int j = 0; j < Ny-1; j++) {
      for(int k = 1; k < Nz-1; k++) {
        T u_xx   = cx * (u(i-1,j,k) - 2*u(i,j,k) + u(i+1,j,k));
        T u_yy   = cy * (u(i,j-1,k) - 2*u(i,j,k) + u(i,j+1,k));
        T u_zz   = cz * (u(i,j,k-1) - 2*u(i,j,k) + u(i,j,k+1));
        L(i,j,k) = u_xx + u_yy + u_zz; }}}

        */
  // todo: compute Laplacian for boundary points....but how? ...maybe using a one-sided 
  // approximation? or maybe using linear extrapolation from the interior points?



}
// see (1), 5.2 for 1D and 10.2 for 2D difference operators
// todo: 
// -implement bi-laplacian (aka biharmonic)...but is this actually a thing in 3D? in 1D and 2D
//  it's used for modeling stiffness...but in 3D? hmm...
// -implement Laplacian for cylidrical and spherical coordinates
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
 studying the relevant modling techniques
-we could actually make the coupling linear as well, by letting the force just be proportional and 
 opposite to the displacement-difference at the contact point - this would correspond to a 
 bi-directional Hooke-spring with equilibirium when it's squeezed into zero space and the strings 
 would be allowed to pass through each other - certainly a very non-physical situation but maybe
 interesting anyway
-to make it more realistic, the strings could also couple to a soundboard at their ends
-i would expect to see interesting musical effects when the lengths (or tensions) of the string are
 tuned to harmonic ratios

*/