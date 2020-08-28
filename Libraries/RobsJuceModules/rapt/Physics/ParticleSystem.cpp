template<class T>
rsVector3D<T> rsParticle<T>::getGravitationalFieldAt(rsVector3D<T> p, T cG)
{
  rsVector3D<T> r = p - pos;     // vector from this particle to p
  T d = r.getEuclideanNorm();    // distance
  r /= d;                        // r is now normalized to unit length
  return -(cG*mass/(d*d)) * r;   // inverse square law
  // see https://en.wikipedia.org/wiki/Gravitational_field
}

template<class T>
rsVector3D<T> rsParticle<T>::getElecricFieldAt(rsVector3D<T> p, T cE)
{
  rsVector3D<T> r = p - pos;    // vector from this particle to p
  T d = r.getEuclideanNorm();   // distance
  r /= d;                       // r is now normalized to unit length
  return (cE*charge/(d*d)) * r; // inverse square law
  // see (3), page 88, Eq. 7.3 - true only for static fields
}

template<class T>
rsVector3D<T> rsParticle<T>::getMagneticFieldAt(rsVector3D<T> p, T cM)
{
  rsVector3D<T> r = p - pos;                // vector from this particle to p
  T d = r.getEuclideanNorm();               // distance
  r /= d;                                   // r is now normalized to unit length
  return (cM*charge/(d*d)) * cross(vel, r); // inverse square law
  // see http://www.phys.uri.edu/gerhard/PHY204/tsl210.pdf
  // true only for static fields
}

// a lot of duplicated code among these 3 functions (only the last line is different), maybe 
// factor out the common code

template<class T>
T rsParticle<T>::getGravitationalPotentialAt(rsVector3D<T> p, T cG)
{
  rsVector3D<T> r = p - pos;
  T d = r.getEuclideanNorm();
  return -cG*mass/d;
  // see https://en.wikipedia.org/wiki/Gravitational_potential#Mathematical_form
}

template<class T>
T rsParticle<T>::getElectricPotentialAt(rsVector3D<T> p, T cE)
{
  rsVector3D<T> r = p - pos;
  T d = r.getEuclideanNorm();
  return cE*charge/d;         
  // (3), page 96, Eq. 7.19
}

template<class T>
rsVector3D<T> rsParticle<T>::getMagneticPotentialAt(rsVector3D<T> p, T cM)
{
  rsVector3D<T> r = p - pos;
  T d = r.getEuclideanNorm();
  return cM*charge*vel/d; 
  // (2), page 15.8, Eq. 15.24 (with j = charge*vel, integral can be ignored)
}

template<class T>
T rsParticle<T>::getGravitationalEnergy(const rsParticle<T>& p, T cG)
{
  T potential = getGravitationalPotentialAt(p.pos, cG);
  return potential * p.mass; // maybe * 0.5
  // 0.5, because we assign one half to one particle and the other half to the other
  // ...naah...we don't use this convention here
}

template<class T>
T rsParticle<T>::getElectricEnergy(const rsParticle<T>& p, T cE)
{
  T potential = getElectricPotentialAt(p.pos, cE);
  return potential * p.charge; // (2), page 8.10, Eq. 8.28 or page 15.7, Eq. 15.21
  // we don't use the 0.5 factor, we assign all energy to passed particle p instead of distributing
  // it between p and this equally
}

template<class T>
T rsParticle<T>::getMagneticEnergy(const rsParticle<T>& p, T cM)
{
  rsVector3D<T> potential = getMagneticPotentialAt(p.pos, cM);
  return p.charge * dot(p.vel, potential); // (2), page 15.6, Eq. 15.20
  // factor 0.5 left out, see comments above
}

//-------------------------------------------------------------------------------------------------

template<class T>
rsParticleSystem<T>::rsParticleSystem(int numParticles)
{
  setNumParticles(numParticles);
}

// Setup:

template<class T>
void rsParticleSystem<T>::setNumParticles(int newNumParticles)
{
  particles.resize(newNumParticles);
  forces.resize(newNumParticles);
  initialPositions.resize(newNumParticles);  
  initialVelocities.resize(newNumParticles);
}

// Inquiry:

template<class T>
T rsParticleSystem<T>::getKineticEnergy()
{
  T E = 0;
  for(size_t i = 0; i < particles.size(); i++)
    E += particles[i].getKineticEnergy();
  //return E;
  return E/stepSize; // plots suggest that we need to scale it that way
}

template<class T>
T rsParticleSystem<T>::getGravitationalPotentialEnergy()
{
  T E = 0;
  for(size_t i = 0; i < particles.size(); i++)
    for(size_t j = i+1; j < particles.size(); j++)
      E += particles[i].getGravitationalEnergy(particles[j], cG);
  return E;
}

template<class T>
T rsParticleSystem<T>::getElectricPotentialEnergy()
{
  T E = 0;
  for(size_t i = 0; i < particles.size(); i++)
    for(size_t j = i+1; j < particles.size(); j++)
      E += particles[i].getElectricEnergy(particles[j], cE);
  return E;
}

template<class T>
T rsParticleSystem<T>::getMagneticPotentialEnergy()
{
  T E = 0;
  for(size_t i = 0; i < particles.size(); i++)
    for(size_t j = i+1; j < particles.size(); j++)
      E += particles[i].getMagneticEnergy(particles[j], cM);
  //return E;
  return E/stepSize; // plots suggest that we need to scale it that way
}

template<class T>
T rsParticleSystem<T>::getPotentialEnergy()
{
  return getGravitationalPotentialEnergy() // preliminary - uses 1/d^2 force law
    + getElectricPotentialEnergy() 
    + getMagneticPotentialEnergy();


  /*
  // This ist still incorrect. In particuar, it doesn't take into account our tweaks to the 
  // force-law. I think, to do so, we have to express potetial energy as integral of the force.
  // see (4), page 96
  T E = 0;
  for(size_t i = 0; i < particles.size(); i++)
  {
    for(size_t j = i+1; j < particles.size(); j++)
    {
      // use the 2nd term in (1), Eq. 13.14
      T r = (particles[j].pos - particles[i].pos).getEuclideanNorm(); // r_ij

      // gravitational term:
      E += -cG * particles[i].mass * particles[j].mass / r;

      // electrical term (verify this - formula just inferred by analogy):
      E +=  cE * particles[i].charge * particles[j].charge / r;

      // todo: figure out, how the energy changes when we use different force-laws - use
      // W = F * s ..or W = integral over F ds
    }
  }
  return E;
  */
}

template<class T>
rsVector3D<T> rsParticleSystem<T>::getTotalMomentum()
{
  rsVector3D<T> p;
  for(size_t i = 0; i < particles.size(); i++)
    p += particles[i].getMomentum();
  return p;
}

// Processing:

template<class T>
T rsParticleSystem<T>::getForceScalerByDistance(T d, T size1, T size2)
{
  //return 1 / (d*d*d);          // physical law

  //return 1 / (c + pow(d,p));     // 1 / (c + d^p) ...seems stable with c=1, limit (d=0): 1/c
  //return 1 / pow(c+d,p);         // 1 / (c + d)^p, limit: 1/c^p
  //return pow((c+1)/(c+d), p);    // ((c+1)/(c+d))^p, limit: ((c+1)/c)^p

  //return (1 + size1 + size2) / (size1 + size2 + pow(d,p)); 
  //return 1 / (size1 + size2 + pow(d,p)); 
  return 1 / (rsMax(size1 + size2, pow(d,p))); 
  //return 1 / (max(size1,  size2, pow(d,p))); 
}

template<class T>
T rsParticleSystem<T>::getForceByDistance(T d, T size1, T size2)
{
  return d * getForceScalerByDistance(d, size1, size2);
}

template<class T>
rsVector3D<T> rsParticleSystem<T>::getForceBetween(const rsParticle<T>& p1, const rsParticle<T>& p2)
{
  // precomputations:
  rsVector3D<T> r = p2.pos - p1.pos;    // vector pointing from p1 to p2

  //// old - physically correct:
  //T d = r.getEuclideanNorm();   // distance between p1 and p2 == |r|
  //T s = 1 / (d*d);              // reciprocal of |r|^2 - used as multiplier in various places
  //r  /= d;                      // r is now normalized to unit length (for k=0)

  // new:
  T d = r.getEuclideanNorm();        // distance between p1 and p2 == |r|
  T s = getForceScalerByDistance(d, p1.size, p2.size);


  // compute the 3 forces:
  rsVector3D<T> f;
  f += s*cG * p1.mass   * p2.mass   * r;                                // gravitational force
  f -= s*cE * p1.charge * p2.charge * r;                                // electric force

  // magnetic force is still experimental:
  //f += s*cM * p1.charge * p2.charge * cross(p1.vel, cross(p2.vel, r));  // magnetic force
  //f += s*cM * p1.charge * p2.charge * cross(p2.vel, cross(p1.vel, r));  // or should it be this way? - seems to make no difference, but it's not generally the same
  f -= s*cM * p1.charge * p2.charge * cross(p1.vel, cross(p2.vel, r));  
    // minus sign gives qualitatively more reasonable results: two particles heading off parallel 
    // into the same direction which are only subject to the magnetic force, attract each other
    // as said in "Ein Jahr für die Physik", page 120 ...nope, the minus sign is alright: the 
    // formula in the book is for the force on p2 due to p1 (we need the force on p1 due to p2)
    // ...but the formula also has v2 x (v1 x r) and we do v1 x (v2 x r) ...how can we figure out
    // if the order of the cross-products is right? (it doesn't seem to make a difference for two
    // particles heading off in the same direction, but i think, in general, it may make a 
    // difference


  return f;  

  // see here for the force equations (especially magnetic):
  // https://physics.stackexchange.com/questions/166318/magnetic-force-between-two-charged-particles
  // http://teacher.nsrl.rochester.edu/phy122/Lecture_Notes/Chapter30/chapter30.html

  // todo: we need some precautions for cases when the p1 and p2 are (almost) at the same position
  // we get division by zero in such cases (because r2 becomes 0)

  // todo: check, if the magnetic force law has the correct sign (later in the stackexchange 
  // discussion, there's a formula that has the cross product in different order)
  // ...we should really check for all formulas, if we compute the forces on p1 due to p2 (as 
  // desired) or the other way around (which is the same value with negative sign)

  // maybe allow for (non-physical) general inverse power force laws instead of the usual inverse
  // square law

  // invent other forces that have no physical counterpart for example based on v1 x v2, p1 x p1,
  // (p1 x v2) x (p2 x v1), p x v, etc. - just take care that f(j,i) = -f(i,j) holds for every
  // such force
}

template<class T>
void rsParticleSystem<T>::updateForces()
{
  size_t N = particles.size();
  size_t i, j;

  for(i = 0; i < N; i++)
    forces[i] = rsVector3D<T>(T(0), T(0), T(0)); // use setComponents(0,0,0) - no copying

  //// naive:
  //for(i = 0; i < N; i++){
  //  for(j = 0; j < N; j++){
  //    if(i != j) forces[i] += getForceBetween(particles[i], particles[j]); }}

  // optimized, using force(j, i) = -force(i, j):
  for(i = 0; i < N; i++){
    for(j = i+1; j < N; j++){
      rsVector3D<T> f = getForceBetween(particles[i], particles[j]);
      forces[i] += f;
      forces[j] -= f; }}

  // todo: test, if both loops give same results (up to roundoff error)
}

template<class T>
void rsParticleSystem<T>::updateVelocities()
{
  for(size_t i = 0; i < particles.size(); i++)
    particles[i].vel += stepSize * forces[i] / particles[i].mass;
}

template<class T>
void rsParticleSystem<T>::updatePositions()
{
  for(size_t i = 0; i < particles.size(); i++)
    particles[i].pos += particles[i].vel;

  // this update equation represents naive, forward Euler integration, i.e.
  // x += dx;
  // maybe try trapezoidal integration instead:
  // x += 0.5*(dx + dxOld);
  // dxOld = dx;
  // ...this formula applies a 2-value moving average filter to the dx signal before accumulating 
  // it into the integrated signal x
  // do the same trapezoidal integration also for velocities
}

template<class T>
void rsParticleSystem<T>::reset()
{
  for(size_t i = 0; i < particles.size(); i++) {
    particles[i].pos = initialPositions[i];
    particles[i].vel = initialVelocities[i]; }
}


/*
Physical Background:

The gravitational law is often expressed as:

F1 = -F2 = G * m1 * m2 * rn / |r|^2  
where: 
F1:        force on particle 1 due to particle 2 (vector)
F2:        force on particle 2 due to particle 1 (vector)
G:         gravitational constant (scalar)
m1, m2:    masses of the two particles (scalars)
p1, p2:    positions of the two particles (vectors)
r = p2-p1: vector from p1 to p2
d := |r|:  Euclidean distance between the two particles (scalar)
rn = r/d:  normalized (unit length) direction vector from p1 to p2 (vector)

It can now also be written as:
F1 = G * m1 * m2 * rn / d^2    | use rn = r/d
   = G * m1 * m2 * r  / d^3 

The nice thing about writing it like the 2nd line is that all of the distance dependency is now
given by the 1/d^3 factor. The problem is when the distance d becomes zero. In the 1st line, we 
would have one division-by-zero when computing rn = r/d and a second one when dividing the whole
thing by d^2. We have isolated the problem into one single 1/d^3 factor. Now, we can tinker with 
the function, like:

F1 = G * m1 * m2 * r * f(d)

when f(d) = 1/d^3, we recover the good old physical law. But we can now try different functions 
like: f(d) = 1 / (c+d^p), f(d) = 1 / (c+d)^p = (1/(c+d))^p, f(d) = ((c+1)/(c+d))^p, etc.
for some parameters c,p. If c = 0, p = 3, we get the old law. With c > 0, we can avoid 
divisions-by-zero and with p != 3, we can have different asymptotic force-dependencies on the 
distance. The asymptotic dependence is given by d * f(d) which reduces to d/d^3 = 1/d^2 for
c = 0, p = 3. With c = 0, p = -1, we get f(d) = d which is like Hooke's law of linear springs. 
Maybe, it's more convenient to use 1 / (c + d^-p) such that we get Hooke's law for p = 1
..oh - wait - is that actually correct?

Questions: 
-Which of the above formulas is best?
-What is a good choice for c?
-Should c depend on the stepSize?
-How does the new force law change the potential energy?

ToDo:
-maybe have different force laws, like Hooke's law F = k*x, where x = d - de
 where de is the equilibrium distance - all particles try to have the same distance

Generally, the force on a particle p1 due to another particle p2 can be expressed by fields like:

F1 = m1*G2 + c1*E2 + c1*(v1 x B2)   | x denotes the cross-product
where:
F1: force on particle 1 (vector)
m1: mass of p1 (scalar)
G2: gravitational field at p1 caused by p2 (vector)
c1: charge of p1 (scalar)
E2: electric field at p1 caused by p2 (vector)
v1: velocity of p1 (vector)
B2: magnetic field at p1 caused by p2 (vector)
We should have functions in rsParticle to compute the gravitational, electric and magnetic
fields at some position vector p. Then, we can compute the forces via the formula above and compare
to the results in rsParticle system in a unit test.
*/

