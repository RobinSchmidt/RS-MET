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
  return E;
}

template<class T>
T rsParticleSystem<T>::getPotentialEnergy()
{
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
      // W = F * s
    }
  }
  return E;
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
rsVector3D<T> rsParticleSystem<T>::getForceBetween(const rsParticle<T>& p1, const rsParticle<T>& p2)
{
  T k = 0;      // for physical force-law with singularity (which spoils numeric simulation)
  //k = 20*stepSize; // test - make this a user parameter
  k = 1.0f;
  T p = 2.0f; // freq seems to be (roughly) inversely proportional to this

  // instead of the physical inverse square force-law F = k / r^2, we may use F = k / (c+r)^p which
  // reduces to the physical law for c=0,p=2 - allows to mitigate sigularity effects and gives
  // more flexibility, maybe c should depend on the stepSize and/or exponent? try to figure 
  // something out that makes the behavior more or less independent from the stepSize

  // precomputations:
  rsVector3D<T> r = p2.pos - p1.pos;    // vector pointing from p1 to p2

  //// old - physically correct for k=0:
  //T r2  = r.getSquaredEuclideanNorm();  // squared distance between p1 and p2 == |r|^2
  //T s   = 1 / (k + r2);                       // reciprocal of |r|^2 - used as multiplier in various places
  //r *= sqrt(s);                         // r is now normalized to unit length (for k=0)

  // new:
  T d = r.getEuclideanNorm();      // distance between p1 and p2 == |r|
  //T s = 1 / (k + pow(d,p));        // reciprocal of k+|r|^p - used as multiplier in various places
  T s = 1 / pow(k+d,p);            // maybe 1 / pow(k+d, p) is better - experiment
  r *= 1/(k+d);                    // r is now normalized to unit length (for k=0)


  // compute the 3 forces:
  rsVector3D<T> f;
  f += s*cG * p1.mass   * p2.mass   * r;                                // gravitational force
  f -= s*cE * p1.charge * p2.charge * r;                                // electric force
  f += s*cM * p1.charge * p2.charge * cross(p1.vel, cross(p2.vel, r));  // magnetic force
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
}

template<class T>
void rsParticleSystem<T>::updateForces()
{
  size_t N = particles.size();
  size_t i, j;

  for(i = 0; i < N; i++)
    forces[i] = 0;

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
    particles[i].vel += forces[i] / particles[i].mass;
}

template<class T>
void rsParticleSystem<T>::updatePositions()
{
  for(size_t i = 0; i < particles.size(); i++)
    particles[i].pos += stepSize * particles[i].vel;
}

template<class T>
void rsParticleSystem<T>::reset()
{
  for(size_t i = 0; i < particles.size(); i++) {
    particles[i].pos = initialPositions[i];
    particles[i].vel = initialVelocities[i]; }
}