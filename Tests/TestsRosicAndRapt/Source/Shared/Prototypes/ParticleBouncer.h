#ifndef ParticleBouncer_h
#define ParticleBouncer_h

/** A sound generator based on the idea of a particle that bounces around inside an enclosure. */

class ParticleBouncer
{

public:

 //------------------------------------------------------------------------------------------------
 // construction/destruction:

  ParticleBouncer();   ///< Constructor.

  ~ParticleBouncer();  ///< Destructor.

 //------------------------------------------------------------------------------------------------
 // audio processing:

  double getSample();

  //-----------------------------------------------------------------------------------------------
  // others:

  void reset();

protected:

};

double ParticleBouncer::getSample()
{
  return 0.0;
}

#endif
