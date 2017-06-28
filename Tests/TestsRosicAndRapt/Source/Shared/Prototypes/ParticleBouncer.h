#ifndef ParticleBouncer_h
#define ParticleBouncer_h

/**   */

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
