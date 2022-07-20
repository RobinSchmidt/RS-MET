#ifndef rosic_SignalMeasures_h
#define rosic_SignalMeasures_h

namespace rosic
{

  /** This is a class which wraps certain measurements of a signal.

  ToDo: move this code into the TrackMeter class, make it a struct.  */

  class SignalMeasures
  {
  public:

    SignalMeasures();  /**< Constructor. */

    float leftLevel;
    float rightLevel;
    float midLevel;
    float sideLevel;
    float crossCorrelation;
  };




} // end namespace rosic

#endif // rosic_SignalMeasures_h
