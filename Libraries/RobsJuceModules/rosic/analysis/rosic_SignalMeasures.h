#ifndef rosic_SignalMeasures_h
#define rosic_SignalMeasures_h

namespace rosic
{

  /**

  This is a class which wraps certain measurements of a signal.

  */

  class SignalMeasures
  {
  public:

    SignalMeasures();  /**< Constructor. */

    double leftLevel;
    double rightLevel;
    double midLevel;
    double sideLevel;
    double crossCorrelation;
  };


} // end namespace rosic

#endif // rosic_SignalMeasures_h
