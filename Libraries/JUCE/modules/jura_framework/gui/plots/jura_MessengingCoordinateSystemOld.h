#ifndef jura_MessengingCoordinateSystemOld_h
#define jura_MessengingCoordinateSystemOld_h

#pragma warning(disable : 4250) // disable dominance warning in MSVC



class MessengingCoordinateSystemOld;

/** This class is an observer for MessengingCoordinateSystemOld objects. Subclasses of this class have to override the callback function
coordinateSystemChanged in order to respond to any changes in the observed MessengingCoordinateSystemOld object(s). */

class CoordinateSystemOldObserver
{

public:

  virtual ~CoordinateSystemOldObserver() {}

  /** This is the callback, you must override in order to respond to any changes in the observed MessengingCoordinateSystemOld
  object(s). */
  virtual void coordinateSystemChanged(MessengingCoordinateSystemOld *coordinateSystemThatHasChanged) = 0;

  juce_UseDebuggingNewOperator;

};


//=======================================================================================================================================

/**

This class extends the CoordinateSystemOld class in auch a way as to broadcast change-messages, whenever a parameter (such as the visible
range) changes. It is intended to be used for plugIns and applications that must be able to restore the appearance of the
CoordinateSystemOld. If you derive your class from this virtual subclass and another (virtual) subclass of CoordinateSystemOld via
multiple inheritance, make sure, that this class comes first - otherwise the attached ChangeListeners will get a pointer with a wrong
address in their changeListenerCallback functions.

Example:
-good: class SpectrumDisplay	: public MessengingCoordinateSystemOld, public CurveFamilyPlot
-bad:  class SpectrumDisplay	: public CurveFamilyPlot, public MessengingCoordinateSystemOld


\todo: make some dedicated class CoordinateSystemOldListener or something - this callback-thing with the pointer to void as argument does
not seem to work reliably (inherited sub-objects may have addresses other than expected)

*/

class MessengingCoordinateSystemOld : virtual public CoordinateSystemOld
{

public:

  //-------------------------------------------------------------------------------------------------------------------------------------
  // construction/destruction:

  MessengingCoordinateSystemOld(const juce::String &name = juce::String("MessangingCoordinateSytem"));
  virtual ~MessengingCoordinateSystemOld();

  //-------------------------------------------------------------------------------------------------------------------------------------
  // range-management:

  virtual void setMaximumRange(double newMinX, double newMaxX, double newMinY, double newMaxY);
  virtual void setMaximumRange(CoordinateSystemRangeOld newMaximumRange);
  virtual void setMaximumRangeX(double newMinX, double newMaxX);
  virtual void setMaximumRangeY(double newMinY, double newMaxY);
  virtual void setMaximumRangeMinX(double newMinX);
  virtual void setMaximumRangeMaxX(double newMaxX);
  virtual void setMaximumRangeMinY(double newMinY);
  virtual void setMaximumRangeMaxY(double newMaxY);

  virtual void setCurrentRange(double newMinX, double newMaxX, double newMinY, double newMaxY);
  virtual void setCurrentRange(CoordinateSystemRangeOld newRange);
  virtual void setCurrentRangeX(double newMinX, double newMaxX);
  virtual void setCurrentRangeY(double newMinY, double newMaxY);
  virtual void setCurrentRangeMinX(double newMinX);
  virtual void setCurrentRangeMaxX(double newMaxX);
  virtual void setCurrentRangeMinY(double newMinY);
  virtual void setCurrentRangeMaxY(double newMaxY);


  //-------------------------------------------------------------------------------------------------------------------------------------
  // others:

  /** Registers an observer which will get notified about changes. */
  virtual void addCoordinateSystemOldObserver(CoordinateSystemOldObserver* observerToAdd);

  /** De-registers an observer. */
  virtual void removeCoordinateSystemOldObserver(CoordinateSystemOldObserver* observerToRemove);

  /** De-registers all observers. */
  virtual void removeAllCoordinateSystemOldObservers();

  /** Calls CoordinateSystemOldObserver::coordinateSystemChanged for all our registered observers. */
  virtual void sendCoordinateSystemChangedMessage(MessengingCoordinateSystemOld *coordinateSystemThatHasChanged);


  juce_UseDebuggingNewOperator;

protected:

  juce::Array<CoordinateSystemOldObserver*, CriticalSection> observers;

};

#endif
