#ifndef jura_ObservablePlot_h
#define jura_ObservablePlot_h

#if defined _MSC_VER
#pragma warning(disable : 4250) // disable dominance warning in MSVC
#endif

// ToDo: make this class obsolete by having an rsPlotSettingsObserver

class rsObservablePlot;

/** This class is an observer for rsObservablePlot objects. Subclasses of this class have to 
override the callback function coordinateSystemChanged in order to respond to any changes in the 
observed rsObservablePlot object(s). */

class rsPlotObserver
{

public:

  virtual ~rsPlotObserver() {}

  /** This is the callback, you must override in order to respond to any changes in the observed 
  rsObservablePlot object(s). */
  virtual void coordinateSystemChanged(rsObservablePlot *coordinateSystemThatHasChanged) = 0;

  juce_UseDebuggingNewOperator;

};


//=================================================================================================

/**

This class extends the rsPlot class in auch a way as to broadcast change-messages, whenever a parameter (such as the visible
range) changes. It is intended to be used for plugIns and applications that must be able to restore the appearance of the
rsPlot. If you derive your class from this virtual subclass and another (virtual) subclass of rsPlot via
multiple inheritance, make sure, that this class comes first - otherwise the attached ChangeListeners will get a pointer with a wrong
address in their changeListenerCallback functions.

Example:
-good: class SpectrumDisplay	: public rsObservablePlot, public CurveFamilyPlot
-bad:  class SpectrumDisplay	: public CurveFamilyPlot, public rsObservablePlot


\todo: make some dedicated class CoordinateSystemOldListener or something - this callback-thing with the pointer to void as argument does
not seem to work reliably (inherited sub-objects may have addresses other than expected)

\todo
-make a class rsPlotRangeObserver - it's better to observe the range itself, we can get rid of this
 virtual inheritance stuff then (which is a pita)

*/

class rsObservablePlot : virtual public rsPlot
{

public:

  //-------------------------------------------------------------------------------------------------------------------------------------
  // construction/destruction:

  rsObservablePlot(const juce::String &name = juce::String("MessangingCoordinateSytem"));
  virtual ~rsObservablePlot();

  //-------------------------------------------------------------------------------------------------------------------------------------
  // range-management:

  virtual void setMaximumRange(double newMinX, double newMaxX, double newMinY, double newMaxY);
  virtual void setMaximumRange(rsPlotRange newMaximumRange);
  virtual void setMaximumRangeX(double newMinX, double newMaxX);
  virtual void setMaximumRangeY(double newMinY, double newMaxY);
  virtual void setMaximumRangeMinX(double newMinX);
  virtual void setMaximumRangeMaxX(double newMaxX);
  virtual void setMaximumRangeMinY(double newMinY);
  virtual void setMaximumRangeMaxY(double newMaxY);

  virtual void setCurrentRange(double newMinX, double newMaxX, double newMinY, double newMaxY);
  virtual void setCurrentRange(rsPlotRange newRange);
  virtual void setCurrentRangeX(double newMinX, double newMaxX);
  virtual void setCurrentRangeY(double newMinY, double newMaxY);
  virtual void setCurrentRangeMinX(double newMinX);
  virtual void setCurrentRangeMaxX(double newMaxX);
  virtual void setCurrentRangeMinY(double newMinY);
  virtual void setCurrentRangeMaxY(double newMaxY);


  //-------------------------------------------------------------------------------------------------------------------------------------
  // others:

  /** Registers an observer which will get notified about changes. */
  virtual void addCoordinateSystemOldObserver(rsPlotObserver* observerToAdd);

  /** De-registers an observer. */
  virtual void removeCoordinateSystemOldObserver(rsPlotObserver* observerToRemove);

  /** De-registers all observers. */
  virtual void removeAllCoordinateSystemOldObservers();

  /** Calls rsPlotObserver::coordinateSystemChanged for all our registered observers. */
  virtual void sendCoordinateSystemChangedMessage(rsObservablePlot *coordinateSystemThatHasChanged);


  juce_UseDebuggingNewOperator;

protected:

  juce::Array<rsPlotObserver*, CriticalSection> observers;

};

#endif
