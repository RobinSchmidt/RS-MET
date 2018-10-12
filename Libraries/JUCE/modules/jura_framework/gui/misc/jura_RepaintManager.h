#ifndef jura_RepaintManager_h
#define jura_RepaintManager_h


class rsRepaintClient
{
public:

  virtual bool needsRepaint() = 0;

};
// hmm...maybe this class is not needed...we'll see

//=================================================================================================

/*
Idea:
-RepaintManager derives from juce::Timer
-in the timer callback we spawn calls to repaint() (on the message thread) to all registered 
 repaintees
-repaintees are juce::Components and register via registerRepaintee
-functions:
 setRepaintInterval/Rate
-maybe have some mechanism to figure out, if a repaint is actually needed or if the gui is up to 
 date ...maybe via a class RepaintClient (that derives from juce::Component) with a virtual function
 needsRepaint - when the repinters spaws the repaints, it tries to casts the repaintees to 
 RepaintClient and if this it successful, it checks needsRepaint and repaints only, if this returns 
 true - if the cast fails, it will just spawn repaint
-Modulatab
 


*/

class JUCE_API rsRepaintManager : public juce::Timer
{

public:

  rsRepaintManager();

  void registerRepaintee(juce::Component* repaintee)
  {
    RAPT::rsAppendIfNotAlreadyThere(repaintees, repaintee);
  }

  void deRegisterRepaintee(juce::Component* repaintee)
  {
    RAPT::rsRemoveFirstOccurrence(repaintees, repaintee);
  }

  void setRepaintRate(int newRateInHz)
  {
    if(newRateInHz <= 0)
      stopTimer();
    else
      startTimerHz(newRateInHz);
  }

  void timerCallback() override;


protected:

  std::vector<juce::Component*> repaintees;

};

#endif  