#ifndef jura_DeletionRequester_h
#define jura_DeletionRequester_h

class rsDeletor;

/** A baseclass for objects that must ultimately request their own deletion from some other object
which can be considered to be the owner of the to-be-deleted object. */

class JUCE_API rsDeletionRequester
{

public:

  /** Constructor. You should pass a valid pointer to an rsDeletorObject which will be used to 
  issue the deletion request for this object. */
  rsDeletionRequester(rsDeletor* deletorToUse) : deletor(deletorToUse) {}

  /** Destructor. */
  virtual ~rsDeletionRequester() {}

  /** Subclasses can call this function to request their own deletion. */
  void requestDeletion();

protected:

  rsDeletor *deletor;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsDeletionRequester)
};

//=================================================================================================

/** A baseclass for objects that can delete other objects (of a subclass of rsDeletionRequester) on
their request. */

class JUCE_API rsDeletor
{

public:

  /** Constructor. */
  rsDeletor() {}

  /** Destructor. */
  virtual ~rsDeletor() 
  {
    // Maybe we should also keep pointers to all our rsDeletionRequesters and delete them in our 
    // destructor. Otherwise, the rsDeletor may be deleted before the referencing 
    // rsDeletionRequesters are being deleted leading to dangling pointers there. On the other 
    // hand, that creates overhead in cases where it can be assured that the deletor outlives its
    // deletees. maybe it should be done in a subclass rsSafeDeletor or something?
  }

  /** Deletes the passed object (immediately). */
  virtual void deleteObject(rsDeletionRequester* objectToDelete)
  {
    delete objectToDelete;
  }


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsDeletor)
};

// todo: Maybe rename this file-pair to DeferredDeletion.h/cpp

//=================================================================================================

/** A subclass of rsDeletor that doesn't immediately delete the object that is passed but instead
moves it into a garbage collector that will get emptied at some later time (using a timer 
callback). Deferring the deletion to some later time is sometimes necessary when an object wants to 
delete itself, because after being marked for deletion, it may still receive callbacks. For 
example, a button that deletes itself (immediately) when clicking on it (on mouseDown), would still 
receive a mouseUp event after that deletion, leading to an access violation. Think of a close 
button that is part of a window. Clicking on may delete the window and with it the button. Such 
things can be circumvented by deferring the deletion (although, in this case, it may also be an
option to delet on mouseUp). */

class JUCE_API rsGarbageCollector : public rsDeletor, public juce::Timer
{

public:

  /** Constructor. You can pass a time-interval in milliseconds at which cleanups (i.e. actual 
  deletion of the garbage) will occur. */
  rsGarbageCollector(int cleanUpInterval = 1000);

  /** Destructor. */
  virtual ~rsGarbageCollector();

  /** Sets the interval (in milliseconds) at which the garbage will actually be deleted. If there 
  is no garbage to delete, the timerCallback will actually not get called (the timer is started 
  when garbage is added and stopped after a cleanup). */
  void setCleanUpInterval(int newInterval);

  /** Moves the object to which the passed pointer points into our garbage bin for later 
  deletion. */
  void deleteObject(rsDeletionRequester* objectToDelete) override;

  /** Overrides the timerCallback in order to epmty our trash can. */
  void timerCallback() override;


protected:

  /** Actually deletes the garbage in our array. */
  void deleteGarbage();

  std::vector<rsDeletionRequester*> garbage; // array of pointers to garbage objects
  int cleanUpInterval = 1000;                // in milliseconds

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsGarbageCollector)
};


#endif  
