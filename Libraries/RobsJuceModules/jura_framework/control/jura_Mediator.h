#ifndef jura_Mediator_h
#define jura_Mediator_h

// todo: maybe move to rosic, prepend rs prefix

class Mediator;

/** This class represents an object that can communicate with other objects through a mediator 
as in the Mediator pattern described in GoF - Design Patterns. */

class MediatedColleague
{

public:

  friend class Mediator;

  //-----------------------------------------------------------------------------------------------
  // \name Construction/Destruction:

  /** Constructor. */
  MediatedColleague() {}

  /** Destructor. */
  virtual ~MediatedColleague();

  //-----------------------------------------------------------------------------------------------
  // \name Setup/Inquiry/Mediation:

  /** Sets the mediator to be used by this colleague. The function will also take care to register 
  this colleague with the new mediator and perhaps de-register from an old mediator (if any). */
  virtual void setMediator(Mediator *newMediator);

  /** Returns the mediator that is used by this colleague. */
  virtual Mediator* getMediator() const { return mediator; }

  /** Notifies the mediator (if any) by calling its 
  colleagueHasSentNotification(MediatedColleague*, int) function with "this" as first argument. The
  second argument is assigned from the functions parameter (i.e. passed through). */
  virtual void notifyMediator(int messageCode = 0, void* messageData = nullptr);

  /** This is the callback function that you must override in your concrete MediatedColleague 
  subclass in order to respond to notifications from the mediator. The mediator spawns this 
  callback on all colleagues (including the colleague that originated the call) whenever one of the 
  colleagues notifies the mediator via notifyMediator. The mediator passes the originating 
  colleague through, together with some optional message code which can be used to further dispatch 
  on some message code and/or with some subclass-defined message data.  */
  virtual void handleMediatorNotification(MediatedColleague *originatingColleague, 
    int messageCode = 0, void* messageData = nullptr) = 0;
  // ToDo:
  // -Maybe pass a void pointer to a data-object. Sometimes, an integer message-code is not enough
  //  and we need some more data...maybe it should not be a void pointer but a pointer to some 
  //  rsMessageData class - this is needed because dynamic_cast doesn't work for void-pointers. We
  //  could do a static cast, but that may be a bit dangerous - we should perhaps be amore 
  //  defensive


protected:

  Mediator *mediator = nullptr;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MediatedColleague)
};

//=================================================================================================

/** This class represents an object that mediates a many-to-many communication between a bunch of 
participating colleagues as described in GoF - Design Patterns. It's a variant of the pattern that
uses message codes to identify the messages in order to be as general as possible. The concrete 
subclasses must define the meanings of the message codes. */

class Mediator
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Construction/Destruction:

  /** Constructor. */
  Mediator();

  /** Destructor. */
  virtual ~Mediator();

  //-----------------------------------------------------------------------------------------------
  // \name Mediation management:

  /** Registers a colleague that wants to get notified whenever one of the participating colleagues 
  spawns a message. */
  virtual void registerColleague(MediatedColleague *colleagueToRegister);

  // De-registers a previously registered colleague. */
  virtual void deRegisterColleague(MediatedColleague *colleagueDeToRegister);

  /** This will call mediatorHasSentNotification(MediatedColleague*, int messageCode) on all 
  participating colleagues with parameters passed through from the function's arguments. */
  virtual void sendNotificationToColleagues(
    MediatedColleague *originatingColleague, int messageCode = 0, void* messageData = nullptr);

  /** This is the callback function that gets called whenever one of the colleagues notifies "this" 
  mediator. The baseclass implementation will just call sendNotificationToColleagues with the 
  parameters passed through. You may want override it in your Mediator subclass in order to 
  implement more specific behavior. */
  virtual void colleagueHasSentNotification(MediatedColleague *originatingColleague, 
    int messageCode = 0, void* messageData = nullptr);


protected:

  std::vector<MediatedColleague*> colleagues;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Mediator)
};

#endif  
