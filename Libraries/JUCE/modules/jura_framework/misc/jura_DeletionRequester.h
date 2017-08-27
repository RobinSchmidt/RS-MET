#ifndef jura_DeletionRequester_h
#define jura_DeletionRequester_h

class rsDeletor;

/** A baseclass for objects that must ultimately request their own deletion from some other object
which can be considered to be the owner of the to-be-deleted object. */

class rsDeletionRequester
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

class rsDeletor
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

  virtual void deleteObject(rsDeletionRequester* objectToDelete)
  {
    delete objectToDelete;
  }


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsDeletor)
};

// todo: Maybe make variations of the Deletor - in particular, a GarbageCollector that doesn't 
// delete immediately but only marks an object for later deletion and receives a (Timer?) callback 
// at regular intervals inside of which the actual deletions occur. Maybe rename rsDeletor to 
// rsImmediateDeletor...or keep rsDeletor as baseclass for rsImmediateDeletor, 
// rsImmediateDeletorSafe, rsGarbageCollector, etc.. rename this file-pair to 
// DeferredDeletion.h/cpp

#endif  
