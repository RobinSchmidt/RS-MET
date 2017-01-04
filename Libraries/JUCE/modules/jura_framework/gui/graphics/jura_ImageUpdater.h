#ifndef jura_ImageUpdater_h
#define jura_ImageUpdater_h

/** A class for keeping track of updates of an juce::Image object. Works in conjunction with the 
ImageUpdater class. The idea is that you may update an image inside some worker thread and the let 
your GUI component be repainted whenever that image was updated. */

class JUCE_API ImageUpdateListener
{
public:
  ImageUpdateListener() {}
  virtual void imageWasUpdated(Image* image) = 0;
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ImageUpdateListener)
};


/** A class for keeping objects informed about updates of an Image object. */

class JUCE_API ImageUpdater
{

public:

  ImageUpdater() {}
  virtual ~ImageUpdater(){}

  virtual void addImageUpdateListener(ImageUpdateListener* listenerToAdd);
  virtual void removeImageUpdateListener(ImageUpdateListener* listenerToRemove);
  virtual void sendImageUpdateNotification(Image* image);

protected:

  Array<ImageUpdateListener*> imageUpdateListeners;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ImageUpdater)
};

#endif  
