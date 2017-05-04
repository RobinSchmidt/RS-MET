
//-------------------------------------------------------------------------------------------------
// construction/destruction:

ThreadedDrawingComponent::ThreadedDrawingComponent(TimeSliceThread* newThreadToUse) 
{
  ScopedLock scopedLock(clientAreaImageLock);
  clientAreaImage            = NULL;
  clientAreaImage            = new Image(Image::ARGB, 1, 1, true);
  clientAreaImageIsDirty     = false;
  repaintDelayInMilliseconds = 100;

  if( newThreadToUse != NULL )
    threadToUse = newThreadToUse;
  else
    threadToUse = DrawingThread::getInstance();
    
  threadToUse->addTimeSliceClient(this);
}

ThreadedDrawingComponent::~ThreadedDrawingComponent()
{
  ScopedLock scopedLock(clientAreaImageLock);

  if( threadToUse != NULL )
    threadToUse->removeTimeSliceClient(this);

  if( clientAreaImage != NULL )
    delete clientAreaImage;
}

//-------------------------------------------------------------------------------------------------
// setup:

void ThreadedDrawingComponent::setDrawingThread(TimeSliceThread *newThreadToUse)
{
  if( threadToUse != NULL )
    threadToUse->removeTimeSliceClient(this);

  threadToUse = newThreadToUse;

  if( threadToUse != NULL )
    threadToUse->addTimeSliceClient(this);
}


//-------------------------------------------------------------------------------------------------
// callbacks:

void ThreadedDrawingComponent::timerCallback()
{
  repaint();
  //stopTimer();
}

void ThreadedDrawingComponent::resized()
{
  Component::resized();
  setDirty();
}

void ThreadedDrawingComponent::paint(Graphics &g)
{
	if( threadToUse == NULL )
	{
		Component::paint(g);
		return;
	}
		
  if( clientAreaImageIsDirty )
  {
    //g.fillAll(Colours::white);
    return;
  }

  bool lockAquired = clientAreaImageLock.tryEnter();
  if( lockAquired )  
  {
    if( clientAreaImage != NULL )
      g.drawImageAt(*clientAreaImage, 0, 0, false);
    else
      g.fillAll(Colours::red);

    // could be that we have just drawn a dirty image, so stop the timer only when the image just 
    // drawn is not dirty:
    if( !clientAreaImageIsDirty ) 
      //stopTimer();
      startTimer(1000);

    clientAreaImageLock.exit();
  }
  else
  {
    // we could not aquire the mutex-lock for the image to be drawn - it is probably currently
    // held by the drawing thread - we draw a white background and spawn a new (asynchronous) 
    // call to this function:
    //g.fillAll(Colours::white);
    startTimer(repaintDelayInMilliseconds);
  }
}

int ThreadedDrawingComponent::useTimeSlice()
{
  //ScopedLock scopedLock(clientAreaImageLock);
  if( clientAreaImageIsDirty )
  {
    renderClientAreaImageInternal();
    //return true;             // indicate that the thread is currently very busy - old
  }
  else
  {
    //return true;
    //return false;            // indicate that the thread is currently not too busy
  }
  return 20; // number of milliseconds after this thread shall be called again - maybe we need to tweak
}

//-------------------------------------------------------------------------------------------------
// others:

void ThreadedDrawingComponent::setDirty(bool shouldSetToDirty)
{
  ScopedLock imageLock(clientAreaImageLock);
  clientAreaImageIsDirty = shouldSetToDirty;
  //startTimer(300);

	if( threadToUse == NULL && clientAreaImageIsDirty )
		repaint();
}

void ThreadedDrawingComponent::drawComponent(Image *imageToDrawOnto)
{
  Graphics g(*imageToDrawOnto);
  int w = imageToDrawOnto->getWidth();
  int h = imageToDrawOnto->getHeight();
  g.setColour(Colours::black);
  //g.drawFittedText(String(T("ThreadedDrawingComponent")), 0, 0, w, h, Justification::centred, 1);
	  // triggers a JUCE-breakpoint when called early on app-startup
}

void ThreadedDrawingComponent::renderClientAreaImageInternal()
{
  ScopedLock imageLock(clientAreaImageLock);

  if( getWidth() < 1 || getHeight() < 1 )
  {
    clientAreaImageIsDirty = false;
    return;
  }

  if(  clientAreaImage->getWidth()  != getWidth() 
    || clientAreaImage->getHeight() != getHeight() )
  {
    allocateClientAreaImage(getWidth(), getHeight());
  }

  // initialize the image as a white canvas and set the dirty flag:
  Graphics g(*clientAreaImage);
  //g.fillAll(Colours::white);
  //clientAreaImageIsDirty = true;

  // call the actual drawing rutine (hich is supposed to be overriden by subclasses):
  drawComponent(clientAreaImage);
 
  // when the drawComponent-function returns, we assume that the clientArwaImage has been drawn, so
  // we set our dirty flag flase and trigger a (delayed) repaint:
  clientAreaImageIsDirty = false;
  startTimer(repaintDelayInMilliseconds);
}

bool ThreadedDrawingComponent::allocateClientAreaImage(int desiredWidth, int desiredHeight)
{
  ScopedLock imageLock(clientAreaImageLock);

  bool result = false;

  desiredWidth  = jmax(1, desiredWidth);
  desiredHeight = jmax(1, desiredHeight);

  // allocate memory for the first time:
  if( clientAreaImage == NULL )
  {
    clientAreaImage = new Image(Image::ARGB, desiredWidth, desiredHeight, true);

    if( clientAreaImage == NULL )
    {
      jassertfalse;
      //showMemoryAllocationErrorBox(String("ThreadedDrawingComponent::allocateClientAreaImage")); // re-activate
      return false;  
    }
    else
      result = true;
  }

  // reallocate memory, if necessary (i.e. the deired size differs from the current size of the 
  // image):
  if(    clientAreaImage->getWidth()  != desiredWidth  
      || clientAreaImage->getHeight() != desiredHeight  )
  {
    // delete the old and create a new Image-object:
    if( clientAreaImage != NULL )
    {
      delete clientAreaImage;
      clientAreaImage = NULL;
    }
    clientAreaImage = new Image(Image::ARGB, desiredWidth, desiredHeight, true);
    result = true;

    if( clientAreaImage == NULL )
    {
      //showMemoryAllocationErrorBox(String("ThreadedDrawingComponent::allocateClientAreaImage")); re-activate
      jassertfalse;
      return false; 
    }
  }

  // when we indeed have allocated new memory, the image associated with this new memory is 
  // certainly not what we want to see:
  if( result = true )
  {
    clientAreaImageIsDirty = true;
    Graphics g(*clientAreaImage);
    g.fillAll(Colours::white);
  }

  return result;
}



