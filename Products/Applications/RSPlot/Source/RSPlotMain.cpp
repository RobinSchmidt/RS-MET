//#include "../../../../third_party_code/juce/juce.h"
//#include "../../../rojue/misc/MetallicLookAndFeel.h"

#include "RSPlotContentComponent.h"

//=============================================================================
/** 

This is the top-level window that we'll pop up. Inside it, we'll create and
show a RSPlotContentComponent component.

*/

class RSPlotWindow : public DocumentWindow
{

public:

  RSPlotWindow() : DocumentWindow(("RSPlot"),
    Colours::lightgrey, 
    DocumentWindow::allButtons, 
    true)
  {
    coordinateSystemDemoContentComponent = 
      new RSPlotContentComponent(String(("RSPlot")));
    setContentOwned(coordinateSystemDemoContentComponent, false);
    //setContentComponent(coordinateSystemDemoContentComponent); // deprecated

    setResizable(true, true);
    setResizeLimits(800, 600, 40000, 30000);
    setVisible(true);

    centreWithSize(800, 600);

    setTitleBarHeight(24);
  }

  ~RSPlotWindow()
  {
    // the content component will be deleted automatically, so no need to do
    // it here
  }

  void closeButtonPressed()
  {
    // When the user presses the close button, we'll tell the app to quit. This 
    // window will be deleted by the app object as it closes down.
    JUCEApplication::quit();
  }

  void resized()
  {
    DocumentWindow::resized();
    coordinateSystemDemoContentComponent->setBounds(
      0, getTitleBarHeight(), getWidth(), getHeight()-getTitleBarHeight());
    //contentComponent->setSize(getWidth(), getHeight());
  }

  void setBounds(int x, int y, int width, int height)
  {
    DocumentWindow::setBounds(x, y, width, height);
    coordinateSystemDemoContentComponent->setBounds(
      0, getTitleBarHeight(), getWidth(), getHeight()-getTitleBarHeight());
    //contentComponent->setSize(getWidth(), getHeight());
  }


protected:

  RSPlotContentComponent* coordinateSystemDemoContentComponent;

};


//=============================================================================
/** 
This is the application object that is started up when Juce starts. It handles
the initialisation and shutdown of the whole application.
*/

class RSPlot : public JUCEApplication
{

  /* Important! NEVER embed objects directly inside your JUCEApplication class! 
  Use ONLY pointers to objects, which you should create during the 
  initialise() method (NOT in the constructor!) and delete in the shutdown() 
  method (NOT in the destructor!)

  This is because the application object gets created before Juce has been 
  properly initialised, so any embedded objects would also get constructed 
  too soon.
  */

  RSPlotWindow* coordinateSystemDemoWindow;

  //RLookAndFeel rLookAndFeel;
  //LookAndFeel lookAndFeel;
  //MetallicLookAndFeel lookAndFeel;

public:

  RSPlot() : coordinateSystemDemoWindow (0)
  {
    // NEVER do anything in here that could involve any Juce function being 
    // called - leave all your startup tasks until the initialise() method.
  }

  ~RSPlot()
  {
    // Your shutdown() method should already have done all the things necessary to
    // clean up this app object, so you should never need to put anything in
    // the destructor.

    // Making any Juce calls in here could be very dangerous...
  }

  void initialise (const String& commandLine)
  {
    //LookAndFeel::setDefaultLookAndFeel(RLookAndFeel::getInstance());

    // just create the main window...
    coordinateSystemDemoWindow = new RSPlotWindow();

    /*  ..and now return, which will fall into to the main event
    dispatch loop, and this will run until something calls
    JUCEAppliction::quit().

    In this case, JUCEAppliction::quit() will be called by the
    hello world window being clicked.
    */
  }

  void shutdown()
  {
    // clear up:
    if (coordinateSystemDemoWindow != 0)
      delete coordinateSystemDemoWindow;
  }

  const String getApplicationName()
  {
    return ("RSPlot");
    //return T("CoordinateSystem3D Demo");
  }

  const String getApplicationVersion()
  {
    return ("1.0");
  }

  bool moreThanOneInstanceAllowed()
  {
    return true;
  }

  void anotherInstanceStarted (const String& commandLine)
  {

  }

  /*
  virtual void unhandledException (const std::exception* e,
    const String& sourceFilename,
    const int lineNumber)
  {
    // this function was overriden by Robin Schmidt as a workaround to avoid
    // a myterious linker error in Visual C++ 2005 Express
  }
  */

};


//==============================================================================
// This macro creates the application's main() function..
START_JUCE_APPLICATION (RSPlot)
