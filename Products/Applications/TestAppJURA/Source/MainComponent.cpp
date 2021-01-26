#ifndef MAINCOMPONENT_H_INCLUDED  // hey, but it's a .cpp file (this was auto generated)
#define MAINCOMPONENT_H_INCLUDED

#include  "Tests/RAPT/Visualization/ImageProcessing/ImagePainterTests.h"
#include  "Tests/UnitTestsView.h"


/** This component lives inside our window, and this is where you should put all your controls 
and content. */

class MainContentComponent : public AudioAppComponent
{

public:

  MainContentComponent()
  {
    initSubComponents();
    createWidgets();
    setSize(800, 600);
    setAudioChannels(2, 2); // number of input and output channels
  }

  ~MainContentComponent()
  {
    // try to get rid of these delet calls:
    delete buttonPainter;
    delete buttonUnitTests;
    shutdownAudio();
  }


  // Audio callbacks:

  void prepareToPlay (int samplesPerBlockExpected, double sampleRate) override
  {
    // This function will be called when the audio device is started, or when
    // its settings (i.e. sample rate, block size, etc) are changed.

    // You can use this function to initialise any resources you might need,
    // but be careful - it will be called on the audio thread, not the GUI thread.
    
    // For more details, see the help for AudioProcessor::prepareToPlay()
  }

  void getNextAudioBlock(const AudioSourceChannelInfo& bufferToFill) override
  {
    // see AudioProcessor::getNextAudioBlock()
    bufferToFill.clearActiveBufferRegion();
  }

  void releaseResources() override
  {
    // This will be called when the audio device stops, or when it is being
    // restarted due to a setting change.
 
    // For more details, see the help for AudioProcessor::releaseResources()
  }


  // Overriden Component callbacks:

  //void paint (Graphics& g) override
  //{
  //    // (Our component is opaque, so we must completely fill the background with a solid colour)
  //  g.fillAll (Colours::black);


  //  // You can add your drawing code here!
  //}


  void resized() override
  {
    int m  = 4;   // margin
    int x  = m;
    int y  = m;
    int w  = 80; // button width
    int h  = 20; // button height
    int dx = w + m;

    // tab buttons:
    buttonPainter->setBounds(  x, y, w, h); x += dx;
    buttonUnitTests->setBounds(x, y, w, h);

    // sub components:
    x = 0;
    y = buttonUnitTests->getBottom() + m;
    w = getWidth()  - x;
    h = getHeight() - y;
    painter.setBounds(      x, y, w, h);
    unitTestsView.setBounds(x, y, w, h);
  }

private:

  void createWidgets()
  {
    buttonPainter = new RRadioButton("Painter");
    buttonPainter->addToRadioButtonGroup(&tabButtonGroup);
    addAndMakeVisible(buttonPainter); // maybe use addWidget (inherit from jura::Editor)

    buttonUnitTests = new RRadioButton("Unit Tests");
    buttonUnitTests->addToRadioButtonGroup(&tabButtonGroup);
    addAndMakeVisible(buttonUnitTests);
  }

  void initSubComponents()
  {
    //addAndMakeVisible(painter);
    addAndMakeVisible(unitTestsView);
  }
  // todo: implement response to buttons to switch between the screens



  // widgets:
  RRadioButtonGroup tabButtonGroup;
  RRadioButton *buttonPainter, *buttonUnitTests;

  // sub-components
  PainterComponent painter;
  UnitTestsView unitTestsView;


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MainContentComponent)
};

// (This function is called by the app startup code to create our main component)
Component* createMainContentComponent() 
{ 
  return new MainContentComponent(); 
}

#endif  // MAINCOMPONENT_H_INCLUDED
