#ifndef jura_RAudioProcessorEditor_h
#define jura_RAudioProcessorEditor_h

//#include "RPlugInEngineEditor.h"
////#include "RAudioProcessor.h"
////class RAudioProcessor;
////class RAudioModuleEditor;
////#include "../../../third_party_code/juce/juce.h"
////#include "RAudioProcessorEditor.h"
////#include "../../../rojue/components/RobsSlider.h"
////#include "../../../rojue/misc/GraphicsTools.h"
////#include "../../../rojue/misc/MetallicLookAndFeel.h"

class JUCE_API RAudioProcessorEditor : public AudioProcessorEditor, public ChangeListener
{

public:

  /** Constructor. When created, this will register itself with the filter for changes. It's safe 
  to assume that the filter won't be deleted before this object is. */
  RAudioProcessorEditor(RAudioProcessor* const ownerPlugIn);

  /** Destructor. */
  ~RAudioProcessorEditor();

  /**<Implements the purely virtual changeListenerCallback()-method of the ChangeListener 
  base-class. */
  virtual void changeListenerCallback(void *objectThatHasChanged);

  /** Standard Juce paint callback. */
  virtual void paint(Graphics& g);


private:

  // this object just wraps the actual editor-object to interface it with the actual 
  // AudioProcessor, so we have the actual editor as member:
  // RPlugInEngineEditor* audioEngineEditor;

  void updateParametersFromFilter();

  // handy wrapper method to avoid having to cast the filter to a RAudioProcessor
  // every time we need it..
  RAudioProcessor* getFilter() const throw()       
  { return (RAudioProcessor*) getAudioProcessor(); }

};

#endif
