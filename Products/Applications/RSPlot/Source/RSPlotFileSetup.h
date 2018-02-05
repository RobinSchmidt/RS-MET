#ifndef __JUCE_RSPLOTFILESETUP_JUCEHEADER__
#define __JUCE_RSPLOTFILESETUP_JUCEHEADER__

//#include "../../../rojue/misc/StateFileManager.h"
//#include "../../../rojue/components/SymbolButton.h"
//#include "../../../rojue/components/widgets/rojue_RButton.h"
//#include "../../../rojue/components/widgets/rojue_RLabel.h"
//#include "../../../rojue/components/widgets/rojue_RTextField.h"
//using namespace rojue;

#include "../JuceLibraryCode/JuceHeader.h"
using namespace juce;
using namespace jura;

class RSPlotFileSetup	: public Component
{

public:

  RSPlotFileSetup();
  /**< Constructor. */

  virtual ~RSPlotFileSetup();
  /**< Destructor. */

  virtual void paint(Graphics &g) override;
  virtual void resized() override;


  //=============================================================================================
  juce_UseDebuggingNewOperator;


  RButton    *loadButton, *saveButton, *plusButton, *minusButton, *pngExportButton, 
             *svgExportButton;

  RTextField *imageExportLabel, *imageWidthLabel, *imageHeightLabel;

  RTextEntryField *fileNameLabel, *imageWidthEditLabel, *imageHeightEditLabel;

  // old:
  //ToggleButton*  loadButton;
  //ToggleButton*  saveButton;
  //RButton*       plusButton;
  //RButton*       minusButton;
  //ToggleButton*  pngExportButton;
  //ToggleButton*  svgExportButton;

  //Label*         imageExportLabel;
  //Label*         imageWidthLabel;
  //Label*         imageHeightLabel;

  //Label*         fileNameLabel;
  //Label*         imageWidthEditLabel;
  //Label*         imageHeightEditLabel;
};

#endif	
