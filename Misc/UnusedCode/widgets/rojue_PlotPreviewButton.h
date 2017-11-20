#ifndef rojue_PlotPreviewButton_h
#define rojue_PlotPreviewButton_h

#include "../coordinate_systems/rojue_CoordinateSystemOld.h"
#include "rojue_RButton.h"

namespace rojue
{

  /**

  This class is an PlotPreviewButton that can be used to switch between various editors and preview some
  plot from the editor.

  */

  class PlotPreviewButton : public RButton
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructs a preview button. */
    PlotPreviewButton(const juce::String& name, const CoordinateSystemOld* plotToPreview = NULL);    

    /** Destructor. */
    virtual ~PlotPreviewButton();

    //---------------------------------------------------------------------------------------------
    // appearance:

    /** Paints the button. */
    virtual void paint(Graphics &g);

    //=============================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    Image* plotPreviewImage;


  };

}

#endif  
