/*******************************************************************************
 The block below describes the properties of this module, and is read by
 the Projucer to automatically generate project code that uses it.
 For details about the syntax and how to create or use a module, see the
 JUCE Module Format.txt file.


 BEGIN_JUCE_MODULE_DECLARATION

  ID:               jura_framework
  vendor:           RS-MET
  version:          0.0.1
  name:             JUCE/RAPT audio software framework
  description:      GUI and I/O framework for developing RAPT based software
  website:          http://www.rs-met.com
  license:          GPL/Commercial

  dependencies:     juce_core, juce_audio_basics, juce_graphics, juce_gui_basics, 
                    juce_audio_formats, juce_audio_processors
  OSXFrameworks:    
  iOSFrameworks:    

 END_JUCE_MODULE_DECLARATION

*******************************************************************************/


#ifndef JURA_FRAMEWORK_H_INCLUDED
#define JURA_FRAMEWORK_H_INCLUDED

#include <juce_core/juce_core.h> 
#include <juce_audio_basics/juce_audio_basics.h> 
#include <juce_audio_formats/juce_audio_formats.h> 
#include <juce_audio_processors/juce_audio_processors.h> 
#include <juce_graphics/juce_graphics.h> 
#include <juce_gui_basics/juce_gui_basics.h> 
using namespace juce;

namespace jura
{

#include "control/jura_ScopedPointerLock.h"
#include "control/jura_Callbacks.h"
#include "control/jura_Parameter.h"
#include "control/jura_AutomatableParameter.h"
#include "control/jura_AutomatableModule.h"
#include "control/jura_StateManager.h"

#include "tools/jura_MiscTools.h"
#include "tools/jura_StringTools.h"
#include "tools/jura_XmlTools.h"

#include "files/jura_AudioFileInfo.h"
#include "files/jura_FileManager.h"
#include "files/jura_StateFileManager.h"
#include "files/jura_AudioFileManager.h"
#include "files/jura_ImageFileManager.h"

#include "gui/graphics/jura_ColourAHSL.h"
#include "gui/graphics/jura_ColourizableBitmap.h"
#include "gui/graphics/jura_ColourScheme.h"
#include "gui/graphics/jura_BitmapFonts.h"
#include "gui/graphics/jura_GraphicsTools.h"

#include "gui/misc/jura_DescribedComponent.h"
#include "gui/misc/jura_ColourSchemeComponent.h"
#include "gui/misc/jura_RectangleComponent.h"

#include "gui/widgets/jura_RWidget.h"
#include "gui/widgets/jura_WidgetSet.h"
#include "gui/widgets/jura_RTextField.h"
#include "gui/widgets/jura_RTextEditor.h" // maybe later it should be moved below RScrollBar...
#include "gui/widgets/jura_RButton.h"     // ..if the texteditor becomes scrollable
#include "gui/widgets/jura_RScrollBar.h"
#include "gui/widgets/jura_RTreeView.h"
#include "gui/widgets/jura_RPopUpComponent.h"
#include "gui/widgets/jura_RPopUpMenu.h"
#include "gui/widgets/jura_RSlider.h"
#include "gui/widgets/jura_RComboBox.h"
#include "gui/widgets/jura_StateLoadSaveWidgetSet.h"
// there are still some special widgets missing - copy them over soon...

#include "gui/plots/jura_CoordinateSystemOld.h"        // these should be renamed - get rid of the "Old"
#include "gui/plots/jura_CoordinateSystemZoomerOld.h"
#include "gui/plots/jura_CurveFamilyPlotOld.h" 
#include "gui/plots/jura_SpectrumDisplayOld.h"

// after the plots, we may add some further plot-based widgets, such as XY-Pads, frequency-response
// editors, etc.

#include "gui/editors/jura_Editor.h"

#include "gui/dialogs/jura_RDialogBox.h"
#include "gui/dialogs/jura_ColourSchemeSetupDialog.h"

}

#endif   // JURA_FRAMEWORK_H_INCLUDED
