#ifdef JURA_FRAMEWORK_H_INCLUDED
 /* When you add this cpp file to your project, you mustn't include it in a file where you've
    already included any other headers - just put it inside a file on its own, possibly with your config
    flags preceding it, but don't include anything else. That also includes avoiding any automatic prefix
    header files that the compiler may be using.
 */
 #error "Incorrect use of JUCE cpp file"
#endif

#include "jura_framework.h"

namespace jura
{
#include "misc/jura_DeletionRequester.cpp"

#include "tools/jura_MiscTools.cpp"
#include "tools/jura_StringTools.cpp"
#include "tools/jura_XmlTools.cpp"

#include "control/jura_Mediator.cpp"
#include "control/jura_Parameter.cpp"
#include "control/jura_ParameterProfiles.cpp"
#include "control/jura_SmoothableParameter.cpp"
#include "control/jura_AutomatableParameter.cpp"
#include "control/jura_MetaParameter.cpp"
#include "control/jura_VoiceManager.cpp"
#include "control/jura_ModulatableParameter.cpp"
#include "control/jura_PreDefinedParameters.cpp"
#include "control/jura_ParameterManager.cpp"
#include "control/jura_StateManager.cpp"

#include "files/jura_AudioFileInfo.cpp"
#include "files/jura_FileManager.cpp"
#include "files/jura_StateFileManager.cpp"
#include "files/jura_AudioFileManager.cpp"
#include "files/jura_ImageFileManager.cpp"

#include "gui/graphics/jura_ColorMap.cpp"
#include "gui/graphics/jura_ColourAHSL.cpp"
#include "gui/graphics/jura_ColourizableBitmap.cpp"
#include "gui/graphics/jura_ColourScheme.cpp"
#include "gui/graphics/jura_BitmapFonts.cpp"
#include "gui/graphics/jura_GraphicsTools.cpp"
#include "gui/graphics/jura_ImageUpdater.cpp"

#include "gui/misc/jura_DescribedComponent.cpp"
#include "gui/misc/jura_ColourSchemeComponent.cpp"
#include "gui/misc/jura_RectangleComponent.cpp"
#include "gui/misc/jura_MessageBoxes.cpp"
#include "gui/misc/jura_RepaintManager.cpp"

#include "gui/widgets/jura_RWidget.cpp"
#include "gui/widgets/jura_RTextField.cpp"
#include "gui/widgets/jura_RTextEditor.cpp"
#include "gui/widgets/jura_RButton.cpp"
#include "gui/widgets/jura_PlotPreviewButton.cpp" 
#include "gui/widgets/jura_RScrollBar.cpp"
#include "gui/widgets/jura_RTreeView.cpp"
#include "gui/widgets/jura_RPopUpComponent.cpp"
#include "gui/widgets/jura_RPopUpMenu.cpp"
#include "gui/widgets/jura_RSlider.cpp"
#include "gui/widgets/jura_RDraggableNumber.cpp"
#include "gui/widgets/jura_RComboBox.cpp"
#include "gui/widgets/jura_RTimeGridComboBox.cpp"
#include "gui/widgets/jura_RSyncIntervalComboBox.cpp"
#include "gui/widgets/jura_NodeEditor.cpp"
#include "gui/widgets/jura_AutomatableWidget.cpp"
#include "gui/widgets/jura_VectorPad.cpp"

#include "gui/misc/jura_ComponentScrollContainer.cpp"

#include "gui/widgets/widget_sets/jura_WidgetSet.cpp"
#include "gui/widgets/widget_sets/jura_StateLoadSaveWidgetSet.cpp"
#include "gui/widgets/widget_sets/jura_FileSelectionBox.cpp" 
#include "gui/widgets/widget_sets/jura_ColorMapLoader.cpp"

#include "gui/plots/jura_PlotRange.cpp"
#include "gui/plots/jura_PlotSettings.cpp"
#include "gui/plots/jura_PlotDrawer.cpp"
#include "gui/plots/jura_Plot.cpp"
#include "gui/plots/jura_PlotEditor.cpp" 
#include "gui/plots/jura_ObservablePlot.cpp"
#include "gui/plots/jura_PlotZoomer.cpp"
#include "gui/plots/jura_DataPlot.cpp"
#include "gui/plots/jura_FunctionPlot.cpp" 
#include "gui/plots/jura_SpectrumPlot.cpp"
#include "gui/plots/jura_WaveformPlot.cpp"

// AudioBufferUser stuff (needed by waveform display)
#include "audio/jura_AudioFileBuffer.cpp"  
#include "audio/jura_AudioFileBufferUser.cpp" 

#include "gui/panels/jura_PanelRange.cpp"
#include "gui/panels/jura_Panel.cpp"
#include "gui/panels/jura_DrawingThread.cpp"
#include "gui/panels/jura_ThreadedDrawingComponent.cpp"
#include "gui/panels/jura_ThreadedDrawingPanel.cpp"
#include "gui/panels/jura_CoordinateSystem.cpp"
#include "gui/panels/jura_CurveFamilyPlot.cpp"
#include "gui/panels/jura_CoordinateSystemZoomer.cpp"
#include "gui/panels/jura_InteractiveCoordinateSystem.cpp"
#include "gui/panels/jura_WaveformDisplay.cpp"
#include "gui/panels/jura_DualWaveformDisplay.cpp"

#include "gui/dialogs/jura_RDialogBox.cpp"
#include "gui/dialogs/jura_RMessageBox.cpp"
#include "gui/dialogs/jura_ColourSchemeSetupDialog.cpp"
#include "gui/dialogs/jura_ImageSavingDialog.cpp"

#include "gui/editors/jura_Editor.cpp"

#include "audio/jura_AudioSampleBufferFunctions.cpp"
#include "audio/jura_ImmediatePlaybackAudioSource.cpp"
#include "audio/jura_AudioModule.cpp" // needs editor stuff
#include "audio/jura_AudioModuleFactory.cpp"
#include "audio/jura_AudioModuleSelector.cpp"
//#include "audio/jura_PolyModule.cpp"
//#include "audio/jura_PolySlot.cpp"
//#include "audio/jura_PolyModulators.cpp"
#include "audio/jura_AudioPlugin.cpp"

#include "misc/jura_Experimental.cpp"

}
