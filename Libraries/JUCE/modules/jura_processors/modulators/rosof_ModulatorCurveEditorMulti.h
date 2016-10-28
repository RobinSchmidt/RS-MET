#ifndef rosof_ModulatorCurveEditorMulti_h
#define rosof_ModulatorCurveEditorMulti_h

#include "rosof_ModulatorCurveEditor.h"

namespace rosof
{

  /**

  This class ..

  */

  class ModulatorCurveEditorMulti	: public ModulatorCurveEditor
  {

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    ModulatorCurveEditorMulti(const juce::String& name = juce::String(T("ModulatorCurveEditorMulti")));   

    /** Destructor. */
    virtual ~ModulatorCurveEditorMulti(); 

    //-------------------------------------------------------------------------------------------------------------------------------------
    // setup:

    /** Adds a rosic::Breakpointmodulator object to the array of edited modulators. */
    virtual void addModulatorToEdit(rosic::BreakpointModulator* newModulatorToEdit);

    /** Selects one of the rosic::Breakpointmodulator objects for editing. */
    virtual void selectModulatorToEdit(int index);

    /** Sets a colour-scheme for one of the widget sets (associated to one of the modulators) - 
    after passing the pointer, the WidgetColourScheme will be owned by this object (inserted into an 
    OwnedArray), so make sure to not delete it yourself after passing it here. */
    //virtual void setWidgetColorScheme(int index, WidgetColourScheme* newColorScheme);

    /** Overrides CoordinateSystem::setCurrentRangeX to update all the plot-data. */
    virtual void setCurrentRangeX(double newMinX, double newMaxX);

    /** Overrides CoordinateSystem::setCurrentRange to update all the plot-data. */
    virtual void setCurrentRange(CoordinateSystemRangeOld newRange);

    /** Clears the array of curve-coulours. */
    virtual void clearCurveColours();

    /** Appends a colour for one of the curves to the end of the array. */
    virtual void appendCurveColour(const ColourAHSL& colourToAppend);

    /** Changes a colour for one of the curves - the index is assumed to exists in the array, if it's out of range, the function does 
    nothing. */
    virtual void changeCurveColour(int index, const ColourAHSL& newColour);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // inquiry:

    /** Overriden to make it possible to specify curve colours from outside. */
    virtual Colour getCurveColour(int index) const; 

    //-------------------------------------------------------------------------------------------------------------------------------------
    // callbacks:

    /** Overrides resized in order to update all curves instead of only the one of the modulator 
    that is currently being edited. */
    virtual void resized();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // others:

    /** Calls updatePlotCurveData(int, BreakpointModulator*, bool) for all modulators in the array
    and optionally redraws the curves and/or the coordinatesystem afterwards. */
    virtual void updateCurveDataForAllPlots(bool redrawCurves, bool redrawCoordinateSystem);

    /** Updates the maximum range of the inherited CoordinateSystem according to the postions of the
    most extreme breakpoints. */
    virtual void updateMaximumRange(bool alsoUpdateCurrentRange = false);

    //=====================================================================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    /** Overrides ModulatorCurveEditor::plotCurveFamily in order to draw all curves. */
    virtual void plotCurveFamily(Graphics &g, Image *targetImage = NULL, XmlElement *targetSVG = NULL);

    juce::Array<BreakpointModulator*,     CriticalSection> modulators;
    juce::OwnedArray<WidgetColourScheme*, CriticalSection> widgetColourSchemes;

    juce::Array<ColourAHSL> curveColours;

  };

}

#endif  
