#ifndef rojue_Panel_h
#define rojue_Panel_h

#include "rojue_PanelRange.h"

namespace rojue
{

  /**

  This class is a component, intended to be used as base-class for all components that need some
  underlying coordinate-system (albeit possibly implicit and invisble). It provides means to set
  up maximum values in the x- and y-direction, to select a currently shown range within this 
  maximum range and has functions for conversion between coordinates in the panel's own intrinsic 
  coordinate system and component-coordinates. 

  */

  class Panel	: virtual public Component
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    Panel(const juce::String &componentName = juce::String(T("Panel")));   

    /** Destructor. */
    virtual ~Panel();            

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets the maximum for the currently visible range. For logarithmic x- and/or y-axis-scaling, 
    make sure that the respective minimum value is greater than zero! */
    virtual void setMaximumRange(double newMinX, double newMaxX, double newMinY, double newMaxY);

    /** Sets the maximum for the currently visible range. */
    virtual void setMaximumRange(PanelRange newMaximumRange);

    /** Sets the maximum visible range for the x-axis. */
    virtual void setMaximumRangeX(double newMinX, double newMaxX);

    /** Sets the maximum visible range for the y-axis. */
    virtual void setMaximumRangeY(double newMinY, double newMaxY);

    /** Sets the minimum value for the range of x. */
    virtual void setMaximumRangeMinX(double newMinX);

    /** Sets the maximum value for the range of x. */
    virtual void setMaximumRangeMaxX(double newMaxX);

    /** Sets the minimum value for the range of y. */
    virtual void setMaximumRangeMinY(double newMinY);

    /** Sets the maximum value for the range of y. */
    virtual void setMaximumRangeMaxY(double newMaxY);

    /** Sets the currently visible range. For logarithmic x- and/or y-axis-scaling, make sure that 
    the respective minimum value is greater than zero! */
    virtual void setCurrentRange(double newMinX, double newMaxX, double newMinY, double newMaxY);

    /** Sets the currently visible range. */
    virtual void setCurrentRange(PanelRange newRange);

    /** Sets the currently visible range for the y-axis. */
    virtual void setCurrentRangeX(double newMinX, double newMaxX);

    /** Sets the currently visible range for the y-axis. */
    virtual void setCurrentRangeY(double newMinY, double newMaxY);

    /** Sets the minimum value of x. */
    virtual void setCurrentRangeMinX(double newMinX);

    /** Sets the maximum value of x. */
    virtual void setCurrentRangeMaxX(double newMaxX);

    /** Sets the minimum value of y. */
    virtual void setCurrentRangeMinY(double newMinY);

    /** Sets the maximum value of y. */
    virtual void setCurrentRangeMaxY(double newMaxY);

    /** Sets a minimum value for the width of the panel measured in the intrinsic coordinate 
    system. */
    virtual void setMinimumWidth(double newMinWidth);

    /** Sets a minimum value for the height of the panel measured in the intrinsic coordinate 
    system. */
    virtual void setMinimumHeight(double newMinHeight);

    /** Sets the interval of the horizontal coarse grid. */
    virtual void setHorizontalCoarseGridInterval(double newGridInterval); 

    /** Sets the interval of the horizontal fine grid. */
    virtual void setHorizontalFineGridInterval(double newGridInterval); 

    /** Sets the interval of the vertical coarse grid. */
    virtual void setVerticalCoarseGridInterval(double newGridInterval); 

    /** Sets the interval of the vertical fine grid. */
    virtual void setVerticalFineGridInterval(double newGridInterval); 

    /** With this function, the automatic re-rendering of the underlying image can be turned on or 
    off. If on (default), everytime you change a parameter which will change the appearance of 
    the Panel, it will be re-rendered. However, when you want to change many parameters at a time, 
    this can be wasteful in terms of CPU-load. In these cases, it can be useful to switch the 
    automatic re-rendering temporarily off. */
    //virtual void setAutoReRendering(bool shouldAutomaticallyReRender);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the maximum for the currently visible range. */
    virtual PanelRange getMaximumRange() const { return maximumRange; }

    /** Returns the minimum value for the range of x. */
    virtual double getMaximumRangeMinX() const { return maximumRange.getMinX(); }

    /** Returns the maximum value for the range of x. */
    virtual double getMaximumRangeMaxX() const { return maximumRange.getMaxX(); }

    /** Returns the minimum value for the range of y. */
    virtual double getMaximumRangeMinY() const { return maximumRange.getMinY(); }

    /** Returns the maximum value for the range of y. */
    virtual double getMaximumRangeMaxY() const { return maximumRange.getMaxY(); }

    /** Returns the currently visible range. */
    virtual PanelRange getCurrentRange() const { return currentRange; }

    /** Returns the minimum value of x. */
    virtual double getCurrentRangeMinX() const { return currentRange.getMinX(); }

    /** Returns the maximum value of x. */
    virtual double getCurrentRangeMaxX() const { return currentRange.getMaxX(); }

    /** Returns the minimum value of y. */
    virtual double getCurrentRangeMinY() const { return currentRange.getMinY(); }

    /** Returns the maximum value of y. */
    virtual double getCurrentRangeMaxY() const { return currentRange.getMaxY(); }

    /** Returns the minimum width to be displayed. */
    virtual double getMinimumWidth() const { return minimumWidth; }

    /** Returns the minimum height to be displayed. */
    virtual double getMinimumHeight() const { return minimumHeight; }

    /** Returns the interval of the horizontal coarse grid. */
    virtual double getHorizontalCoarseGridInterval() const { return horizontalCoarseGridInterval; }

    /** Returns the interval of the horizontal fine grid. */
    virtual double getHorizontalFineGridInterval() const { return horizontalFineGridInterval; }

    /** Returns the interval of the vertical coarse grid. */
    virtual double getVerticalCoarseGridInterval() const { return verticalCoarseGridInterval; }

    /** Returns the interval of the vertical fine grid. */
    virtual double getVerticalFineGridInterval() const { return verticalFineGridInterval; }

    //---------------------------------------------------------------------------------------------
    // coordinate-transformations:

    /** Function for converting the x- and y-coordinate values into the corresponding coordinates 
    in the component (double precision version).*/
    virtual void transformToComponentsCoordinates(double &x, double &y) const;

    /** Function for converting the x- and y-coordinate values into the corresponding coordinates 
    in the component (single precision version).*/
    virtual void transformToComponentsCoordinates(float &x, float &y) const;

    /** Function for converting the x- and y-coordinate values measured in the components 
    coordinate system to the corresponding coordinates of our plot (double precision version). */
    virtual void transformFromComponentsCoordinates(double &x, double &y) const;

    /** Function for converting the x- and y-coordinate values measured in the components 
    coordinate system to the corresponding coordinates of our plot (single precision version). */
    virtual void transformFromComponentsCoordinates(float &x, float &y) const;

    //---------------------------------------------------------------------------------------------
    // callbacks:

    /** Overrides the resized()-function of the component base-class. */
    virtual void resized();

    //---------------------------------------------------------------------------------------------
    // others:

    /** This function should be called by member functions that change some parameter such that a 
    re-drawing is required. It will only set the dirty flag but not itself initiate any 
    re-drawing. */
    virtual void setDirty(bool shouldSetToDirty = true);

    //=============================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    /** Overrides drawComponent inherited from ThreadedDrawingComponent in order to do the actual 
    drawing operations. The implementation here will only draw a placeholder. */
    //virtual void drawComponent(Image* imageToDrawOnto);

    /** Constrains the currentRange to not not exceed to maximumRange and to ensure a minimumWidth 
    and minimumHeight. By default, this will also call setDirty() but that may be supressed by 
    passing false to the parameter. */
    virtual void constrainCurrentRange(bool callSetDirty = true);

    /** Updates the scale-factors which are needed when transforming from the Panel's 
    coordinates to Component's coordinates and vice versa.  */
    virtual void updateScaleFactors();

    /** The current- and maximum-range object for the panel (in intrinsic coordinates). */
    PanelRange currentRange, maximumRange; 

    /** Some variables to ensure a minimum width and height of the intrinsic coordinate system. */
    double minimumWidth, minimumHeight;

    double scaleX;
    double scaleY;
    double horizontalCoarseGridInterval;
    double horizontalFineGridInterval;
    double verticalCoarseGridInterval;
    double verticalFineGridInterval;

    static const double minGridInterval;

  };

}

#endif  
