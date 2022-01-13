#ifndef jura_ColourScheme_h
#define jura_ColourScheme_h

/** This is a baseclass fo ColourSchemes for various GUI objects - like, for example, widgets, 
plots, editors, etc. The basclass defines some macro-parameters that can be applied to all these 
objects such as the general look (bright-on-dark vs. dark-on-bright), a hue offset, a saturation 
multiplier, etc. Whe the user changes one of these macro parameters, a call to updateColours() 
will be triggered - you need to override this function in your subclass to compute the actual 
colours.

default:
sky: dark-on-bright, hue= 0.65, saturation=0.60
coffee: b-o-d, h=0.09, s=0.40
mint: d-o-b, h=0.47, s=0.30
radar: b-o-d, h=0.35, s=0.40
rainstorm: bod, h=0.74, s=0.15
red wine: bod, h=0.02. s=0.6
wasp: bod, 0.14, 0.8
cool: bod, 0.55, 0.4  */

class JUCE_API ColourScheme
{

public:

  /** An enumeration of the general appearances. */
  enum appearances
  {
    DARK_ON_BRIGHT = 0,
    BRIGHT_ON_DARK
  };

  /** An enumeration of some predefined defaul color-schemes. */
  enum defaultSchemes
  {
    CHROME = 0,
    SKY,
    COFFEE,
    MINT,
    RADAR,
    RAINSTORM,
    REDWINE,
    WASP,
    COOL
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  ColourScheme()
  {
    appearance           = DARK_ON_BRIGHT;
    //appearance           = BRIGHT_ON_DARK;
    centralHue           = 0.f;
    saturationMultiplier = 0.f;
    brightnessGamma      = 1.f;
  }

  /** Destructor. */
  virtual ~ColourScheme() {}

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Selects one of the predefined appearances. @see appearances */
  virtual void setAppearance(int newAppearance) { appearance = newAppearance; updateColours(); }

  /** Selects one of the predefined appearances from a String. */
  virtual void setAppearanceFromString(const juce::String& newAppearanceString);

  /** Sets a central hue - all other hues will be defined with respect to this one. */
  virtual void setCentralHue(float newHue) { centralHue = newHue; updateColours(); }

  /** Sets a global multiplier for all saturations. */
  virtual void setSaturationMultiplier(float newMultiplier)
  {
    saturationMultiplier = newMultiplier;
    updateColours();
  }

  /** Sets a gamma value for all luminances. */
  virtual void setBrightnessGamma(float newGamma) { brightnessGamma = newGamma; updateColours(); }

  /** Clears the array of hue-offsets. */
  virtual void clearHueOffsets() { hueOffsets.clear(); }

  /** Appends a hue-offset into our array of hue-offsets (which are used for various parts of the
  editor). */
  virtual void appendHueOffset(float offsetToAppend) { hueOffsets.add(offsetToAppend); }

  /** Sets up one of the hue-offsets with the given index. If the given index does not yet exist in
  our array, the new offset will be appended, possibly filling up with zeros at indices in between
  the current max-index and the new index. */
  virtual void setHueOffset(int index, float newOffset)
  {
    if(index >= hueOffsets.size())
    {
      while(index > hueOffsets.size())
        appendHueOffset(0.f);
      appendHueOffset(newOffset);
    }
    else
      hueOffsets.set(index, newOffset);
  }

  // setBrightness, setContrast, setHueSpread

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the index of the currently selected appearance. @see appearances */
  virtual int getAppearance() const { return appearance; }

  /** Returns a String representing the currently selected appearance. */
  virtual juce::String getAppearanceString() const;

  /** Returns the central hue. @see setCentralHue */
  virtual float getCentralHue() const { return centralHue; }

  /** Returns the multiplier for all saturations. @see setSaturationMultiplier */
  virtual float getSaturationMultiplier() const { return saturationMultiplier; }

  /** Returns the gamma value for all luminances. @see setBrightnessGamma */
  virtual float getBrightnessGamma() const { return brightnessGamma; }

  /** Returns the number of defined hue-offsets. These offsets may be used for various parts of
  the editor or sub-editors. */
  virtual int getNumHueOffsets() const { return hueOffsets.size(); }

  /** Returns the hue-offset at the given index. These offsets may be used for various parts of
  the editor or sub-editors. */
  virtual float getHueOffset(int index) const
  {
    if(index >= hueOffsets.size())
      return 0.f;
    else
      return hueOffsets[index];
  }

  juce_UseDebuggingNewOperator;

protected:

  /** Method that will be called when the user has changed one of our macro-parameters - subclasses
  should update their colours then accordingly. */
  virtual void updateColours() = 0;

  int   appearance;
  float centralHue, saturationMultiplier, brightnessGamma; // contrast, brightness, gamma, hueSpread
  juce::Array<float> hueOffsets;                           // offsets for hues for various parts of the editor

};

//=================================================================================================

/** A ColourScheme for editors defining colours for the background-gradient, headline, etc.  */

class JUCE_API EditorColourScheme : public ColourScheme
{

public:

  EditorColourScheme() 
  { 
    //appearance = BRIGHT_ON_DARK;
    appearance = DARK_ON_BRIGHT;
    updateColours(); 
  }

  Colour topLeft, topRight, bottomLeft, bottomRight, outline, headline, headlineOutline, text;

  juce_UseDebuggingNewOperator;

protected:

  virtual void updateColours();

};

//=================================================================================================

/** A ColourScheme for widgets defining colours for the background, outline, etc.  */

class JUCE_API WidgetColourScheme : public ColourScheme
{

public:

  WidgetColourScheme() 
  { 
    appearance = BRIGHT_ON_DARK;
    //appearance = DARK_ON_BRIGHT;
    updateColours(); 
  }

  Colour background, outline, handle, text, weakHighlight, strongHighlight,
    special; 
  // \todo define mouseOverOutline colour, get rid of special, introduce weakHighlight, 
  // strongHighlight

  juce_UseDebuggingNewOperator;

protected:

  virtual void updateColours();

};

//=================================================================================================

/** A ColourScheme for plots defining colours for the background, curves, etc.  */

class JUCE_API PlotColourScheme : public ColourScheme
{

public:

  enum curveColouringStrategies
  {
    UNIFORM,        ///< all curves have the same colour
    ALTERNATING     ///< alternating around the central hue with increasing deviations
  };

  PlotColourScheme()
  {
    appearance = BRIGHT_ON_DARK;
    curveColouringStrategy = ALTERNATING;
    curveHueSpread         = 0.125;
    updateColours();
  }

  /** Selects one of the strategies to colour multiple curves in the plot. This will determine what
  Colour will be returned by subsequent call to getCurveColour @see curveColouringStrategies */
  void setCurveColouringStrategy(int newStrategy)
  {
    curveColouringStrategy = newStrategy;
    updateColours();
  }

  /** Returns the colour for the curve with the given index. */
  Colour getCurveColour(int index) const;


  Colour getCurveColourUniform(int index) const;
  Colour getCurveColourAlternating(int index) const;


  Colour     topLeft, topRight, bottomLeft, bottomRight, outline, text, axes, coarseGrid, fineGrid; //, curves;
  ColourAHSL curvesAHSL;

  juce_UseDebuggingNewOperator;

protected:

  virtual void updateColours();

  int   curveColouringStrategy;
  float curveHueSpread;

};

#endif
