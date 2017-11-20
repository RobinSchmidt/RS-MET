#ifndef rojue_ColourAHSL_h
#define rojue_ColourAHSL_h

#include "../includesForRojue.h"

namespace rojue
{

  /** This class represents a colour via the constituents opacity (alpha), hue, saturation and luminance (brightness). This representation
  is more intuitive to the user than the ARGB (alpha, red, green, blue) representation and thus lends itself better to designing 
  color-schemes for GUIs  */

  class ColourAHSL
  {

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. You can pass the initial values for hue, saturation, luminance and alpha here. If you pass nothing, it defaults
    to opaque black. Hue and saturation will default to 0.f such that when you later increase saturation and luminance, it will become some 
    shade of red (if you increase luminance only, if will turn into a shade of gray). */
    ColourAHSL(float hue = 0.f, float saturation = 0.f, float luminance = 0.f, float alpha = 1.f);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // setup:

    /** Sets the opacity (alpha) value. Values outside the range 0...1 are clipped. */
    void setAlpha(float newAlpha) { alpha = jlimit(0.f, 1.f, newAlpha); }

    /** Sets the hue value. Values outside the range 0...1 are periodically wrapped to make a color-circle. */
    void setHue(float newHue) { hue = fmod(newHue, 1.f); }

    /** Sets the saturation value. Values outside the range 0...1 are clipped. */
    void setSaturation(float newSaturation) { saturation = jlimit(0.f, 1.f, newSaturation); }

    /** Sets the luminance (brightness) value. Values outside the range 0...1 are clipped. */
    void setLuminance(float newLuminance) { luminance = jlimit(0.f, 1.f, newLuminance); }

    //-------------------------------------------------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the opacity (alpha) value. */
    float getAlpha() const { return alpha; }

    /** Returns the hue value. */
    float getHue() const { return hue; }

    /** Returns the saturation value. */
    float getSaturation() const { return saturation; }

    /** Returns the luminance (brightness) value. */
    float getLuminance() const { return luminance; }

    /** Converts this colour into the standard JUCE colour format (ARGB) and returns it. */
    Colour getAsJuceColour() const;

    //-------------------------------------------------------------------------------------------------------------------------------------
    // modifications:

    /** Returns a colour with a different hue with respect to this colour (but with all other parameters the same). */
    ColourAHSL withHue(float hueToUse) const 
    { 
      return ColourAHSL(hueToUse, saturation, luminance, alpha); 
    }

    /** Returns a colour with some hue-offset with respect to this colour. */
    ColourAHSL withHueOffset(float hueOffset) const 
    { 
      return ColourAHSL(hue+hueOffset, saturation, luminance, alpha); 
    }

    /** Returns a colour that has its saturation multiplied with some factor with respect to this colour. Note that when the factor is 
    larger than unity and this colour is already strongly saturated, it might be necesarry to clip the resulting saturation 1.0 */
    ColourAHSL withMultipliedSaturation(float multiplier) const 
    { 
      return ColourAHSL(hue, multiplier*saturation, luminance, alpha); 
    }

    /** Returns a colour that has its luminance raised to some power with respect to this colour. The power is commonly referred to as
    'gamma' in the image processing community. Above unity, it expands the contrast for bright colors and compresses it for dark colors, 
    below unity, the situation is vice versa. */
    ColourAHSL withLuminanceGamma(float gamma) const
    {
      return ColourAHSL(hue, saturation, pow(luminance, gamma), alpha); 
    }

    /** Returns a colour with several modifiers applied. For the meaning of each, see the corresponding functions for setting one at a 
    time. */
    ColourAHSL withModifiersApplied(float hueOffset, float saturationMulitplier, float luminanceGamma)
    {
      return ColourAHSL(hue+hueOffset, saturationMulitplier*saturation, pow(luminance, luminanceGamma), alpha);
    }

    //-------------------------------------------------------------------------------------------------------------------------------------
    // static conversion functions:

    /** Converts a color from HSL representation into RGB representation. */
    static void hslToRgb(float h, float s, float l, float &r, float &g, float &b);

    /** Converts a color from RGB representation into HSL representation. */
    static void rgbToHsl(float r, float g, float b, float &h, float &s, float &l);

    //=====================================================================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    float alpha, hue, saturation, luminance;

  };

}

#endif  