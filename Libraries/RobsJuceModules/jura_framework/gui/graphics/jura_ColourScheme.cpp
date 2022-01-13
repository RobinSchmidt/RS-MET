//#include "rojue_ColourScheme.h"
//using namespace rojue;

void ColourScheme::setAppearanceFromString(const juce::String& newAppearanceString)
{
  if( newAppearanceString == String("DarkOnBright") )
    setAppearance(DARK_ON_BRIGHT);
  else if( newAppearanceString == String("BrightOnDark") )
    setAppearance(BRIGHT_ON_DARK);
  else
    jassertfalse;  // invalid string
}

juce::String ColourScheme::getAppearanceString() const
{
  switch( appearance )
  {
  case DARK_ON_BRIGHT: return String("DarkOnBright");
  case BRIGHT_ON_DARK: return String("BrightOnDark");
  default: return String();
  }
}

void EditorColourScheme::updateColours()
{
  ColourAHSL backgroundAHSL, topLeftAHSL, topRightAHSL, bottomLeftAHSL, bottomRightAHSL,
      outlineAHSL, headlineAHSL, headlineOutlineAHSL, textAHSL;

  switch( appearance )
  {
  case DARK_ON_BRIGHT:
    {
      topLeftAHSL.setLuminance(        1.0);
      topRightAHSL.setLuminance(       0.8125);
      bottomLeftAHSL.setLuminance(     0.8125);
      //topRightAHSL.setLuminance(       0.75);
      //bottomLeftAHSL.setLuminance(     0.75);
      bottomRightAHSL.setLuminance(    1.0);
      outlineAHSL.setLuminance(        0.5);
      headlineAHSL.setLuminance(       0.3125);
      headlineOutlineAHSL.setLuminance(1.0);
      textAHSL.setLuminance(           0.25);
    }
    break;
  case BRIGHT_ON_DARK:
    {
      topLeftAHSL.setLuminance(        0.375);
      //topRightAHSL.setLuminance(       0.25);
      //bottomLeftAHSL.setLuminance(     0.25);
      topRightAHSL.setLuminance(       0.1875);
      bottomLeftAHSL.setLuminance(     0.1875);
      bottomRightAHSL.setLuminance(    0.375);
      outlineAHSL.setLuminance(        0.5);
      headlineAHSL.setLuminance(       0.8125);
      headlineOutlineAHSL.setLuminance(0.0);
      textAHSL.setLuminance(           0.75);
    }
    break;
  }

  // saturations (preliminary):
  topLeftAHSL.setSaturation(        1.0);
  topRightAHSL.setSaturation(       1.0);
  bottomLeftAHSL.setSaturation(     1.0);
  bottomRightAHSL.setSaturation(    1.0);
  outlineAHSL.setSaturation(        1.0);
  headlineAHSL.setSaturation(       1.0);
  headlineOutlineAHSL.setSaturation(1.0);
  textAHSL.setSaturation(1.0);

  // apply global modifiers (macro parameters):
  float ch = centralHue;
  float sm = saturationMultiplier;
  float bg = brightnessGamma;
  topLeftAHSL         = topLeftAHSL.withModifiersApplied(        ch, sm, bg);
  topRightAHSL        = topRightAHSL.withModifiersApplied(       ch, sm, bg);
  bottomLeftAHSL      = bottomLeftAHSL.withModifiersApplied(     ch, sm, bg);
  bottomRightAHSL     = bottomRightAHSL.withModifiersApplied(    ch, sm, bg);
  outlineAHSL         = outlineAHSL.withModifiersApplied(        ch, sm, bg);
  headlineAHSL        = headlineAHSL.withModifiersApplied(       ch, sm, bg);
  headlineOutlineAHSL = headlineOutlineAHSL.withModifiersApplied(ch, sm, bg);
  textAHSL            = textAHSL.withModifiersApplied(ch, sm, bg);

  // convert to ARGB and store in members::
  topLeft         = topLeftAHSL.getAsJuceColour();
  topRight        = topRightAHSL.getAsJuceColour();
  bottomLeft      = bottomLeftAHSL.getAsJuceColour();
  bottomRight     = bottomRightAHSL.getAsJuceColour();
  outline         = outlineAHSL.getAsJuceColour();
  headline        = headlineAHSL.getAsJuceColour();
  headlineOutline = headlineOutlineAHSL.getAsJuceColour();
  text            = textAHSL.getAsJuceColour();
}

void WidgetColourScheme::updateColours()
{
  ColourAHSL backgroundAHSL, outlineAHSL, handleAHSL, textAHSL, specialAHSL;

  switch( appearance )
  {
  case DARK_ON_BRIGHT:
    {
      backgroundAHSL.setLuminance(0.9375);
      outlineAHSL.setLuminance(   0.5);
      //handleAHSL.setLuminance(    0.625);
      handleAHSL.setLuminance(    0.6875);
      //handleAHSL.setLuminance(    0.75);
      textAHSL.setLuminance(      0.f);
      textAHSL.setAlpha(          0.8125f);
      specialAHSL.setLuminance(   1.f); // preliminary
    }
    break;
  case BRIGHT_ON_DARK:
    {
      backgroundAHSL.setLuminance(0.0);
      outlineAHSL.setLuminance(   0.5);
      handleAHSL.setLuminance(    0.435f);
      textAHSL.setLuminance(      1.f);
      textAHSL.setAlpha(          0.8125f);
      specialAHSL.setLuminance(   1.f); // preliminary

      //outlineAHSL.setHue(-0.0625);  // experimental
      //handleAHSL.setHue( +0.0625);
    }
    break;
  }

  // saturations (preliminary):
  backgroundAHSL.setSaturation(1.0);
  outlineAHSL.setSaturation(   1.0);
  handleAHSL.setSaturation(    1.0);
  textAHSL.setSaturation(      1.0);
  specialAHSL.setSaturation(   1.0);

  float ch = centralHue;
  float sm = saturationMultiplier;
  float bg = brightnessGamma;
  backgroundAHSL = backgroundAHSL.withModifiersApplied(ch, sm, bg);
  outlineAHSL    = outlineAHSL.withModifiersApplied(   ch, sm, bg);
  handleAHSL     = handleAHSL.withModifiersApplied(    ch, sm, bg);
  textAHSL       = textAHSL.withModifiersApplied(      ch, sm, bg);
  specialAHSL    = specialAHSL.withModifiersApplied(   ch, sm, bg);

  background = backgroundAHSL.getAsJuceColour();
  outline    = outlineAHSL.getAsJuceColour();
  handle     = handleAHSL.getAsJuceColour();
  text       = textAHSL.getAsJuceColour();
  special    = specialAHSL.getAsJuceColour();
}

void PlotColourScheme::updateColours()
{
  ColourAHSL topLeftAHSL, topRightAHSL, bottomLeftAHSL, bottomRightAHSL, outlineAHSL, textAHSL,
    axesAHSL, coarseGridAHSL, fineGridAHSL; //, curvesAHSL;

  switch( appearance )
  {
  case DARK_ON_BRIGHT:
    {
      topLeftAHSL.setLuminance(    1.0);
      bottomRightAHSL.setLuminance(1.0);
      topRightAHSL.setLuminance(   1.0);
      bottomLeftAHSL.setLuminance( 1.0);
      //topRightAHSL.setLuminance(   0.9375);  // this is for the background with the gradient - too slow
      //bottomLeftAHSL.setLuminance( 0.9375);

      textAHSL.setLuminance(       0.0);
      outlineAHSL.setLuminance(    0.5);
      axesAHSL.setLuminance(       0.1875);
      coarseGridAHSL.setLuminance( 0.75);
      fineGridAHSL.setLuminance(   0.875);
      //curvesAHSL.setLuminance(     0.f);
      curvesAHSL.setLuminance(     0.25f);
      curvesAHSL.setAlpha(         0.8125f);
    }
    break;
  case BRIGHT_ON_DARK:
    {
      //topLeftAHSL.setLuminance(    0.125);  // this is for the background with the gradient - too slow
      //bottomRightAHSL.setLuminance(0.125);
      topLeftAHSL.setLuminance(    0.0);
      bottomRightAHSL.setLuminance(0.0);
      topRightAHSL.setLuminance(   0.0);
      bottomLeftAHSL.setLuminance( 0.0);

      //topLeftAHSL.setLuminance(    0.1875);
      //topRightAHSL.setLuminance(   0.0);
      //bottomLeftAHSL.setLuminance( 0.0);
      //bottomRightAHSL.setLuminance(0.1875);
      outlineAHSL.setLuminance(    0.5);
      textAHSL.setLuminance(       0.8125);
      axesAHSL.setLuminance(       0.8125);
      coarseGridAHSL.setLuminance( 0.3125);
      fineGridAHSL.setLuminance(   0.1875);
      //curvesAHSL.setLuminance(     1.f);
      curvesAHSL.setLuminance(     0.875f);
      curvesAHSL.setAlpha(         0.8125f);
    }
    break;
  }

  // saturations (preliminary):
  topLeftAHSL.setSaturation(    1.0);
  topRightAHSL.setSaturation(   1.0);
  bottomLeftAHSL.setSaturation( 1.0);
  bottomRightAHSL.setSaturation(1.0);
  outlineAHSL.setSaturation(    1.0);
  textAHSL.setSaturation(       1.0);
  axesAHSL.setSaturation(       1.0);
  coarseGridAHSL.setSaturation( 1.0);
  fineGridAHSL.setSaturation(   1.0);
  curvesAHSL.setSaturation(     1.0);

  curvesAHSL.setHue(0.f); // needs to be done because this is a member

  float ch = centralHue;
  float sm = saturationMultiplier;
  float bg = brightnessGamma;
  topLeftAHSL     = topLeftAHSL.withModifiersApplied(    ch, sm, bg);
  topRightAHSL    = topRightAHSL.withModifiersApplied(   ch, sm, bg);
  bottomLeftAHSL  = bottomLeftAHSL.withModifiersApplied( ch, sm, bg);
  bottomRightAHSL = bottomRightAHSL.withModifiersApplied(ch, sm, bg);
  outlineAHSL     = outlineAHSL.withModifiersApplied(    ch, sm, bg);
  textAHSL        = textAHSL.withModifiersApplied(       ch, sm, bg);
  axesAHSL        = axesAHSL.withModifiersApplied(       ch, sm, bg);
  coarseGridAHSL  = coarseGridAHSL.withModifiersApplied( ch, sm, bg);
  fineGridAHSL    = fineGridAHSL.withModifiersApplied(   ch, sm, bg);
  curvesAHSL      = curvesAHSL.withModifiersApplied(     ch, sm, bg);

  topLeft     = topLeftAHSL.getAsJuceColour();
  topRight    = topRightAHSL.getAsJuceColour();
  bottomLeft  = bottomLeftAHSL.getAsJuceColour();
  bottomRight = bottomRightAHSL.getAsJuceColour();
  outline     = outlineAHSL.getAsJuceColour();
  text        = textAHSL.getAsJuceColour();
  axes        = axesAHSL.getAsJuceColour();
  coarseGrid  = coarseGridAHSL.getAsJuceColour();
  fineGrid    = fineGridAHSL.getAsJuceColour();
  //curves      = curvesAHSL.getAsJuceColour();
}

Colour PlotColourScheme::getCurveColour(int index) const
{
  switch( curveColouringStrategy )
  {
  case UNIFORM:     return getCurveColourUniform(index);
  case ALTERNATING: return getCurveColourAlternating(index);

  default:
    {
      jassertfalse;
      return Colours::red;
    }
  }
}

Colour PlotColourScheme::getCurveColourUniform(int index) const
{
  return curvesAHSL.getAsJuceColour();
}

Colour PlotColourScheme::getCurveColourAlternating(int index) const
{
  float sign   = pow(-1.f, index);
  float factor = (float) (index/2 + 1);
  //float hue = curvesAHSL.getHue() + sign*factor*curveHueSpread;
  //ColourAHSL tmpColourAHSL = curvesAHSL.withHueOffset(sign*factor*curveHueSpread);
  return curvesAHSL.withHueOffset(sign*factor*curveHueSpread).getAsJuceColour();
}


/*

ToDo:

-maybe introduce a slider to scale the amount of the gradient for the background from
 -100%...+100%
-maybe have a similar slider for the plots
-maybe the sliders/widgets can also have so sort of gradient (but probably not the bilinear one,
 maybe more some sort of pseudo 3D gradient - maybe vertically from dark to bright to dark, the 
 buttons maybe with darker outer regions...or maybe generally, the fill-ins should have darker
 outer regions


*/