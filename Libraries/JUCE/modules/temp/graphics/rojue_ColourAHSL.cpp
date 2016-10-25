#include "rojue_ColourAHSL.h"
using namespace rojue;

//-----------------------------------------------------------------------------------------------------------------------------------------    
// construction/destruction:

ColourAHSL::ColourAHSL(float hue, float saturation, float luminance, float alpha)
{
  this->alpha      = jlimit(0.f, 1.f, alpha);
  this->saturation = jlimit(0.f, 1.f, saturation);
  this->luminance  = jlimit(0.f, 1.f, luminance);

  while( hue < 0.f )
    hue += 1.f;
  while( hue > 1.f )
    hue -= 1.f;
  this->hue = hue;
}

//-----------------------------------------------------------------------------------------------------------------------------------------    
// conversions:

void ColourAHSL::hslToRgb(float h, float s, float l, float &r, float &g, float &b)
{
  float rSat, gSat, bSat, rTmp, gTmp, bTmp;

  //...i think we need to convert our hue (0...1) to 0...360 first...
  h *= 360.f;

  while (h < 0)
    h += 360;
  while (h > 360)
    h -= 360;


  if(h < 120) 
  {
    rSat = (120 - h) / 60.f;
    gSat = h / 60.f;
    bSat = 0;
  } 
  else if (h < 240) 
  {
    rSat = 0;
    gSat = (240 - h) / 60.f;
    bSat = (h - 120) / 60.f;
  } 
  else 
  {
    rSat = (h - 240) / 60.f;
    gSat = 0;
    bSat = (360 - h) / 60.f;
  }

  rSat = jmin(rSat, 1.f);
  gSat = jmin(gSat, 1.f);
  bSat = jmin(bSat, 1.f);

  rTmp = 2 * s * rSat + (1 - s);
  gTmp = 2 * s * gSat + (1 - s);
  bTmp = 2 * s * bSat + (1 - s);

  if (l < 0.5) 
  {
    r = l * rTmp;
    g = l * gTmp;
    b = l * bTmp;
  } 
  else 
  {
    r = (1 - l) * rTmp + 2 * l - 1;
    g = (1 - l) * gTmp + 2 * l - 1;
    b = (1 - l) * bTmp + 2 * l - 1;
  }
}

void ColourAHSL::rgbToHsl(float r, float g, float b, float &h, float &s, float &l)
{
   float min   = jmin(r, g, b);
   float max   = jmax(r, g, b);
   float delta = max - min;

   l = 0.5f * (min + max);

   s = 0.f;
   if(l > 0.f && l < 1.f)
      s = delta / (l < 0.5 ? (2*l) : (2-2*l));

   h = 0.f;
   if(delta > 0.f) 
   {
      if(max == r && max != g)   
        h += (g - b) / delta;
      if(max == g && max != b)      
        h += (2 + (b - r) / delta);
      if(max == b && max != r)
        h += (4 + (r - g) / delta);
      h *= 60;

      //...i think we need to convert this hue (0...360) to 0...1 now...
      h /= 360.f;

   }
}

//-----------------------------------------------------------------------------------------------------------------------------------------    
// inquiry:

Colour ColourAHSL::getAsJuceColour() const
{
  float rf, gf, bf;
  hslToRgb(hue, saturation, luminance, rf, gf, bf);
  uint8 r = (uint8) (255.f * rf);
  uint8 g = (uint8) (255.f * gf);
  uint8 b = (uint8) (255.f * bf);
  return Colour(r, g, b, alpha);
}


