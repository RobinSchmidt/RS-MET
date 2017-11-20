#ifndef RS_COLORSPACES_H
#define RS_COLORSPACES_H

namespace RSLib
{

  /**

  This class represents a color with components red, green, blue and an opacity (aka alpha) value
  using an 8-bit unsigned integer (i.e. a byte) for each component. The order in which the
  components appear is memory may be platform dependent in order to match the ordering that the
  respective native blitting function needs.

  */

  class RSLib_API rsColorRGBA
  {

  public:

    enum colorIndices
    {
      black,
      blue,
      cyan,
      green,
      magenta,
      red,
      transparent,
      white,
      yellow
      //....
    };


    /** \name Construction/Destruction */

    /** Default constructor. You may initialize the components as you like or accept the default
    values which are all zero, defining a transparent black as default color. */
    rsColorRGBA(rsUint8 red = 0, rsUint8 green = 0, rsUint8 blue = 0, rsUint8 alpha = 0);

    /** Creates a color from an index. @see colorIndices. */
    rsColorRGBA(int index);

    /** Creates a color from a floating point number that should represent a grayscale value in the
    range 0...1. */
    rsColorRGBA(float grayValue);


    /** \name Setup */

    /** Assigns this color to an avarage value between its original value and the passed color. */
    void avarageWith(rsColorRGBA c)
    {
      r = (rsUint8) (((rsUint16)r + (rsUint16)c.r) >> 1);
      g = (rsUint8) (((rsUint16)g + (rsUint16)c.g) >> 1);
      b = (rsUint8) (((rsUint16)b + (rsUint16)c.b) >> 1);
      a = (rsUint8) (((rsUint16)a + (rsUint16)c.a) >> 1);
    }

    /** Blends this color with the color given in c based on the opacity value of c. The resulting
    color will have r,g,b component that are weighted averages of this color's (old) components and
    the corresponding components of c. The opacity of this color will remain as is. The function is
    useful for compositing images with alpha-blending. */
    void blendWith(rsColorRGBA c)
    {
      rsUint16 w2 = c.a;      // weight for other color's components
      rsUint16 w1 = 256 - w2; // weight for this color's components
      r = (rsUint8) (((rsUint16)r*w1 + (rsUint16)c.r*w2) >> 8);
      g = (rsUint8) (((rsUint16)g*w1 + (rsUint16)c.g*w2) >> 8);
      b = (rsUint8) (((rsUint16)b*w1 + (rsUint16)c.b*w2) >> 8);
    }

    /** Sets the component values of this color as a weighted average of two other colors. The
    weight should be between 0 and 256 where 0 means fully c1 and 256 means fully c2.
    Specifically, we use w as weight for c2 and 255-w as weight for c1 and renormalize
    afterwards. */
    void setAsWeightedAverage(const rsColorRGBA &c1, const rsColorRGBA &c2, rsUint16 w)
    {
      rsUint16 w1 = 256 - w;
      //rsUint16 w1 = 255 - w;

      r = (rsUint8) (((rsUint16)c1.r*w1 + (rsUint16)c2.r*w) >> 8);
      g = (rsUint8) (((rsUint16)c1.g*w1 + (rsUint16)c2.g*w) >> 8);
      b = (rsUint8) (((rsUint16)c1.b*w1 + (rsUint16)c2.b*w) >> 8);
      a = (rsUint8) (((rsUint16)c1.a*w1 + (rsUint16)c2.a*w) >> 8);

      // \todo: check the case for w = 255, with black/white seems like the components have a value
      // of 254 when it should be 255 - maybe this is due to upward/downward rounding when shifting
      // the bits (it should round upward but rounds downward or something)
    }

    // \todo: maybe inline these functions later (they are supposed to be applied per-pixel)

    // Converts the color into a grayscale value by ...
    //void convertToGrayScale(bool useColorWeighting);

    // \todo write functions to combine colors in different ways (weighted average,
    // alpha-blend, etc.)


    /** \name Related Color Creation */

    /** Returns a color that has the same r, g, b components as "this" color and the given alpha
    value. */
    rsColorRGBA withAlpha(rsUint8 alpha)
    {
      rsColorRGBA c = *this;
      c.a = alpha;
      return c;
    }


    /** \name Factory Functions */

    /** Creates a grayscale values where all the r, g, b components are equal to the passes
    grayValue. Optionally an opacity (alpha) value may be passed - if none is passed, a fully
    opaque color will be created. */
    static rsColorRGBA fromGrayValue(rsUint8 grayValue, rsUint8 alpha = 255)
    {
      rsColorRGBA c;
      c.r = c.g = c.b = grayValue;
      c.a = alpha;
      return c;
    }

    // fromInt32, fromString, etc...


    /** \name Data */

    #if defined RS_SYSTEM_WIN32
      rsUint8 b, g, r, a;  // use rsUint8
    #elif defined RS_SYSTEM_LINUX
      rsUint8 a, r, g, b;
    #elif defined RS_SYSTEM_OSX
      rsUint8 b, g, r, a;
    #else
      #error Platform not supported
    #endif

  };

  //===============================================================================================

  /**

  This class represents a color with components red, green, blue and an opacity (aka alpha) value
  using a single precision floating point number for each component.

  */

  class RSLib_API rsColorFloatRGBA
  {

  public:

    /** \name Construction/Destruction */

    /** Default constructor. You may initialize the components as you like or accept the default
    values which are all zero, defining a transparent black as default color. */
    rsColorFloatRGBA(float red = 0.f, float green = 0.f,float blue = 0.f, float alpha = 0.f);


    /** \name Setup */

    /** Assigns this color to an avarage value between its original value and the passed color. */
    inline void avarageWith(rsColorFloatRGBA colorToAverageWith)
    {
      r = 0.5f * (r + colorToAverageWith.r);
      g = 0.5f * (g + colorToAverageWith.g);
      b = 0.5f * (b + colorToAverageWith.b);
      a = 0.5f * (a + colorToAverageWith.a);
    }

    /** Sets the components of thsi color. */
    inline void setComponents(float newRed, float newGreen, float newBlue, float newAlpha)
    {
      r = newRed;
      g = newGreen;
      b = newBlue;
      a = newAlpha;
    }

    // Converts the color into a grayscale value by ...
    //void convertToGrayScale(bool useColorWeighting);

    // \todo write functions to combine colors in different ways (weighted average, alpha-blend,
    // etc.)


    /** \name Inquiry */

    inline float getRed()   const { return r; }
    inline float getGreen() const { return g; }
    inline float getBlue()  const { return b; }
    inline float getAlpha() const { return a; }


    /** \name Static Member Functions */

    static rsColorFloatRGBA getAverageColor(const rsColorFloatRGBA &c1, const rsColorFloatRGBA &c2)
    {
      return rsColorFloatRGBA(0.5f*(c1.r+c2.r),
                              0.5f*(c1.g+c2.g),
                              0.5f*(c1.b+c2.b),
                              0.5f*(c1.a+c2.a));
    }

    /** Returns a color that is the weighted average between the two given colors. c1 is weighted
    with 1-weight and c2 is weighted with weight itself. */
    static rsColorFloatRGBA getWeightedAverageColor(const rsColorFloatRGBA &c1,
                                                    const rsColorFloatRGBA &c2,
                                                    const float weight)
    {
      float w1 = 1.f - weight;
      float w2 = weight;
      return rsColorFloatRGBA(w1*c1.r+w2*c2.r,
                              w1*c1.g+w2*c2.g,
                              w1*c1.b+w2*c2.b,
                              w1*c1.a+w2*c2.a);
    }


    /** \name Data */

    float r, g, b, a;

  };

  // todo: write classes for other color spaces like HLSA, etc.

}

#endif
