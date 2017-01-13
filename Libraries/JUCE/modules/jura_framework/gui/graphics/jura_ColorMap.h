#ifndef jura_ColorMap_h
#define jura_ColorMap_h

/** A class for representing a colormap that is used for visualizing 2D data. It is based on a 
juce::ColourGradient object, but allows for (much!) faster access of the individual colors in the
gradient by creating an internal array of pre-interpolated colors. So, basically, it's a lookup 
table for a ColourGradient object. 

\todo: provide methods for conversion from/to xml
*/

class JUCE_API ColorMap
{

public:

  /** Enumeration of the hard coded predefined maps. */
  enum defaultMaps
  {
    grayScale,
    fire,
    ice
  };

  /** Constructor. */
  ColorMap();


  /** \name Setup */

  /** Fills the array of colors based on the passed ColourGradient object. */
  void setFromColourGradient(const ColourGradient &g);

  /** Fills the array of colors with one of the default, hard coded maps. 
  @see defaultMaps */
  void setDefaultMap(int index);

  /** Sets the size of the array. Larger sizes will give smoother color transitions. */
  void setSize(int newSize);


  /** \name Color retrieval */

  /** Returns the color as uint32 value at the given normalized index which must be in the range 
  0..1 (ends inclusive). For efficiency, the range is not checked, so the caller must make sure 
  that it is indeed in this range, otherwise an access violation will occur. */
  inline uint32 getColorAsUint32(float normalizedIndex)
  {
    return colors[int(normalizedIndex*lastIndex)];
  }


  /** \name Misc */

  /** Returns one of the default color gradients as juce::ColourGradient object. */
  ColourGradient getDefaultGradient(int index);

protected:

  /** Updates our internally stored array of colors according to the gradient member. */
  void updateArray();


  std::vector<uint32> colors;
  int lastIndex;

  ColourGradient gradient; // the juce color gradient object on which this map is based


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ColorMap)
};

#endif  
