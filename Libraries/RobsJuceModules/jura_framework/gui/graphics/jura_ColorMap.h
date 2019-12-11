#ifndef jura_ColorMap_h
#define jura_ColorMap_h

/** A class for representing a colormap that is used for visualizing 2D data. It is based on a 
juce::ColourGradient object, but allows for (much!) faster access of the individual colors in the
gradient by creating an internal array of pre-interpolated colors. So, basically, it's a lookup 
table for a ColourGradient object. It derives from ChangeBroadcaster and sends out a 
changeListenerCallback whenever it has changed, so objects that need to keep track of such changes
can derive from ChangeListener and register themselves to receive callbacks.

\todo: maybe get rid of specifying default-colormaps that are not based on a juce::ColourGradient
or: make an own ColorGradientAHSL class based on ColorAHSL and in the xml sve/load function make a
distinction - maybe the xml format must specify the color-format (ARGB, AHSL, maybe others)

References:
http://www.research.ibm.com/people/l/lloydt/color/color.HTM
https://www.mathworks.com/tagteam/81137_92238v00_RainbowColorMap_57312.pdf
http://davidjohnstone.net/pages/lch-lab-colour-gradient-picker

todo: include the colormaps desribed here:
https://bids.github.io/colormap/
especially parula, magma, viridis, inferno are attractive

*/

class JUCE_API ColorMap : public ChangeBroadcaster
{

public:

  /** Enumeration of the hard coded predefined maps. */
  enum defaultMaps
  {
    gray,
    fire,
    ice, 
    rainbow
  };
  // \todo: provide bipolar maps - for example going from bright blue to black in the middle and 
  // then to bright red

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

  /** Sets up the colormap for the given xml element. */ 
  void setFromXml(const XmlElement& xml);

  /** Returns an xml representation of this colormap. */
  XmlElement* getAsXml();


  /** \name Color retrieval */

  /** Returns the color as uint32 value at the given normalized position which must be in the range 
  0..1 (ends inclusive). For efficiency, the range is not checked, so the caller must make sure 
  that it is indeed in this range, otherwise an access violation will occur. */
  inline uint32 getColorAsUint32(float position) const
  {
    return colors[int(position*lastIndex)];
  }


  /** \name Misc */

  /** Returns one of the default color gradients as juce::ColourGradient object. */
  ColourGradient getDefaultGradient(int index);

  /** Returns the juce::ColourGradient object that this color map is based on. */
  ColourGradient getAsColourGradient() { return gradient; }

protected:

  /** Fills our array of default map names. Called internally in the constructor. */
  void fillDefaultMapNameArray();

  /** Updates our internally stored array of colors according to the gradient member. */
  void updateArray();

  std::vector<uint32> colors;
  int lastIndex;

  int defaultMapIndex;          // if not -1, we are using one of the predefined maps
  StringArray defaultMapNames;
    // this could actually be a static member - we would have to make fillDefaultMapNameArray
    // a static member function then and inside it, we could check, if the array is already filled
    // and fill it only if necessarry

  ColourGradient gradient; // the juce color gradient object on which this map is based


  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ColorMap)
};

#endif  
