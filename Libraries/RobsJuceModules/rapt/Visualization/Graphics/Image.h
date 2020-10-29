#ifndef RAPT_IMAGE_H_INCLUDED
#define RAPT_IMAGE_H_INCLUDED

/** This is a class for representing images. An image is defined as an rectangular array of pixels
("picture elements"), each of which having a specified color. The class is templated on the type of 
the color, so it can be used to represent monochrome, RGB, RGBA, HLS, HLSA and whatever types of 
images. 

\todo: 
-move some more functions into cpp file
-implement access operators (x, y), [row][column] 
-maybe factor out a common baseclass rsArray2D from rsImage and rsMatrix
*/

template<class TPix>  // pixel type
class rsImage
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Construction/Destruction */

  /** Constructor. Allocates memory for the pixels. */
  rsImage(int initialWidth = 1, int initialHeight = 1) 
    : width(initialWidth), height(initialHeight)
  {
    allocateMemory();
  }

  /** Constructor. Allocates memory for the pixels and initializes pixel values. */
  rsImage(int initialWidth, int initialHeight, const TPix &initialPixelColor)
    : width(initialWidth), height(initialHeight)
  {
    allocateMemory();
    fillAll(initialPixelColor);
  }

  /** Constructor. Initializes width and height and copies the image-data from initialData into
  our member data area. */
  rsImage(int initialWidth, int initialHeight, const TPix *initialData)
    : width(initialWidth), height(initialHeight)
  {
    allocateMemory();
    memcpy(data, initialData, getByteSize());
  }

  /** Copy constructor. Creates a deep copy of the pixel data in this image. */
  rsImage(const rsImage& other)
  {
    width  = other.width;
    height = other.height;
    allocateMemory();
    memcpy(data, other.data, getByteSize());
  }

  /** Destructor. Frees dynamically allocated memory. */
  virtual ~rsImage()
  {
    freeMemory();
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets a new size for the image. The contents of the old image is lost when doing this.
  \todo have an optional boolean parameter to retain the contents of the old image */
  virtual void setSize(int newWidth, int newHeight)
  {
    if(width != newWidth || height != newHeight)
    {
      width  = newWidth;
      height = newHeight;
      freeMemory();
      allocateMemory();
    }
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the width of the image in pixels. */
  inline int getWidth() const { return width; }

  /** Returns the height of the image in pixels. */
  inline int getHeight() const { return height; }

  /** Returns the number of pixels in the image. */
  inline int getNumPixels() const { return getWidth() * getHeight(); }

  /** Returns the size of the whole image in bytes. */
  inline int getByteSize() const { return width*height*sizeof(TPix); }
    // maybe use height*getLineStride()

  /** Returns the number of data-bytes that each line contains (exluding alignment/fill bytes,
  if any). */
  inline int getLineStrideInBytes() const { return width*sizeof(TPix); }

  /** Returns a pointer to the pixel at the specified location. */
  inline TPix* getPixelPointer(int x, int y) const { return &data[y*width+x]; }
  // maybe use operator (x, y) instead of function..., but the operator should return a reference

  /** Returns the color of the pixel at (x,y). */
  inline TPix getPixelColor(int x, int y) const { return data[y*width+x]; }

  /** Returns true when (x,y) represent valid pixel-coordinates in this image. Mainly for
  debug purposes. */
  inline bool arePixelCoordinatesValid(int x, int y)
  { return rsIsInRange(x, 0, width-1) && rsIsInRange(y, 0, height-1); }

  /** Returns true, iff this image as the given shape. */
  inline bool hasShape(int width, int height) const
  {
    return this->width == width && this->height == height;
  }

  /** Returns true, iff this image has the same shape as the given other image. */
  bool hasSameShapeAs(const rsImage<TPix>& other) const
  {
    return this->hasShape(other.width, other.height);
  }

  /** Compares all pixel values of this image to those of another image and returns true, if they 
  are all equal up to some given tolerance. It assumes that the other image has the same width
  and height as this. */
  inline bool areAllPixelsEqualTo(const rsImage<TPix>* otherImage, TPix tolerance = TPix(0))
  {
    rsAssert(width  == otherImage->width);
    rsAssert(height == otherImage->height);
    return rsArrayTools::almostEqual(data, otherImage->data, width*height, tolerance);

    /*
    // old:
    TPix err = rsArrayTools::maxDeviation(data, otherImage->data, width*height);
    return rsAbs(err) <= tolerance;
    //return abs(err) <= tolerance;
    // this implementation is not suitable for TPix = rsPixelRGB (in private repo) - it's made for 
    // (signed) floating point numbers as pixels
    // maybe it has to be done in terms of an isCloseTo function that can be implemented also for 
    // rsPixelRGB
    */
  }




  /** Converts the image to a flat array of type std::vector. */
  std::vector<TPix> toStdVector() { return toVector(data, width*height); }

  //-----------------------------------------------------------------------------------------------
  /** \name Manipulations */

  /** Sets the color of the pixel at (x,y). */
  inline void setPixelColor(int x, int y, const TPix &newPixelColor)
  {
    rsAssert(rsIsInRange(x, 0, width-1)); // x-coordinate out of range
    rsAssert(rsIsInRange(y, 0, height-1)); // y-coordinate out of range
    data[y*width+x] = newPixelColor;
  }

  /** Fills the whole picture with a solid color. */
  void fillAll(const TPix &colorToFillWith)
  {
    for(int y = 0; y < getHeight(); y++)
    {
      for(int x = 0; x < getWidth(); x++)
        setPixelColor(x, y, colorToFillWith);
    }
    // ...maybe optimize using memset
  }

  /** Clears the image by setting all pixels to the given color. */
  inline void clear(TPix color = TPix(0)) { fillAll(color);}

  /** Copies the pixel data form the source image into this image. Assumes that the source has the 
  same shape as this image. */
  inline void copyPixelDataFrom(const rsImage<TPix>& source) 
  { 
    rsAssert(this->hasSameShapeAs(source));
    rsArrayTools::copy(source.data, this->data, getNumPixels());
  }
  // todo: maybe if the source has a different shape, change the shape of this image





  template<class TPix2>
  inline void convertPixelDataFrom(const rsImage<TPix2>& source) 
  { 
    rsAssert(source.getWidth() == width && source.getHeight() == height);
    rsArrayTools::convert(source.getPixelPointer(0,0), this->data, getNumPixels());
  }


  /** Flips the image vertically such that top becomes bottom and vice versa. */
  //void flipTopForBottom();

  // increaseContrast, applyFilter(const ImageFilter &filterToApply), etc...
  // maybe make a class ImageProcessor that takes a pointer to an image and performs the
  // operations on that pointed image



  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Allows read/write acces to the pixel at position x,y. */
  inline TPix& operator()(int x, int y) { return data[y*width+x]; }

  inline const TPix& operator()(int x, int y) const { return data[y*width+x]; }




  /** Returns true, iff the two given images have the same shape. */
  static bool haveSameShape(const rsImage<TPix>& img1, const rsImage<TPix>& img2)
  { return img1.getWidth() == img2.getWidth() && img1.getHeight() == img2.getHeight(); }
  // obsolete: use img1.hasSameShapeAs(img2) instead

  /** Returns true, iff the three given images have the same shape. */
  static bool haveSameShape(
    const rsImage<TPix>& img1, const rsImage<TPix>& img2, const rsImage<TPix>& img3)
  { return haveSameShape(img1, img2) && haveSameShape(img1, img3); }


protected:

  /** \name Memory Management */

  /** Allocates the memory for the picture. */
  inline void allocateMemory()
  {
    data = new TPix[width*height];
    fillAll(TPix(0));
  }

  /** Frees the allocated memory for the picture. */
  inline void freeMemory()
  {
    delete[] data;
  }

  // data members:

  int width, height;
  TPix *data;         // contiguous block of memory for the whole image

};

//=================================================================================================

/** A subclass of Image that allows for more efficient resizing. The regular image class does a 
memory reallocation, everytime you call setSize. This subclass here defines a maximum width and 
height below which a resizing operation does not lead to memory reallocation and overrides setSize 
accordingly. This is useful for images that must be dynamically resized at runtime often. */

template<class TPix>  // pixel type
class rsImageResizable : public rsImage<TPix>
{

public:

  /** \name Construction/Destruction */

  /** Constructor. Allocates memory for the pixels. */
  rsImageResizable(int initialWidth = 1, int initialHeight = 1);


  /** \name Setup */

  /** Sets a new size for the image. This will cause a memory reallocation only if the new width 
  or height are greater than our respective maximum width or height values. */
  virtual void setSize(int newWidth, int newHeight) override;

  /** Sets a new maximum width and height. */
  void setMaxSize(int newMaxWidth, int newMaxHeight);


  /** \name Inquiry */

  /** Returns maximum width (without re-allocation of memory taking place). */
  int getMaxWidth() { return maxWidth; }

  /** Returns maximum height (without re-allocation of memory taking place). */
  int getMaxHeight() { return maxHeight; }


protected:

  int maxWidth, maxHeight;

};

#endif