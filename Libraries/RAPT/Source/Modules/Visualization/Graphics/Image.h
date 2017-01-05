#ifndef RAPT_IMAGE_H_INCLUDED
#define RAPT_IMAGE_H_INCLUDED

/** This is a class for representing images. An image is defined as an rectangular array of pixels
("picture elements"), each of which having a specified color. The class is templated on the type of 
the color, so it can be used to represent monochrome, RGB, RGBA, HLS, HLSA and whatever types of 
images. 

\todo: 
move some more functions into cpp file
implement access operators (x, y), [row][column] */

template<class TPix>  // pixel type
class Image
{

public:

  /** \name Construction/Destruction */

  /** Constructor. Allocates memory for the pixels. */
  Image(int initialWidth, int initialHeight);

  /** Constructor. Allocates memory for the pixels and initializes pixel values. */
  Image(int initialWidth, int initialHeight, const TPix &initialPixelColor);

  /** Constructor. Initializes width and height and copies the image-data from initialData into
  our member data area. */
  Image(int initialWidth, int initialHeight, const TPix *initialData);

  /** Copy constructor. Creates a deep copy of the pixel data in this image. */
  Image(const Image& other);

  /** Destructor. Frees dynamically allocated memory. */
  virtual ~Image();


  /** \name Setup */

  /** Sets a new size for the image. The contents of the old image is lost when doing this.
  \todo have an optional boolean parameter to retain the contents of the old image */
  void setSize(int newWidth, int newHeight);


  /** \name Inquiry */

  /** Returns the width of the image in pixels. */
  inline int getWidth() const
  {
    return width;
  }

  /** Returns the height of the image in pixels. */
  inline int getHeight() const
  {
    return height;
  }

  /** Returns the number of pixels in the image. */
  inline int getNumPixels() const
  {
    return getWidth() * getHeight();
  }

  /** Returns the size of the whole image in bytes. */
  inline int getByteSize() const
  {
    return width*height*sizeof(TPix); // use height*getLineStride()
  }

  /** Returns the number of data-bytes that each line contains (exluding alignment/fill bytes,
  if any). */
  inline int getLineStrideInBytes() const
  {
    return width*sizeof(TPix);
  }

  /** Returns a pointer to the pixel at the specified location. */
  inline TPix* getPointerToPixel(int x, int y) const
  {
    return &data[y*width+x];
  }

  /** Returns the color of the pixel at (x,y). */
  inline TPix getPixelColor(int x, int y) const
  {
    return data[y*width+x];
  }

  /** Returns true when (x,y) represent valid pixel-coordinates in this image. Mainly for
  debug purposes. */
  inline bool arePixelCoordinatesValid(int x, int y)
  {
    return rsIsInRange(x, 0, width-1) && rsIsInRange(y, 0, height-1);
  }

  /** \name Manipulations */

  /** Sets the color of the pixel at (x,y). */
  inline void setPixelColor(int x, int y, const TPix &newPixelColor)
  {
    rsAssert(rsIsInRange(x, 0, width-1)); // x-coordinate out of range
    rsAssert(rsIsInRange(y, 0, height-1)); // y-coordinate out of range
    data[y*width+x] = newPixelColor;
  }

  /** Fills the whole picture with a solid color. */
  void fillAll(const TPix &colorToFillWith);

  /** Flips the image vertically such that top becomes bottom and vice versa. */
  //void flipTopForBottom();

  // increaseContrast, applyFilter(const ImageFilter &filterToApply), etc...
  // maybe make a class ImageProcessor that takes a pointer to an image and performs the
  // operations on that pointed image

protected:

  /** \name Memory Management */

  /** Allocates the memory for the picture. */
  inline void allocateMemory()
  {
    data = new TPix[width*height];
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

#endif