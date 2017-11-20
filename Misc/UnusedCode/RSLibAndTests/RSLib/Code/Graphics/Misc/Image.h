#ifndef RS_IMAGE_H
#define RS_IMAGE_H

namespace RSLib
{

  /**

  This is a class for representing images. An image is defined as an rectangular array of pixels
  ("picture elements"), each of which having a specified color. The class is templated on the type
  of the color, so it can be used to represent monochrome, RGB, RGBA, HLS, HLSA and whatever
  types of images. Explicit instantiations for the most common ones (monochrome, RGBA) are provided
  by the classes rsImageGray and rsImageRGBA.

  */

  template<class ColorType>
  class rsImage
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. Allocates memory for the pixels.
    */
    rsImage(int initialWidth, int initialHeight)
      : width(initialWidth)
      , height(initialHeight)

    {
      allocateMemory();
    }

    /** Constructor. Allocates memory for the pixels and optionally clears  
    \todo: remove boolean parameter - init always all pixels
    */
    rsImage(int initialWidth, int initialHeight, bool initPixelColors,
      const ColorType &initialPixelColor)
      : width(initialWidth)
      , height(initialHeight)

    {
      allocateMemory();
      if( initPixelColors == true )
        fillAll(initialPixelColor);
    }

    /** Constructor. Initializes width and height and copies the image-data from initialData into
    our member data area. */
    rsImage(int initialWidth, int initialHeight, const ColorType *initialData)
      : width(initialWidth)
      , height(initialHeight)
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

    /** \name Setup */

    /** Sets a new size for the image. The contents of the old image is lost when doing this.
    \todo have an optional boolean parameter to retain the contents of the old image */
    void setSize(int newWidth, int newHeight)
    {
      if( width != newWidth || height != newHeight )
      {
        width  = newWidth;
        height = newHeight;
        freeMemory();
        allocateMemory();
      }
    }

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
    int getByteSize() const
    {
      return width*height*sizeof(ColorType); // use height*getLineStride()
    }

    /** Returns the number of data-bytes that each line contains (exluding alignment/fill bytes,
    if any). */
    int getLineStrideInBytes() const
    {
      return width*sizeof(ColorType);
    }

    /** Returns a pointer to the pixel at the specified location. */
    inline ColorType* getPointerToPixel(int x, int y) const
    {
      return &data[y*width+x];
    }

    /** Returns the color of the pixel at (x,y). */
    inline ColorType getPixelColor(int x, int y) const
    {
      return data[y*width+x];
    }

    /** Returns true when (x,y) represent valid pixel-coordinates in this image. Mainly for
    debug purposes. */
    bool arePixelCoordinatesValid(int x, int y)
    {
      return rsIsInRange(x, 0, width-1) && rsIsInRange(y, 0, height-1);
    }

    /** \name Manipulations */

    /** Sets the color of the pixel at (x,y). */
    inline void setPixelColor(int x, int y, const ColorType &newPixelColor)
    {
      rsAssert( rsIsInRange(x, 0, width-1) ); // x-coordinate out of range
      rsAssert( rsIsInRange(y, 0, height-1)); // y-coordinate out of range
      data[y*width+x] = newPixelColor;
    }

    /** Fills the whole picture with a solid color. */
    void fillAll(const ColorType &colorToFillWith)
    {
      for(int y=0; y<getHeight(); y++)
      {
        for(int x=0; x<getWidth(); x++)
          setPixelColor(x, y, colorToFillWith);
      }
      // ...maybe optimize using memset
    }

    /** Flips the image vertically such that top becomes bottom and vice versa. */
    void flipTopForBottom()
    {
      int bytesPerLine = getLineStrideInBytes();
      void *tmpLine = malloc(bytesPerLine);
      for(int i = 0; i < height/2; i++)
        swapDataBuffers(&data[i*width], &data[(height-i-1)*width], tmpLine, bytesPerLine);
      free(tmpLine);
    }

    // increaseContrast, applyFilter(const ImageFilter &filterToApply), etc...
    // maybe make a class ImageProcessor that takes a pointer to an image and performs the
    // operations on that pointed image

  protected:

    /** \name Memory Management */

    /** Allocates the memory for the picture. */
    void allocateMemory()
    {
      data = new ColorType[width*height];
    }

    /** Frees the allocated memory for the picture. */
    void freeMemory()
    {
      delete[] data;
    }

    // data members:

    int width, height;  // maybe use rsUint32
    ColorType *data;    // contiguous block of memory for the whole image

  };

  // explicit instantiations to be compiled into the dll:
  //RSLib_TEMPLATE_INSTANCE template class RSLib_API rsImage<rsColorRGBA>;
  //RSLib_TEMPLATE_INSTANCE template class RSLib_API rsImage<rsUint8>;

  // typedefs for explicit instantiations:
  typedef rsImage<rsColorRGBA> RSLib_API rsImageRGBA;
  typedef rsImage<rsUint8> RSLib_API rsImageGray;

  // \todo: check if the instantiation is present in the dll only
  // maybe write a macro that facilitates these explicit intstantiations

  //===============================================================================================

  /**

  This is a class for representing a rectangular region of an rsImage object.

  \todo currently this class uses a template instantiation of rsImageRGBA - maybe it can be
  generalized to be itself a template class

  \todo provide getters/setters for x,y,w,h

  */

  class RSLib_API rsImageRegionRGBA : public rsRectangle2D<int>
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. */
    rsImageRegionRGBA(rsImageRGBA *image = NULL, int x = 0, int y = 0,
      int width = 0, int height = 0);


    /** \name Inquiry */

    /** Returns the width of the image in pixels. */
    inline int getWidth() const
    {
      return w;
    }

    /** Returns the height of the image in pixels. */
    inline int getHeight() const
    {
      return h;
    }

    /** Returns the number of pixels in the image region. */
    inline int getNumPixels() const
    {
      return getWidth() * getHeight();
    }

    /** Returns the distance in bytes between successive lines. */
    int getLineStrideInBytes() const
    {
      return theImage->getLineStrideInBytes();
    }

    /** Returns a pointer to the pixel at the specified location with respect to the origin of the
    region. */
    inline rsColorRGBA* getPointerToPixel(int x, int y) const
    {
      return theImage->getPointerToPixel(this->x + x, this->y + y);
    }

    /** Returns the color of the pixel at (x, y). */
    inline rsColorRGBA getPixelColor(int x, int y) const
    {
      return theImage->getPixelColor(this->x + x, this->y + y);
    }


    /** \name Manipulations */

    /** Sets the color of the pixel at (x, y). */
    inline void setPixelColor(int x, int y, const rsColorRGBA &newPixelColor)
    {
      rsAssert( rsIsInRange(x, 0, w-1) ); // x-coordinate out of range
      rsAssert( rsIsInRange(y, 0, h-1));  // y-coordinate out of range
      theImage->setPixelColor(this->x + x, this->y + y, newPixelColor);
    }

    /** Fills the whole picture with a solid color. */
    void fillAll(const rsColorRGBA &colorToFillWith)
    {
      for(int iy = 0; iy < getHeight(); iy++)
      {
        for(int ix = 0; ix < getWidth(); ix++)
          setPixelColor(ix, iy, colorToFillWith);
      }

      // \todo optimize by retrieving a data-pointer and manipulating the data directly,
      // possibly using memset
    }

  protected:

    /** \name Misc */

    /** Makes sure that the region does not extend beyond the bounds of the underlying image */
    void clipRegionToUnderlyingImage();


    /** \name Data */

    rsImageRGBA *theImage;

  };

}

#endif
