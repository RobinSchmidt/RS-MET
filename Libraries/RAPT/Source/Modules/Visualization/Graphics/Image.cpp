template<class TPix>
Image<TPix>::Image(int initialWidth, int initialHeight) 
  : width(initialWidth), height(initialHeight)
{
  allocateMemory();
}

template<class TPix>
Image<TPix>::Image(int initialWidth, int initialHeight, const TPix &initialPixelColor)
  : width(initialWidth), height(initialHeight)
{
  allocateMemory();
  fillAll(initialPixelColor);
}

template<class TPix>
Image<TPix>::Image(int initialWidth, int initialHeight, const TPix *initialData)
  : width(initialWidth), height(initialHeight)
{
  allocateMemory();
  memcpy(data, initialData, getByteSize());
}

template<class TPix>
Image<TPix>::Image(const Image<TPix>& other)
{
  width  = other.width;
  height = other.height;
  allocateMemory();
  memcpy(data, other.data, getByteSize());
}

template<class TPix>
Image<TPix>::~Image()
{
  freeMemory();
}

// setup:

template<class TPix>
void Image<TPix>::setSize(int newWidth, int newHeight)
{
  if(width != newWidth || height != newHeight)
  {
    width  = newWidth;
    height = newHeight;
    freeMemory();
    allocateMemory();
  }
}

template<class TPix>
void Image<TPix>::fillAll(const TPix &colorToFillWith)
{
  for(int y=0; y<getHeight(); y++)
  {
    for(int x=0; x<getWidth(); x++)
      setPixelColor(x, y, colorToFillWith);
  }
  // ...maybe optimize using memset
}

//template<class TPix>
//void Image<TPix>::flipTopForBottom()
//{
//  int bytesPerLine = getLineStrideInBytes();
//  void *tmpLine = malloc(bytesPerLine);
//  for(int i = 0; i < height/2; i++)
//    swapDataBuffers(&data[i*width], &data[(height-i-1)*width], tmpLine, bytesPerLine);
//  free(tmpLine);
//}

//-------------------------------------------------------------------------------------------------

template<class TPix>
ImageResizable<TPix>::ImageResizable(int initialWidth, int initialHeight)
  : Image<TPix>(initialWidth, initialHeight)
{
  maxWidth  = width;
  maxHeight = height;
}

template<class TPix>
void ImageResizable<TPix>::setSize(int newWidth, int newHeight)
{
  if(newWidth > maxWidth || newHeight > maxHeight)  // memory reallocation, only if necessary
    setMaxSize(rsMax(maxWidth, newWidth), rsMax(maxHeight, newHeight));
  width  = newWidth;
  height = newHeight;
}

template<class TPix>
void ImageResizable<TPix>::setMaxSize(int newMaxWidth, int newMaxHeight)
{
  maxWidth  = newMaxWidth;
  maxHeight = newMaxHeight;
  Image<TPix>::setSize(maxWidth, maxHeight); // for memory reallocation
}