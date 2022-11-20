/*
template<class TPix>
rsImage<TPix>::rsImage(int initialWidth, int initialHeight)
  : width(initialWidth), height(initialHeight)
{
  allocateMemory();
}

template<class TPix>
rsImage<TPix>::rsImage(int initialWidth, int initialHeight, const TPix &initialPixelColor)
  : width(initialWidth), height(initialHeight)
{
  allocateMemory();
  fillAll(initialPixelColor);
}

template<class TPix>
rsImage<TPix>::rsImage(int initialWidth, int initialHeight, const TPix *initialData)
  : width(initialWidth), height(initialHeight)
{
  allocateMemory();
  memcpy(data, initialData, getByteSize());
}

template<class TPix>
rsImage<TPix>::rsImage(const rsImage<TPix>& other)
{
  width  = other.width;
  height = other.height;
  allocateMemory();
  memcpy(data, other.data, getByteSize());
}

template<class TPix>
rsImage<TPix>::~rsImage()
{
  freeMemory();
}
*/

// setup:

/*
template<class TPix>
void rsImage<TPix>::setSize(int newWidth, int newHeight)
{
  if(width != newWidth || height != newHeight)
  {
    width  = newWidth;
    height = newHeight;
    freeMemory();
    allocateMemory();
  }
}
*/

/*
template<class TPix>
void rsImage<TPix>::fillAll(const TPix &colorToFillWith)
{
  for(int y=0; y<getHeight(); y++)
  {
    for(int x=0; x<getWidth(); x++)
      setPixelColor(x, y, colorToFillWith);
  }
  // ...maybe optimize using memset
}
*/

//template<class TPix>
//void Image<TPix>::flipTopForBottom()
//{
//  int bytesPerLine = getLineStrideInBytes();
//  void *tmpLine = malloc(bytesPerLine);
//  for(int i = 0; i < height/2; i++)
//    swapDataBuffers(&data[i*width], &data[(height-i-1)*width], tmpLine, bytesPerLine);
//  free(tmpLine);
//}


template<class TPix>
void rsImage<TPix>::blendWith(const rsImage<TPix>& img2, TPix w1, TPix w2)
{
  int w = rsMin(getWidth(),  img2.getWidth());
  int h = rsMin(getHeight(), img2.getHeight());
  for(int y = 0; y < h; y++)
    for(int x = 0; x < w; x++)
      (*this)(x, y) = w1 * (*this)(x, y) + w2 * img2(x, y);

  // ToDo:
  // -Maybe use an optimized version of the code when both images have the same dimensions. In this
  //  case addressing the individual pixels via img2(x, y) etc is not necessary
  // -Needs unit test.
}

//-------------------------------------------------------------------------------------------------

template<class TPix>
rsImageResizable<TPix>::rsImageResizable(int initialWidth, int initialHeight)
  : rsImage<TPix>(initialWidth, initialHeight)
{
  maxWidth  = initialWidth;
  maxHeight = initialHeight;
//  maxWidth  = width;
//  maxHeight = height;
}

template<class TPix>
void rsImageResizable<TPix>::setSize(int newWidth, int newHeight)
{
  if(newWidth > maxWidth || newHeight > maxHeight)  // memory reallocation, only if necessary
    setMaxSize(rsMax(maxWidth, newWidth), rsMax(maxHeight, newHeight));
//  width  = newWidth;
//  height = newHeight;
  rsImage<TPix>::width = newWidth;
  rsImage<TPix>::height = newHeight;
}

template<class TPix>
void rsImageResizable<TPix>::setMaxSize(int newMaxWidth, int newMaxHeight)
{
  maxWidth  = newMaxWidth;
  maxHeight = newMaxHeight;
  rsImage<TPix>::setSize(maxWidth, maxHeight); // for memory reallocation
}
