template<class TPix, class TWgt, class TCor>
ImagePainter<TPix, TWgt, TCor>::ImagePainter(Image<TPix> *imageToPaintOn, Image<TWgt> *brushToUse)
{
  image = imageToPaintOn;
  brush = brushToUse;
}

// setup

template<class TPix, class TWgt, class TCor>
void ImagePainter<TPix, TWgt, TCor>::setImageToPaintOn(Image<TPix> *imageToPaintOn)
{
  image = imageToPaintOn;
}

template<class TPix, class TWgt, class TCor>
void ImagePainter<TPix, TWgt, TCor>::setBrushToUse(Image<TWgt> *brushToUse)
{
  brush = brushToUse;
}

// painting

template<class TPix, class TWgt, class TCor>
void ImagePainter<TPix, TWgt, TCor>::paintDot(TCor x, TCor y, TPix color)
{
  int wi = image->getWidth();
  int hi = image->getHeight();
  int wb = brush->getWidth();
  int hb = brush->getHeight();

  // ...something to do...
  // we need a (nested) loop over all (x,y) pixels in the brush, multiply the brush value there 
  // with the color accumulate it into the corresponding location in the target image using 
  // bilinear deinterpolation to write into fractional positions
  // or: loop over a rectangle of size wb,wh in the target image and read out the brush at 
  // corresponding positions using bilinear interpolation ...maybe provide both methods so we
  // can compare the results (should they be equal? idk) and also the performance
}