
template<class TPix, class TWgt, class TCor>
ImageDrawer<TPix, TWgt, TCor>::ImageDrawer(Image<TPix> *imageToDrawOn)
{
  setImageToDrawOn(imageToDrawOn);

  blendMode = 0;
  blendFunction = linearBlend;
}

// setup:

template<class TPix, class TWgt, class TCor>
void ImageDrawer<TPix, TWgt, TCor>::setImageToDrawOn(Image<TPix> *imageToDrawOn)
{
  image = imageToDrawOn;
}

template<class TPix, class TWgt, class TCor>
void ImageDrawer<TPix, TWgt, TCor>::setBlendMode(int newMode)
{
  rsAssert(newMode >= 0);
  rsAssert(newMode < NUM_BLEND_MODES);
  blendMode = newMode;
  switch(blendMode)
  {
  case BLEND_LINEAR:       blendFunction = linearBlend;    break;
  case BLEND_ADD_CLIP:     blendFunction = addAndClip;     break;
  case BLEND_ADD_SATURATE: blendFunction = addAndSaturate; break;
  }
}

// blend functions:

template<class TPix, class TWgt, class TCor>
void ImageDrawer<TPix, TWgt, TCor>::linearBlend(TPix &pixel, TPix color, TWgt blend)
{
  pixel = (1-blend)*pixel + blend*color;
}

template<class TPix, class TWgt, class TCor>
void ImageDrawer<TPix, TWgt, TCor>::addAndClip(TPix &pixel, TPix color, TWgt blend)
{
  pixel = rsMin(TPix(1), pixel + blend*color);
}

template<class TPix, class TWgt, class TCor>
void ImageDrawer<TPix, TWgt, TCor>::addAndSaturate(TPix &pixel, TPix color, TWgt blend)
{
  color *= blend;
  pixel  = (pixel + color) / (TPix(1) + color);
}
