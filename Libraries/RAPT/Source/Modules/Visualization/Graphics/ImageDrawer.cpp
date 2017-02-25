template<class TPix, class TWgt, class TCor>
ImageDrawer<TPix, TWgt, TCor>::ImageDrawer(Image<TPix> *imageToDrawOn)
{
  setImageToDrawOn(imageToDrawOn);
  setBlendMode(0);
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

//=================================================================================================

template<class TPix, class TWgt, class TCor>
void LineDrawer<TPix, TWgt, TCor>::setLineProfile(int newProfile)
{
  rsAssert(newProfile >= 0);
  rsAssert(newProfile < NUM_LINE_PROFILES);
  profileIndex = newProfile;

  // include a switch and assign the lineProfile function pointer here - provide a couple of
  // predefined profiles in a similar way as is done the the blend-modes

}
