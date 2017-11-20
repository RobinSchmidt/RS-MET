using namespace RSLib;

// Construction/Destruction:

rsImageRegionRGBA::rsImageRegionRGBA(rsImageRGBA *image, int x_, int y_, int width, int height)
: theImage(image)
, rsRectangle2D<int>(x_, y_, width, height)
{
  clipRegionToUnderlyingImage();
}

// Misc:

void rsImageRegionRGBA::clipRegionToUnderlyingImage()
{
  x = rsClipToRange(x, 0, theImage->getWidth()-1);
  w = rsClipToRange(w, 0, theImage->getWidth()-x);
  y = rsClipToRange(y, 0, theImage->getHeight()-1);
  h = rsClipToRange(h, 0, theImage->getHeight()-y);
}