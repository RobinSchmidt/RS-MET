using namespace RAPT;

void fillWithCheckerBoardPattern(rsImageF& image)
{
  image.clear();
  for(int x = 0; x < image.getWidth(); x++) {
    for(int y = 0; y < image.getHeight(); y++) {
      if( (x+y) % 2 == 0) // x+y even -> fill pixel
        image(x, y) = 1;
    }
  }
}

void fillWithCross(rsImageF& image)
{
  image.clear();
  int cx = image.getWidth()  / 2;
  int cy = image.getHeight() / 2;
  for(int x = 0; x < image.getWidth(); x++)
    image(x, cy) = 1;
  for(int y = 0; y < image.getHeight(); y++)
    image(cx, y) = 1;
}

bool testAlphaMaskPaintingAntiAliased(rsImagePainterFFF& painter, float x, float y)
{
  // We clear the image and paint a "dot" into the image (that is pointed to in the passed 
  // ImagePainter) and then check, if all pixels outside the rectangle of the dot are still
  // zero after that.
  // todo: maybe check also, if the pixels inside the rectangle have expected values
  // ...but this is complicated, maybe first, we should check if they are nonzero (but this 
  // assumes that the dot mask is nonzero everywhere) - maybe we can do it by providing a prototype
  // paint function that loops through the source pixels in the mask (instead of target pixels in 
  // the image), computes the non-integer position where to put the pixel and calla 
  // plotPixel(float x, float y) function for each target point, which performs bounds checking
  // internally - and then compare the pictures drawn with the regular and prototype function.

  bool result = true;

  // initializations:
  int wi, hi, wm, hm, xs, xe, ys, ye;
  rsImageF     *img = painter.getImage();
  rsAlphaMaskF *msk = painter.getAlphaMask();
  wi = img->getWidth();       // image width
  hi = img->getHeight();      // image height
  wm = msk->getWidth();       // mask width
  hm = msk->getHeight();      // mask height
  xs = (int)floor(x-0.5*wm);  // x start (of possibly nonzero values)
  xe = (int)ceil( x+0.5*wm);  // x end
  ys = (int)floor(y-0.5*hm);  // y start
  ye = (int)ceil( y+0.5*hm);  // y end

  // clear the image and paint the dot:
  img->clear();
  painter.paintDotViaMask(x, y, 1.f);

  // loop through the pixels of the image and check if those which should still be zero are zero 
  // indeed:
  for(int yi = 0; yi < hi; yi++)
  {
    for(int xi = 0; xi < wi; xi++)
    {
      if(xi < xs || xi > xe || yi < ys || yi > ye)
        result &= (*img)(xi, yi) == 0.f;
      else
      {
        //result &= (*img)(xi, yi) != 0.f;
        // whether or not the pixel should be nonzero depends on the mask ...maybe we could figure
        // out the correct target value and check against that...
      }
    }
  }

  return result;
}

bool imagePainterUnitTest()
{
  bool result = true;

  rsImageF image;
  rsAlphaMaskF mask;  // maybe use a regular image as mask
  rsImagePainterFFF painter(&image, &mask);


  int imageWidth  = 30;
  int imageHeight = 30; // 50x50 image with 3x3 mask gives an access violation
  int maskSize    = 7;

  // maybe, we should 1st use the simplest case: 1x1 mask

  image.setSize(imageWidth, imageHeight);
  image.clear();

  mask.setSize(maskSize);
  result &= mask.getWidth()  == maskSize;
  result &= mask.getHeight() == maskSize;
  //fillWithCheckerBoardPattern(mask);
  //fillWithCross(mask);
  mask.fillAll(1.f);  // full white

  //painter.paintDotViaMask(0.25f, 0.75f, 1);
  //painter.paintDotViaMask(2.25f, 3.75f, 1);
  //painter.paintDotViaMask(3.25f, 3.75f, 1);

  // draw in center and at all 4 edges:
  float dx = 0.5;
  float dy = 0.5;
  float w  = (float)imageWidth;
  float w2 = w/2;
  float h  = (float)imageHeight;
  float h2 = h/2;
  float b  = 0.75f; // brightness

  // do unit tets for various cases:
  dx = 0.5;
  dy = 0.5;
  result &= testAlphaMaskPaintingAntiAliased(painter,     dx,     dy);   // top-left
  result &= testAlphaMaskPaintingAntiAliased(painter, w2 +dx,     dy);   // top-center
  result &= testAlphaMaskPaintingAntiAliased(painter, w-1+dx,     dy);   // top-right
  result &= testAlphaMaskPaintingAntiAliased(painter,     dx, h2 +dy);   // center-left
  result &= testAlphaMaskPaintingAntiAliased(painter, w-1+dx, h2 +dy);   // center-right
  result &= testAlphaMaskPaintingAntiAliased(painter,     dx, h-1+dy);   // bottom-left
  result &= testAlphaMaskPaintingAntiAliased(painter, w2 +dx, h-1+dy);   // bottom-center
  result &= testAlphaMaskPaintingAntiAliased(painter, w-1+dx, h-1+dy);   // bottom-right
  result &= testAlphaMaskPaintingAntiAliased(painter, w2 +dx, h2 +dy);   // center
  rsAssert(result);

  // try some cases where the dot is painted outside the image:
  dx = -3.0; // when the mask size is 5, -3.0 still works, but -3.1 doesn't
  dy = -3.0;
  result &= testAlphaMaskPaintingAntiAliased(painter,     dx,     dy);   // top-left
  result &= testAlphaMaskPaintingAntiAliased(painter, w2 +dx,     dy);   // top-center
  result &= testAlphaMaskPaintingAntiAliased(painter, w-1+dx,     dy);   // top-right
  result &= testAlphaMaskPaintingAntiAliased(painter,     dx, h2 +dy);   // center-left
  result &= testAlphaMaskPaintingAntiAliased(painter, w-1+dx, h2 +dy);   // center-right
  result &= testAlphaMaskPaintingAntiAliased(painter,     dx, h-1+dy);   // bottom-left
  result &= testAlphaMaskPaintingAntiAliased(painter, w2 +dx, h-1+dy);   // bottom-center
  result &= testAlphaMaskPaintingAntiAliased(painter, w-1+dx, h-1+dy);   // bottom-right
  result &= testAlphaMaskPaintingAntiAliased(painter, w2 +dx, h2 +dy);   // center
  rsAssert(result);



  // paint some dots and write image to file for visual inspection:
  dx = 0.5;
  dy = 0.5;
  image.clear();
  painter.paintDotViaMask(    dx,     dy, b);   // top-left
  painter.paintDotViaMask(w2 +dx,     dy, b);   // top-center
  painter.paintDotViaMask(w-1+dx,     dy, b);   // top-right
  painter.paintDotViaMask(    dx, h2 +dy, b);   // center-left
  painter.paintDotViaMask(w-1+dx, h2 +dy, b);   // center-right
  painter.paintDotViaMask(    dx, h-1+dy, b);   // bottom-left
  painter.paintDotViaMask(w2 +dx, h-1+dy, b);   // bottom-center
  painter.paintDotViaMask(w-1+dx, h-1+dy, b);   // bottom-right
  painter.paintDotViaMask(w2 +dx, h2 +dy, b);   // center

  //painter.paintDotViaMask(10.2f, 10.6f, 1);
  // we use the 1000 here, because the painter uses this strange saturating function - maybe, we 
  // should introduce a blend-mode: mix, add, add-and-clip, add-and-saturate, multiply, ...
  // for testing here, we should use either alpha-blend or add-and-clip (should give same results)

  //writeImageToFilePPM(mask,  "PaintTestMask.ppm");
  //writeImageToFilePPM(image, "PaintTestImage.ppm");

  // ...

  return result;
}

bool colorUnitTest()
{
  bool ok = true;

  using Real  = float;
  using Color = rsColor<Real>;

  // test back-and-forth conversion between HSL and RGB

  int N = 20;        // number of samples along the 3 dimensions (HSL or RGB etc.)

  Real x, y, z;      // original values
  Real a, b, c;      // converted values
  Real X, Y, Z;      // back-converted values
  Real tol = Real(1.e-5); // We need a quite high tolerance! Are the formulas numerically bad?

  /*
  x = 0.0f;
  y = 0.1f;
  z = 0.0f;
  Color::hsl2rgb(x, y, z, &a, &b, &c);
  Color::rgb2hsl(a, b, c, &X, &Y, &Z);
  //Color::rgb2hsl(x, y, z, &a, &b, &c);
  */

  auto check = [](Real x, Real y, Real z, Real X, Real Y, Real Z, Real tol)
  {
    bool ok = true;
    ok &= rsIsCloseTo(x, X, tol);
    ok &= rsIsCloseTo(y, Y, tol);
    ok &= rsIsCloseTo(z, Z, tol);
    return ok;
  };


  // Test HSL/RGB conversions:
  for(int i = 1; i < N; i++) {
    x = Real(i) / Real(N);
    for(int j = 1; j < N; j++) {
      y = Real(j) / Real(N);
      for(int k = 1; k < N; k++) {
        z = Real(k) / Real(N);
        Color::hsl2rgb(x, y, z, &a, &b, &c);
        Color::rgb2hsl(a, b, c, &X, &Y, &Z);
        ok &= check(x,y,z, X,Y,Z, tol); }}}
  // Seems like the roundtrip HSL -> RGB -> HSL does not work when L = 0. It maps to black and we 
  // loose the hue and saturation information. The same thing happens when L = 1: it maps to white
  // and we also can't recover any hue information from that. That's why the loops start at 1 and
  // run only up to N-1, such that we avoid these extreme cases.
  // ToDo: add tests for the edge cases as well

  // Test CIELab to CIEXYZ:
  x = 0.3234436f;
  y = 0.734534f;
  z = 0.413245f;

  Color::lab2xyz(x, y, z, &a, &b, &c);
  Color::xyz2lab(a, b, c, &X, &Y, &Z);
  ok &= check(x,y,z, X,Y,Z, tol);

  Color::xyz2lab(x, y, z, &a, &b, &c);
  Color::lab2xyz(a, b, c, &X, &Y, &Z);
  ok &= check(x,y,z, X,Y,Z, tol);


  x = 0.0f;
  y = 0.0f;
  z = 1.0f;
  float L;
  Color::xyz2lab(x, y, z, &L, &a, &b);
  // (1,0,0) -> (  0,  431,    0)
  // (0,1,0) -> (100, -431,  127)
  // (0,0,1) -> (  0,    0, -172)

  L =  70;
  a =  50;
  b = -30;
  Color::lab2xyz(L, a, b, &x, &y, &z);
  // compare to https://www.nixsensor.com/free-color-converter/
  // it says that L is in 0..100 and a,b in -128...+128


  // Test conversion of RGB to hex colors:
  char hex[8];
  using uchar = unsigned char;
  using Str   = std::string;
  Color::rgb2hex((uchar)  0, (uchar)  0, (uchar)  0, hex); ok &= Str(hex) == Str("#000000");
  Color::rgb2hex((uchar)255, (uchar)255, (uchar)255, hex); ok &= Str(hex) == Str("#FFFFFF");
  Color::rgb2hex((uchar)254, (uchar)254, (uchar)254, hex); ok &= Str(hex) == Str("#FEFEFE");
  Color::rgb2hex((uchar)125, (uchar) 42, (uchar)143, hex); ok &= Str(hex) == Str("#7D2A8F");
  Color::rgb2hex((uchar)106, (uchar) 82, (uchar)160, hex); ok &= Str(hex) == Str("#6A52A0");
  // See: https://www.rapidtables.com/convert/color/rgb-to-hex.html



  // https://en.wikipedia.org/wiki/HSL_and_HSV#Hue_and_chroma
  // https://stackoverflow.com/questions/11980292/how-to-wrap-around-a-range

  return ok;
}
