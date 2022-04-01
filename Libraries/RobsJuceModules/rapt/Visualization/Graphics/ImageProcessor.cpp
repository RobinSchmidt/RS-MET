template<class T>
void rsImageProcessor<T>::gammaCorrection(rsImage<T>& img, T gamma)
{
  T* p = img.getPixelPointer(0, 0);
  for(int i = 0; i < img.getNumPixels(); i++)
    p[i] = pow(p[i], gamma);
}

template<class T>
void rsImageProcessor<T>::invert(rsImage<T>& img)
{
  T* p = img.getPixelPointer(0, 0);
  for(int i = 0; i < img.getNumPixels(); i++)
    p[i] = T(1) - p[i];
}

template<class T>
void rsImageProcessor<T>::normalize(rsImage<T>& img)
{
  T* p = img.getPixelPointer(0, 0);
  int N = img.getNumPixels();
  T min = rsArrayTools::minValue(p, N);
  for(int i = 0; i < N; i++)
    p[i] -= min;
  T max = rsArrayTools::maxValue(p, N);
  for(int i = 0; i < N; i++)
    p[i] /= max;
}
// this *may* be better numerically (less prone to roundoff errors) than the "fast" version -
// needs test

template<class T>
void rsImageProcessor<T>::normalizeFast(rsImage<T>& img)
{
  T* p = img.getPixelPointer(0, 0);
  int N = img.getNumPixels();

  T min = rsArrayTools::minValue(p, N);
  T max = rsArrayTools::maxValue(p, N);
  // todo: write function that extracts min and max in a single pass (may be more efficient due to
  // less data transfer from memory to processor)

  T scl = 1.f / (max-min);
  for(int i = 0; i < N; i++)
    p[i] = scl * (p[i] - min);
}

template<class T>
void rsImageProcessor<T>::normalizeJointly(rsImage<T>& img1, rsImage<T>& img2)
{
  //rsAssert(img2.hasSameShapeAs(img1));  // activate later
  using AT = rsArrayTools;
  int N = img1.getNumPixels();
  T* p1 = img1.getPixelPointer(0, 0);
  T* p2 = img2.getPixelPointer(0, 0);
  T min = rsMin(AT::minValue(p1, N), AT::minValue(p2, N));
  for(int i = 0; i < N; i++) {
    p1[i] -= min;
    p2[i] -= min; }
  T max = rsMax(AT::maxValue(p1, N), AT::maxValue(p2, N));
  for(int i = 0; i < N; i++) {
    p1[i] /= max;
    p2[i] /= max; }
}

template<class T>
rsImage<T> rsImageProcessor<T>::scaleUp(const rsImage<T>& img, int scl)
{
  int w = img.getWidth();
  int h = img.getHeight();
  rsImage<T> result(scl*w, scl*h);
  for(int x = 0; x < w; x++)  {
    for(int y = 0; y < h; y++) {
      for(int i = 0; i < scl; i++) {
        for(int j = 0; j < scl; j++) {
          result(scl*x+i, scl*y+j) = img(x, y); }}}}
  return result;
}
// todo:
// -allow different scaling factors for x and y
// -let the outer loop run over y and the inner over x

template<class T>
void rsImageProcessor<T>::sineShape(rsImage<T>& img)
{
  T* p = img.getPixelPointer(0, 0);
  for(int i = 0; i < img.getNumPixels(); i++)
    p[i] = T(0.5)*(sin(T(PI)*(p[i]-T(0.5)))+T(1));
}

//=================================================================================================
// rsImageContourPlotter

template<class TPix, class TVal>
rsImageContourPlotter<TPix, TVal>::rsImageContourPlotter()
{
  painter.setDeTwist(true);
}

template<class TPix, class TVal>
rsImage<TPix> rsImageContourPlotter<TPix, TVal>::getContourLines(const rsImage<TPix>& z,
  const std::vector<TVal>& levels,  const std::vector<TPix>& colors, bool antiAlias)
{
  rsImage<TPix> img(z.getWidth(), z.getHeight());
  for(size_t i = 0; i < levels.size(); i++)
    drawContour(z, levels[i], img, colors[i % colors.size()], antiAlias);
  return img;
}

template<class TPix, class TVal>
rsImage<TPix> rsImageContourPlotter<TPix, TVal>::getContourFills(const rsImage<TPix>& z,
  const std::vector<TVal>& levels, const std::vector<TPix>& colors, bool antiAlias)
{
  rsImage<TPix> img(z.getWidth(), z.getHeight());
  size_t j = 0; // color index
  size_t nc = colors.size();
  for(size_t i = 0; i < levels.size()-1; i++) {
    fillBetweenContours(z, levels[i], levels[i+1], img, colors[j % nc], antiAlias);
    j++; }
  return img;
}

// internal subroutines:

template<class TPix, class TVal>
void rsImageContourPlotter<TPix, TVal>::drawContour(
  const rsImage<TVal>& z, TVal level, rsImage<TPix>& target, TPix color, bool antiAlias)
{
  painter.setImageToPaintOn(&target);
  for(int i = 0; i < z.getWidth()-1; i++) {
    for(int j = 0; j < z.getWidth()-1; j++) {
      TVal z00 = z(i,   j  );
      TVal z01 = z(i,   j+1);
      TVal z10 = z(i+1, j  );
      TVal z11 = z(i+1, j+1);
      TVal min = rsMin(z00, z01, z10, z11);
      TVal max = rsMax(z00, z01, z10, z11);
      if(min < level && max >= level) {
        TVal x(0), y(0), w(1); // x,y offsets and weight for color
        if(antiAlias)
          contourSubPixelPosition(z00, z01, z10, z11, level, &x, &y, &w);
        painter.paintDot(TVal(i) + x, TVal(j) + y, w * color); }}}
}
// if we do not anti-alias, we need not to call the expensive paintDot and can use the cheaper
// painter.plot instead ...i think
// maybe don't loop over all pixels and follow the contours instead - but then there's no guarantee that
// nothing is missed

template<class TPix, class TVal>
void rsImageContourPlotter<TPix, TVal>::fillBetweenContours(const rsImage<TVal>& z,
  TVal lo, TVal hi, rsImage<TPix>& target, TPix fillColor, bool antiAlias)
{
  painter.setImageToPaintOn(&target);
  for(int i = 0; i < z.getWidth()-1; i++) {
    for(int j = 0; j < z.getWidth()-1; j++) {
      TVal z00 = z(i,   j  );
      TVal z01 = z(i,   j+1);
      TVal z10 = z(i+1, j  );
      TVal z11 = z(i+1, j+1);
      TVal min = rsMin(z00, z01, z10, z11);
      TVal max = rsMax(z00, z01, z10, z11);
      if(min >= lo && max < hi)
        painter.plot(i, j, fillColor);
      else  {
        // We are either on the contour or totally outside the drawing area
        if(!antiAlias) 
        {
          // What are we doing here? I think, the same as in the branch below but using a fixed
          // interpolation weight of 0.5 for blending between 0 (black) and the desired fillColor
          // instead of using weight baed on pixel coverage.
          /*
          if(min < lo && max >= lo)                      // on low contour
            painter.plot(i, j, fillColor * TPix(0.5));
          else if(min < hi && max >= hi)                 // on hi contour
            painter.plot(i, j, fillColor * TPix(0.5)); 
            */

          // Test:
          if(min < lo && max >= lo)                      // on low contour
            painter.plot(i, j, fillColor);
          else if(min < hi && max >= hi)                 // on hi contour
            painter.plot(i, j, TPix(0.0)); 

          /*
          if(min < lo && max >= lo)                      // on low contour
            painter.plot(i, j,TPix(0.0));
          else if(min < hi && max >= hi)                 // on hi contour
            painter.plot(i, j, fillColor); 
            */
          // OK - yes - these two tested variants may also make sense - check with the contours()
          // experiment. Maybe the bool antiAlias should be an int allowing to switch between these
          // 4 modes: coverage, average, high, low - these are the values assigned to the contours,
          // i.e. the 1-pixel wide boundaries between the contour fill colors
        }
        else {
          TVal c; // coverage
          if(min < lo && max >= lo) {
            c = contourPixelCoverage(z00, z01, z10, z11, lo);
            painter.plot(i, j, TPix(TVal(1)-c)*fillColor); }
          else if(min < hi && max >= hi) {
            c = contourPixelCoverage(z00, z01, z10, z11, hi);
            painter.plot(i, j, TPix(c)*fillColor);  }}}}}
}
// if instead of using:
//   if(min >= lo && max <  hi)
// we would use
//   if(min >  lo && max <  hi) -> leaves extra pixels blank (test with circle)
//   if(min >= lo && max <= hi) -> colors extra pixels in
//   if(min >  lo && max <= hi) -> no etra blank or colored pixels but ugly jaggies
// so the chosen variant seems best. this can be tested using the circles (and maybe commenting
// out the code that handles the contour lines - i think it was set to somewhere around 11 or 12
// levels...not sure anymore)

// can we refactor these two functions to avoid the duplicaztion?

template<class TPix, class TVal>
void rsImageContourPlotter<TPix, TVal>::contourSubPixelPosition(
  TVal z00, TVal z01, TVal z10, TVal z11, TVal c, TVal* x, TVal* y, TVal* weight)
{
  // Get line equation coeffs and evaluate the equation at midpoint to get the center of the
  // contour segment. The weight is given by the length divided by sqrt(2) such that diagonals get
  // weight 1.0:
  TVal x0, x1, y0, y1;
  contourSegmentCoeffs(z00, z01, z10, z11, c, x0, y0, x1, y1);
  TVal dx = (x1-x0);
  TVal dy = (y1-y0);
  *x = x0 + TVal(0.5) * dx;
  *y = y0 + TVal(0.5) * dy;
  *weight = sqrt(dx*dx + dy*dy) / sqrt(TVal(2));  // full weight only for diagonals
  //*weight = max(dx, dy);  // alternative - looks worse: screw-effect stronger and there are holes
}
// simpler idea:
// -compute z0 = (z00 + z01) / 2, z1 = (z10 + z11) / 2
//  z0 is the average value on the left, z1 on the right
// -the x-value/offset is determined by how much this average is above/below the target
//  level
// -similar for y
// -or is it the other way around?
// -might be even better than the center of the line

template<class TPix, class TVal>
TVal rsImageContourPlotter<TPix, TVal>::contourPixelCoverage(
  TVal z00, TVal z01, TVal z10, TVal z11, TVal c)
{
  TVal x0, x1, y0, y1;
  TVal A(0);     // covered area
  TVal h(0.5);   // half
  TVal I(1);     // one
  int branch = contourSegmentCoeffs(z00, z01, z10, z11, c, x0, y0, x1, y1);
  switch(branch) {
  case 0: { A = h *    x1  * y0;     if(z00 >= c) A = I-A; } break; // top-left
  case 1: { A = h * (I-x0) * y1;     if(z10 >= c) A = I-A; } break; // top-right
  case 2: { A = h *    x1  * (I-y0); if(z01 >= c) A = I-A; } break; // bottom-left
  case 3: { A = h * (I-x0) * (I-y1); if(z11 >= c) A = I-A; } break; // bottom-right
  case 4: {                              // horizontalish
    if(y0 < y1)  A = y1 + h * (y0-y1);   //   sloping up
    else         A = y0 + h * (y1-y0);   //   sloping down
    if(z00 >= c || z10 >= c)
      A = I-A; } break;
  case 5: {                              // verticalish
    if(x0 < x1)  A = x1 + h * (x0-x1);   //   leaning left
    else         A = x0 + h * (x1-x0);   //   leaning right
    if(z00 >= c || z01 >= c)
      A = I-A; } break; }
  return A;
}
// These simplified formulas (compared to the general formula for traingle areas) work only
// because we know in which order contourSegmentCoeffs returns the coeffs. Maybe we should make it
// swappable whether to use >= or < - sometimes we may want to invert the result - when drawing the
// bin-fills, we sometimes want to fill with the inverted weight ..i think - figure out - if so,
// maybe use a boolean and or let the user pass a comparison function cmp(z00, c), etc... or call
// it like inside(z00, c) or outside(z00, c)

template<class TPix, class TVal>
int rsImageContourPlotter<TPix, TVal>::contourSegmentCoeffs(
  TVal z00, TVal z01, TVal z10, TVal z11, TVal c, TVal& x0, TVal& y0, TVal& x1, TVal& y1)
{
  int branch;
  if((z00 < c && z01 >= c) || (z00 >= c && z01 < c)) {        // segment goes through left border
    x0 = TVal(0);
    y0 = rsLinToLin(c, z00, z01, TVal(0), TVal(1));
    if((z00 < c && z10 >= c) || (z00 >= c && z10 < c)) {      // segment goes through top border
      branch = 0;                                             //   -> top-left
      x1 = rsLinToLin(c, z00, z10, TVal(0), TVal(1));
      y1 = TVal(0); }
    else if((z01 < c && z11 >= c) || (z01 >= c && z11 < c)) { // segment goes through bottom border
      branch = 2;                                             //   -> bottom-left
      x1 = rsLinToLin(c, z01, z11, TVal(0), TVal(1));
      y1 = TVal(1); }
    else {                                                    // segment goes through right border
      branch = 4;                                             //   -> horizontalish
      x1 = TVal(1);
      y1 = rsLinToLin(c, z10, z11, TVal(0), TVal(1)); }}
  else {                                                      // doesn't go through left border
    if((z00 < c && z10 >= c) || (z00 >= c && z10 < c)) {      // goes through top border
      x0 = rsLinToLin(c, z00, z10, TVal(0), TVal(1));
      y0 = TVal(0);
      if((z10 < c && z11 >= c) || (z10 >= c && z11 < c)) {    // goes through right border
        branch = 1;                                           //   -> top-right
        x1 = TVal(1);
        y1 = rsLinToLin(c, z10, z11, TVal(0), TVal(1)); }
      else  {                                                 // goes through bottom border
        branch = 5;                                           //   -> verticalish
        x1 = rsLinToLin(c, z01, z11, TVal(0), TVal(1));
        y1 = TVal(1); }}
    else  {                                                   // doesn't go through top border
      branch = 3;                                             //   -> bottom-right
      x0 = rsLinToLin(c, z01, z11, TVal(0), TVal(1));
      y0 = TVal(1);
      x1 = TVal(1);
      y1 = rsLinToLin(c, z10, z11, TVal(0), TVal(1)); }}
  return branch;
}
// optimize the calls to rsLinToLin to get rid of divisions where possible - keep this code as
// prototype for unit testing the optimized code - i think, it's not possible, but we may get rid
// of some of the multiplications because outMax-outMin = 1 - make a function rsLinTo01, have a
// similar rs01ToLin
// note that the order of (x0,y0),(x1,y1) can't be changed without breaking contourPixelCoverage
// maybe the logical statements can be simplified by checking things like
// if (z00-c)*(z01-c) < 0,  >= 0 instead of the complicated and-or statements - but keep this
// version for unit tests


/*

Ideas:
-apply filtering to the hue channel - before that, we need to unwrap the hue. that's similar to 
 phase unwrapping in audio, but instead of just using the (already unwrapped) left neighbor as 
 reference, we need to consider left, top and maybe topleft neighbors. Maybe we should take the one
 as reference, which is closest in hue to the current. or maybe also take the other channels
 (luminance, saturation) into account. or maybe take an average. top row and right column are 
 treated specially (using only left and top neighbor respectively)

-Compute mean and variance per pixel:
 -mean is just obtained by gaussian blur
 -variance: subtract mean from original -> square it -> blur it -> take sqrt

-HDR-like effect (i hope):
 -assume the input is a floating-point rendering, i.e. dark pixels have nevertheless an 
  uncompromised relative precision due to the float format, so we don't really need the multiple
  different brightness versions of the same image that arise in HDR photography
 -adjust (per pixel) brightness and contrast according to some function that is based on the 
  (per pixel) values of mean and variance, i.e. apply a nonlinear mapping to each pixel's 
  brightness that is controlled by the mean and variance
 -maybe a similar process should be applied to saturation, too (but not to hue)


Should perhaps go into some rsImageAnalyzer class:

-Implement blob-coloring:
 https://en.wikipedia.org/wiki/Blob_detection
 https://www.youtube.com/watch?v=vTUsGzXFmuc
-Define a data-structure for representing image regions (can be used for "blobs"). A region is 
 defined to be a set of pixels, typically (but not necessarily) connected. It can be implemented 
 just as a vector of pairs of integer pixel coordinates x,y.
-On such regions, define functions to extract certain features (see "Taschenbuch der Informatik", 
 pg. 540):
 -Area: number of pixels
 -Circumfence: number of boundary pixels, defined as pixels who have neighbors that are not part of 
  the region (should we use a 4- or 8-neighborhood? Maybe have both variants). Computing it perhaps
  needs an O(N^2) algo where N is the number of pixels? Or can we do better? Perhaps, if the pixels 
  are sorted in some way, we can indeed do better: left or right neighbors, if present, should be 
  immediately adjacent in the array. top or bottom neighbors, if present, should be at most one 
  "row-stride" away (there is no fixe row-stride though, but we may be able to give bounds). Maybe
  it could be helpful, if the region data-structure stores row-indices and within each row, the 
  indices of columns - that could even be more memory efficient (at least for compact regions)
 -Compactness: circumfence^2 / area
 -Connectedness: each pixel has at least one neighbor which is also inside the region
 -Moments: m_{kl} = sum_{x,y} x^k  y^l  b(x,y) where x,y are the pixel coordinates and b(x,y) is
  the pixel's brightness
 -Center of gravity: x_c = m_{10} / m_{00}, y_c = m_{01} / m_{00}
 -Traslation invariant moments: t_{kl} = sum_{x,y} (x-x_c)^k (y-y_c)^l b(x,y). They don't depend on
  the position of the region within the image
 -Translation- and scale-invariant moments: s_{kl} = t_{kl} / pow(t_{00}, 1 + (k+l)/2)
 -Translation-, scale- and rotation-invariant moments: r_1 = t_{20} - t_{02}, 
  r_2 = r_1^2 + 4 t_{11}^2, r_3 = (t_{30} - 3 t_{12})^2 + (3 t_{21} - t_{03})^2
 -Histogram: function from the pixel-values (brightnesses) to a relative frequency of occurence in
  0..1. When pixel values are float, we may need a function from brightness intervals (instead of 
  brightness values) to relative frequencies. From the histogram, we may extract further features
  such as mean brightness, brightness variance...maybe measures like skew, kurtosis, etc. could be
  interesting...maybe some sort of measure of uni-/bi-/tri-/etc-modality? ...although, presumably,
  any region extraction algo will have selected more or less homogenuous regions anyway, so these
  features may mae more sense for an image as a whole?

-Try to classify pixels as background, foreground and boundary (as in anti-aliased vector graphics
 rendering) and try to change the background color without altering the foreground: Example: green
 filled triangle on red background -> switch background to blue -> boundary pixles should change 
 from dark yellowish to purple. Maybe there can be more pixel classes than these 3.

-For blurring, see this: https://www.youtube.com/watch?v=LKnqECcg6Gw
 We should square the pixel brighnesses, do the blur, then take the sqrt. Maybe we can let the user
 specify an exponent or maybe even a pair of functions for pre- and post processing. On the other
 hand, the caller could actually do pre- and post-procesing themselves. Maybe for common exponents
 such as 1 or 2, we can use a cheaper implementation. 


*/