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
  rsImageF result(scl*w, scl*h);
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
rsImage<TPix> rsImageContourPlotter<TPix, TVal>::getContourLines(const rsImage<TPix>& z, 
  const std::vector<TVal>& levels,  const std::vector<TPix>& colors, bool antiAlias)
{
  rsImageF c(z.getWidth(), z.getHeight());
  for(size_t i = 0; i < levels.size(); i++)
    drawContour(z, levels[i], c, colors[i % colors.size()], antiAlias);
  return c;
}

template<class TPix, class TVal>
rsImage<TPix> rsImageContourPlotter<TPix, TVal>::getContourFills(const rsImage<TPix>& z, 
  const std::vector<TVal>& levels, const std::vector<TPix>& colors, bool antiAlias)
{
  rsImageF imgBins(z.getWidth(), z.getHeight());  // fills
  size_t j = 0; // color index
  size_t nc = colors.size();
  for(size_t i = 0; i < levels.size()-1; i++) {
    fillBetweenContours(z, levels[i], levels[i+1], imgBins, colors[j % nc], antiAlias);
    j++; }
  return imgBins;
}


// internal subroutines:

template<class TPix, class TVal>
void rsImageContourPlotter<TPix, TVal>::drawContour(
  const rsImage<TVal>& z, TVal level, rsImage<TPix>& target, TPix color, bool antiAlias)
{
  rsImagePainter<TPix, TVal, TVal> painter(&target);
  painter.setDeTwist(true);  
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
  rsImagePainter<TPix, TVal, TVal> painter(&target);
  for(int i = 0; i < z.getWidth()-1; i++) {
    for(int j = 0; j < z.getWidth()-1; j++) {
      TVal z00 = z(i, j);
      TVal z01 = z(i, j+1);
      TVal z10 = z(i+1, j);
      TVal z11 = z(i+1, j+1);
      TVal min = rsMin(z00, z01, z10, z11);
      TVal max = rsMax(z00, z01, z10, z11);
      if(min >= lo && max < hi)     // this seems to be artifact-free
        painter.plot(i, j, fillColor);
      else  {
        // we are either on the contour or totally outside the drawing area
        if(!antiAlias)  {                  // ok - this looks right
          if(min < lo && max >= lo)                      // on low contour
            painter.plot(i, j, fillColor * TPix(0.5));
          else if(min < hi && max >= hi)                 // on hi contour
            painter.plot(i, j, fillColor * TPix(0.5)); }
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
//   if(min >= lo && max < hi) 
// we would use
//   if(min > lo && max < hi)   -> leaves extra pixels blank (test with circle)
//   if(min >= lo && max <= hi) -> colors extra pixels in
//   if(min > lo && max <= hi)  -> no etra blank or colored pixels but ugly jaggies
// so the chosen variant seems best. this can be tested using the circles (and maybe commenting
// out the code that handles the contour lines - i think it was set to somewhere around 11 or 12 
// levels...not sure anymore)



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
  TVal A = 0.f;  // covered area
  TVal h(0.5);   // half
  TVal I(1.0);   // one
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
// These simplified formulas (compared to the general formula for traingel areas) work only 
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
    x0 = 0.f;
    y0 = rsLinToLin(c, z00, z01, 0.f, 1.f);
    if((z00 < c && z10 >= c) || (z00 >= c && z10 < c)) {      // segment goes through top border
      branch = 0;                                             //   -> top-left
      x1 = rsLinToLin(c, z00, z10, 0.f, 1.f);
      y1 = 0.f; }
    else if((z01 < c && z11 >= c) || (z01 >= c && z11 < c)) { // segment goes through bottom border
      branch = 2;                                             //   -> bottom-left
      x1 = rsLinToLin(c, z01, z11, 0.f, 1.f);
      y1 = 1.f; }
    else {                                                    // segment goes through right border
      branch = 4;                                             //   -> horizontalish
      x1 = 1.f;                                               
      y1 = rsLinToLin(c, z10, z11, 0.f, 1.f); }}
  else {                                                      // doesn't go through left border
    if((z00 < c && z10 >= c) || (z00 >= c && z10 < c)) {      // goes through top border
      x0 = rsLinToLin(c, z00, z10, 0.f, 1.f);
      y0 = 0.f;
      if((z10 < c && z11 >= c) || (z10 >= c && z11 < c)) {    // goes through right border
        branch = 1;                                           //   -> top-right
        x1 = 1.f;
        y1 = rsLinToLin(c, z10, z11, 0.f, 1.f); }
      else  {                                                 // goes through bottom border
        branch = 5;                                           //   -> verticalish
        x1 = rsLinToLin(c, z01, z11, 0.f, 1.f);            
        y1 = 1.f; }}
    else  {                                                   // doesn't go through top border 
      branch = 3;                                             //   -> bottom-right
      x0 = rsLinToLin(c, z01, z11, 0.f, 1.f);
      y0 = 1.f;
      x1 = 1.f;
      y1 = rsLinToLin(c, z10, z11, 0.f, 1.f); }}
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
