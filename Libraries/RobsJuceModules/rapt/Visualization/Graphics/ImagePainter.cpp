template<class TPix, class TWgt, class TCor>
rsImagePainter<TPix, TWgt, TCor>::rsImagePainter(rsImage<TPix> *imageToPaintOn, rsAlphaMask<TWgt> *maskToUse)
{
  //antiAlias = true;
  //useMask = false;
  setNeighbourWeightsForSimpleDot(0, 0);
  setImageToPaintOn(imageToPaintOn);
  setAlphaMaskForDot(maskToUse);
}

// setup

template<class TPix, class TWgt, class TCor>
void rsImagePainter<TPix, TWgt, TCor>::setImageToPaintOn(rsImage<TPix> *imageToPaintOn)
{
  image = imageToPaintOn;
  //if(image != nullptr)
  //{
  //  wi = image->getWidth();
  //  hi = image->getHeight();
  //}
  //else
  //  wi = hi = 0;
}

template<class TPix, class TWgt, class TCor>
void rsImagePainter<TPix, TWgt, TCor>::setAlphaMaskForDot(rsAlphaMask<TWgt> *maskToUse)
{
  mask = maskToUse;
  //wb = mask->getWidth();
  //hb = mask->getHeight();
  // todo: handle nullptr case
}

template<class TPix, class TWgt, class TCor>
void rsImagePainter<TPix, TWgt, TCor>::setNeighbourWeightsForSimpleDot(TWgt straight, TWgt diagonal)
{
  straightNeighbourWeight = straight;
  diagonalNeighbourWeight = diagonal;
}

/*
template<class TPix, class TWgt, class TCor>
void rsImagePainter<TPix, TWgt, TCor>::setAntiAlias(bool shouldAntiAlias)
{
  antiAlias = shouldAntiAlias;
}

template<class TPix, class TWgt, class TCor>
void rsImagePainter<TPix, TWgt, TCor>::setUseAlphaMask(bool shouldUseMask)
{
  useMask = shouldUseMask;
}
*/

// painting

template<class TPix, class TWgt, class TCor>
void rsImagePainter<TPix, TWgt, TCor>::paintDot(TCor x, TCor y, TPix color)
{
  //rsAssert(rsIsFiniteNumber(x));
  //rsAssert(rsIsFiniteNumber(y));
  //rsAssert(rsIsFiniteNumber(color));

  // todo: get rid of this dispatcher code - use a member function pointer instead.
  if(useMask)
  {
    if(antiAlias)
      paintDotViaMask(x, y, color);
    else
      paintDotViaMask((int) ::round(x), (int) ::round(y), color);
  }
  else
  {
    if(antiAlias)
      paintDot3x3(x, y, color, straightNeighbourWeight, diagonalNeighbourWeight);
    else
      paintDot3x3((int) ::round(x), (int) ::round(y), color, 
        straightNeighbourWeight, diagonalNeighbourWeight);
  }
}

template<class TPix, class TWgt, class TCor>
void rsImagePainter<TPix, TWgt, TCor>::paintDot3x3(int x, int y, TPix color, TWgt weightStraight, 
  TWgt weightDiagonal)
{
  int wi = image->getWidth();
  int hi = image->getHeight();

  if(x >= 0 && x < wi && y >= 0 && y < hi)
    plot(x, y, color);

  // apply thickness:
  if(weightStraight > 0.f && x >= 1 && x < wi-1 && y >= 1 && y < hi-1)
  {
    TPix a, ta, sa;
    a  = color;
    ta = a * (TPix)weightStraight;
    sa = a * (TPix)weightDiagonal;

    plot(x-1, y-1, sa);
    plot(x,   y-1, ta);
    plot(x+1, y-1, sa);

    plot(x-1, y,   ta);
    plot(x+1, y,   ta);

    plot(x-1, y+1, sa);
    plot(x,   y+1, ta);
    plot(x+1, y+1, sa);
  }
}

template<class TPix, class TWgt, class TCor>
void rsImagePainter<TPix, TWgt, TCor>::paintDot3x3(TCor x, TCor y, TPix color, TWgt weightStraight, 
  TWgt weightDiagonal)
{
  int wi = image->getWidth();
  int hi = image->getHeight();

  int xi = (int)floor(x);  // integer part of x
  int yi = (int)floor(y);  // integer part of y
  x -= xi;                 // fractional part of x
  y -= yi;                 // fractional part of y

  // compute weights for bilinear deinterpolation/spreading (maybe factor out):
  TPix a, b, c, d;
  d = TPix(x*y);
  c = TPix(y)-d;
  b = TPix(x)-d;
  a = d+TPix(1-x-y);
  // formulas have been obtained by the condition that a+b+c+d = 1...right? but maybe their square 
  // should sum to unity? areas with larger spreading appear darker

  //TPix sum = a+b+c+d;  // test - should be unity

  // optionally normalize make the sum of squares of a,b,c,d to a constant (currently, sum of the 
  // values a,b,c,d themselves is unity) - the factor 0.5 makes it visually similar to no 
  // normalization (maybe try sqrt(0.5) - it seems a bit darker compared to non-anti-aliased mode
  // - maybe compute the ratio between sum-of-pixel values in non-aliased and anti-aliased drawing 
  // forsome "typical" curve and use that factor):
  if(deTwist) {
    TPix s = TPix(0.5) / sqrt(a*a + b*b + c*c + d*d);
    //TPix s = TPix(0.5 + rsSquare(0.5-x) + rsSquare(0.5-y)) / sqrt(a*a + b*b + c*c + d*d); // test - not good!
    a *= s; b *= s; c *= s; d *= s; }
  // maybe try different formulas that normalize in a ways such that a^2+b^2 = y, c^2+d^2 = 1-y,
  // a^2+c^2 = 1-x, b^2+d^2 = x ...or something - but maybe that's not possible

  //sum = a+b+c+d;  // 1, when (x,y) = (0.5,0.5) and 0.5 when (x,y) = (0,0)
  // can we find a formula that returns 1 in "both" cases? what about 
  // sqrt(2) / sqrt(2*(a^2+b^2+c^2+d^2)) - no - what about suing x*y and/or x+y in the numerator in 
  // order to increase the value when (x,y) = (0,0) - we wnat the numerator to behave like:
  // n(0.5,0.5) = 0.5, n(0,0) = 2, ..n(x,y) = 0.5 + 0.5*x*x + 0.5*y*y or
  // n(x,y) = 0.5 + (0.5-x)^2 + (0.5-y)^2 - yes - this formula seems to have the desired behavior, 
  // so try s = TPix(0.5 + (0.5-x)^2 + (0.5-y)^2) / sqrt(a*a + b*b + c*c + d*d);
  // ..but what is the behavior in between? maybe this formula defeats the purpose of de-twisting 
  // in the sense that the numerator undoes the desired nomrlaization -> experiment! maybe plot
  // s(x,y) in the unit-square
  // ...ok - tried it the new formula brings back the twisties

  // compute final coeffs and accumulate values into the 4 pixels:
  a *= color; b *= color; c *= color; d *= color;
  if(xi >= 0 && xi < wi-1 && yi >= 0 && yi < hi-1) {
    plot(xi,   yi,   a);
    plot(xi+1, yi,   b);
    plot(xi,   yi+1, c);
    plot(xi+1, yi+1, d); }

  // apply thickness:
  if(weightStraight > 0.f && xi >= 1 && xi < wi-2 && yi >= 1 && yi < hi-2)
  {
    TPix t, s, sa, sb, sc, sd, ta, tb, tc, td;
    t = (TPix)weightStraight;      // weight for direct neighbour pixels
    s = (TPix)weightDiagonal;      // weight for diagonal neighbour pixels

    sa = s*a;
    sb = s*b;
    sc = s*c;
    sd = s*d;

    ta = t*a;
    tb = t*b;
    tc = t*c;
    td = t*d;

    plot(xi-1, yi-1, sa);
    plot(xi,   yi-1, ta+sb);
    plot(xi+1, yi-1, tb+sa);
    plot(xi+2, yi-1, sb);

    plot(xi-1, yi,   ta+sc);
    plot(xi,   yi,   sd+tb+tc);
    plot(xi+1, yi,   sc+ta+td);
    plot(xi+2, yi,   tb+sd);

    plot(xi-1, yi+1, tc+sa);
    plot(xi,   yi+1, sb+ta+td);
    plot(xi+1, yi+1, sa+tb+tc);
    plot(xi+2, yi+1, td+sb);

    plot(xi-1, yi+2, sc);
    plot(xi,   yi+2, tc+sd);
    plot(xi+1, yi+2, td+sc);
    plot(xi+2, yi+2, sd);
  }
}

template<class TPix, class TWgt, class TCor>
void rsImagePainter<TPix, TWgt, TCor>::paintDotViaMask(int x, int y, TPix color)
{
  int wi = image->getWidth();
  int hi = image->getHeight();
  int wm = mask->getWidth();
  int hm = mask->getHeight();

  // write coordinates in target image:
  int xs = x  - wm/2;    // start x-coordinate
  y      = y  - hm/2;    // start y coordinate
  int xe = xs + wm-1;    // end x coordinate
  int ye = y  + hm-1;    // end y coordinate

  // read coordinates in mask image:
  int mxs = 0;          // start x
  int my  = 0;          // start y

  // checks to not write beyond image bounds:
  if(xs < 0)
  {
    mxs = -xs;
    xs  = 0;
  }
  if(y < 0)
  {
    my = -y;
    y  =  0;
  }
  if(xe >= wi)
    xe = wi-1;
  if(ye >= hi)
    ye = hi-1;

  // the actual painting loop over the pixels in target image and brush image:
  int mx;
  while(y <= ye)
  {
    x  = xs;
    mx = mxs;
    while(x <= xe)
    {
      //// debug:
      //rsAssert(x  >= 0 && x  < wi);
      //rsAssert(y  >= 0 && y  < hi);
      //rsAssert(mx >= 0 && mx < wm);
      //rsAssert(my >= 0 && my < hm);

      plot(x, y, color, (*mask)(mx, my));
      x++;
      mx++;
    }
    y++;
    my++;
  }
}

//template<class TPix, class TWgt, class TCor>
//void rsImagePainter<TPix, TWgt, TCor>::paintDotViaMask(TCor x, TCor y, TPix color)
//{
//  int wi = image->getWidth();
//  int hi = image->getHeight();
//  int wm = mask->getWidth();
//  int hm = mask->getHeight();
//
//  int  xi = (int)floor(x);  // integer part of x
//  int  yi = (int)floor(y);  // integer part of y
//  TCor xf = x - xi;         // fractional part of x
//  TCor yf = y - yi;         // fractional part of y
//
//                            // compute pixel coverages (these are the same as the weights for bilinear deinterpolation):
//  TWgt a, b, c, d;
//  d = TWgt(xf*yf);
//  c = TWgt(yf)-d;
//  b = TWgt(xf)-d;
//  a = d+TWgt(1-xf-yf);
//
//  // compute write coordinates in target image:
//  int xs = xi - wm/2;    // start x-coordinate in image
//  int ys = yi - hm/2;    // start y coordinate in image
//  int xe = xs + wm;      // end x coordinate in image
//  int ye = ys + hm;      // end y coordinate in image
//
//  if(xe < 0 || xs >= wi || ye < 0 || ys >= hi)
//    return; // dot completely off the canvas (...but do we need this check?)
//
//  // read coordinates in mask image:
//  int xms = 0;   // start x in mask
//  int yms = 0;   // start y in mask
//
//  // flags to indicate whether or not we need to draw the 4 edges of the mask:
//  bool leftEdge, rightEdge, topEdge, bottomEdge;
//  leftEdge   = true;
//  rightEdge  = true;
//  topEdge    = true;
//  bottomEdge = true;
//
//  // checks to not write beyond image bounds:
//  if(xs < 0)
//  {
//    xms = -xs-1;
//    xs  = 0;
//    leftEdge = false;
//  }
//  if(ys < 0)
//  {
//    yms = -ys-1;
//    ys  =  0;
//    topEdge = false;
//  }
//  if(xe >= wi)
//  {
//    xe = wi-1;
//    rightEdge = false;
//  }
//  if(ye >= hi)
//  {
//    ye = hi-1;
//    bottomEdge = false;
//  }
//
//  int xm, ym;    // x- and y-index in mask
//  TWgt w;        // weight for current pixel
//
//  //// paint edges and corners:
//  if(leftEdge)
//  {
//    ym = yms;
//    for(yi = ys+1; yi <= ye-1; yi++)
//    {
//      w = a * (*mask)(0, ym+1) + c * (*mask)(0, ym); 
//      plot(xs, yi, color, w);
//      ym++;
//    }
//    if(topEdge)
//      plot(xs, ys, color, a * (*mask)(0, 0));  // top left corner
//  }
//  if(topEdge)
//  {
//    xm = xms;
//    for(xi = xs+1; xi <= xe-1; xi++) 
//    {
//      w = a * (*mask)(xm+1, 0) + b * (*mask)(xm, 0); 
//      plot(xi, ys, color, w);
//      xm++;
//    }
//    if(rightEdge)
//      plot(xe, ys, color, b * (*mask)(wm-1, 0));  // top right corner
//  }
//  if(rightEdge)
//  {
//    ym = yms;
//    for(yi = ys+1; yi <= ye-1; yi++) 
//    {
//      w = b * (*mask)(wm-1, ym+1) + d * (*mask)(wm-1, ym); 
//      plot(xe, yi, color, w);
//      ym++;
//    }
//  }
//  if(bottomEdge)
//  {
//    xm = xms;
//    for(xi = xs+1; xi <= xe-1; xi++) 
//    {
//      w = c * (*mask)(xm+1, hm-1) + d * (*mask)(xm, hm-1); 
//      plot(xi, ye, color, w);
//      xm++;
//    }
//    if(leftEdge)
//      plot(xs, ye, color, c * (*mask)(0,    hm-1));  // bottom left corner
//    if(rightEdge)
//      plot(xe, ye, color, d * (*mask)(wm-1, hm-1));  // bottom right corner
//  }
//
//  // adjust start/endp pixels when edges were drawn:
//  if(leftEdge)   { xs++; xms++; }
//  if(topEdge)    { ys++; yms++; }
//  if(rightEdge)  { xe--;        }
//  if(bottomEdge) { ye--;        }
//
//  // paint interior rectangle:
//  ym = yms;
//  for(yi = ys; yi <= ye-1; yi++)
//  {
//    xm = xms;
//    for(xi = xs; xi <= xe-1; xi++)
//    {
//      w = a * (*mask)(xm+1, ym+1) + b * (*mask)(xm, ym+1)
//        + c * (*mask)(xm+1, ym)   + d * (*mask)(xm, ym);   // check this carefully!
//      plot(xi, yi, color, w); 
//      xm++;
//    }
//    ym++;
//  }
//}

template<class TPix, class TWgt, class TCor>
void rsImagePainter<TPix, TWgt, TCor>::paintDotViaMask(TCor x, TCor y, TPix color)
{
  // This almost works, but there's still something wrong with the edge code - there's a 1 pixel 
  // wide zone border in the image onto which nothing is drawn

  int wi = image->getWidth();
  int hi = image->getHeight();
  int wm = mask->getWidth();
  int hm = mask->getHeight();

  int  xi = (int)floor(x);  // integer part of x
  int  yi = (int)floor(y);  // integer part of y
  TCor xf = x - xi;         // fractional part of x
  TCor yf = y - yi;         // fractional part of y

  // compute pixel coverages (these are the same as the weights for bilinear deinterpolation):
  TWgt a, b, c, d;
  d = TWgt(xf*yf);
  c = TWgt(yf)-d;
  b = TWgt(xf)-d;
  a = d+TWgt(1-xf-yf);

  // compute write coordinates in target image:
  int xs = xi - wm/2;    // start x-coordinate in image
  int ys = yi - hm/2;    // start y coordinate in image
  int xe = xs + wm;      // end x coordinate in image
  int ye = ys + hm;      // end y coordinate in image

  if(xe < 0 || xs >= wi || ye < 0 || ys >= hi)
    return; // dot completely off the canvas (...but do we need this check?)

  // read coordinates in mask image:
  int xms = 0;   // start x in mask
  int yms = 0;   // start y in mask

  // flags to indicate whether or not we need to draw the 4 edges of the mask:
  bool leftEdge, rightEdge, topEdge, bottomEdge;
  leftEdge   = true;
  rightEdge  = true;
  topEdge    = true;
  bottomEdge = true;

  // checks to not write beyond image bounds:
  if(xs < 0)
  {
    xms = -xs;
    xs  = 0;
    leftEdge = false;
  }
  if(ys < 0)
  {
    yms = -ys;
    ys  =  0;
    topEdge = false;
  }
  if(xe >= wi)
  {
    xe = wi-1;
    rightEdge = false;
  }
  if(ye >= hi)
  {
    ye = hi-1;
    bottomEdge = false;
  }

  int xm, ym;    // x- and y-index in mask
  TWgt w;        // weight for current pixel

  // paint edges and corners:
  if(leftEdge)
  {
    ym = yms;
    for(yi = ys+1; yi <= ye-1; yi++)
    {
      w = a * (*mask)(0, ym+1) + c * (*mask)(0, ym); 
      plot(xs, yi, color, w);
      ym++;
    }
    if(topEdge)
      plot(xs, ys, color, a * (*mask)(0, 0));  // top left corner
  }
  if(topEdge)
  {
    xm = xms;
    for(xi = xs+1; xi <= xe-1; xi++) 
    {
      w = a * (*mask)(xm+1, 0) + b * (*mask)(xm, 0); 
      plot(xi, ys, color, w);
      xm++;
    }
    if(rightEdge)
      plot(xe, ys, color, b * (*mask)(wm-1, 0));  // top right corner
  }
  if(rightEdge)
  {
    ym = yms;
    for(yi = ys+1; yi <= ye-1; yi++) 
    {
      w = b * (*mask)(wm-1, ym+1) + d * (*mask)(wm-1, ym); 
      plot(xe, yi, color, w);
      ym++;
    }
  }
  if(bottomEdge)
  {
    xm = xms;
    for(xi = xs+1; xi <= xe-1; xi++) 
    {
      w = c * (*mask)(xm+1, hm-1) + d * (*mask)(xm, hm-1); 
      plot(xi, ye, color, w);
      xm++;
    }
    if(leftEdge)
      plot(xs, ye, color, c * (*mask)(0,    hm-1));  // bottom left corner
    if(rightEdge)
      plot(xe, ye, color, d * (*mask)(wm-1, hm-1));  // bottom right corner
  }

  // paint interior rectangle:
  ym = yms;
  for(yi = ys+1; yi <= ye-1; yi++)
  {
    xm = xms;
    for(xi = xs+1; xi <= xe-1; xi++)
    {
      w = a * (*mask)(xm+1, ym+1) + b * (*mask)(xm, ym+1)
        + c * (*mask)(xm+1, ym)   + d * (*mask)(xm, ym);   // check this carefully!
      plot(xi, yi, color, w); 
      xm++;
    }
    ym++;
  }
}

//template<class TPix, class TWgt, class TCor>
//void rsImagePainter<TPix, TWgt, TCor>::drawDottedLine(TCor x1, TCor y1, TCor x2, TCor y2, TPix color, 
//  TCor density, int maxNumDots, bool scaleByNumDots, TCor minDotDistance)
//{
//  // maybe we should factor out this function into the rsPhaseScopeBuffer class
//
//  TCor dx = x2-x1;
//  TCor dy = y2-y1;
//  TCor pixelDistance = sqrt(dx*dx + dy*dy);
//  int  numDots = rsMax(1, (int)floor(density*pixelDistance/minDotDistance));
//  if(maxNumDots > 0)
//    numDots = rsMin(numDots, maxNumDots);
//
//  TPix scaledColor = color;
//  if(scaleByNumDots)
//    scaledColor = scaledColor / (TPix)numDots;
//
//  TPix c1 = scaledColor;
//  TPix c2 = scaledColor;
//  drawLineDotted(x1, y1, x2, y2, c1, c2, numDots);
//    // preliminary - later we want to use a color gradient, i.e. pass two different colors
//}

template<class TPix, class TWgt, class TCor>
void rsImagePainter<TPix, TWgt, TCor>::drawLineDotted(TCor x1, TCor y1, TCor x2, TCor y2, 
  TPix c1, TPix c2, int numDots)
{
  TCor dx = x2-x1;
  TCor dy = y2-y1;
  TPix dc = c2-c1;

  TCor scaler = (TCor)(1.0 / numDots);
  TCor k;
  for(int i = 1; i <= numDots; i++)
  {
    k = scaler * i;  // == i / numDots
    paintDot(x1 + k*dx, y1 + k*dy, c1 + TPix(k)*dc);
  }

  // maybe we should rename this function to dottedLineTo (we don't draw the very 1st point of the
  // line because we assume it has already been drawn as endpoint of a previous line)

  // i think, we should start the loop at i=0 and use scaler = 1.0 / (numDots-1)
  //TCor scaler = (TCor)(1.0 / (numDots-1));
  //TCor k;
  //for(int i = 0; i < numDots; i++)
  //{
  //  k = scaler * i;  // == i / (numDots-1)
  //  paintDot(x1 + k*dx, y1 + k*dy, scaledColor);
  //}
}

template<class TPix, class TWgt, class TCor>
void rsImagePainter<TPix, TWgt, TCor>::drawDots(TCor* x, TCor* y, TPix* c, int numDots)
{
  for(int i = 0; i < numDots; i++)
    paintDot(x[i], y[i], c[i]);
}

/*
template<class TPix, class TWgt, class TCor>
void rsImagePainter<TPix, TWgt, TCor>::drawDottedSpline(TCor x1, TCor x1s, TCor y1, TCor y1s, 
  TCor x2, TCor x2s, TCor y2, TCor y2s, TPix c1, TPix c2, int numDots)
{
  TCor a[4], b[4];   // coeffs for x(t), y(t)
  cubicSplineArcCoeffs2D(x1, x1s, y1, y1s, x2, x2s, y2, y2s, a, b);

  bool desityCompensation = false;  // make parameter
  if(desityCompensation)
    drawDottedSpline2(a, b, c1, c2);
  else
    drawDottedSpline1(a, b, c1, c2, numDots);
}

template<class TPix, class TWgt, class TCor>
void rsImagePainter<TPix, TWgt, TCor>::drawDottedSpline1(TCor *a, TCor *b, TPix c1, TPix c2, 
  int numDots)
{
  TPix dc = c2-c1;  // color difference
  TCor scaler = (TCor)(1.0 / numDots); // not 1/(numDots-1) because last dot of this call is drawn
  TCor t, x, y;                        // as first dot in next call? ..avoids drawing it twice?
  for(int i = 1; i <= numDots; i++)
  {
    t = scaler * i;  // == i / numDots
    x = rsPolynomial<TCor>::evaluatePolynomialAt(t, a, 3);
    y = rsPolynomial<TCor>::evaluatePolynomialAt(t, b, 3);
    paintDot(x, y, c1 + TPix(t)*dc);
  }
}
// todo: fill t-array and call the function below - avoid code duplication
// obsolete soon

template<class TPix, class TWgt, class TCor>
void rsImagePainter<TPix, TWgt, TCor>::drawDottedSpline2(TCor *a, TCor *b, TPix c1, TPix c2, 
  TCor* t, int numDots)
{
  TPix dc = c2-c1;                      // color difference
  TCor scaler = (TCor)(1.0 / numDots);  // not 1/(numDots-1) because last dot of this call is drawn
  TCor ti, x, y;                        // as first dot in next call? ..avoids drawing it twice?
  for(int i = 0; i < numDots; i++) {
    ti = t[i];
    x = rsPolynomial<TCor>::evaluatePolynomialAt(ti, a, 3);
    y = rsPolynomial<TCor>::evaluatePolynomialAt(ti, b, 3);
    paintDot(x, y, c1 + TPix(ti)*dc); 
    // todo: make a function that just fills 3 arrays for x,y,c instead of directly painting, so it
    // can be used to delegate the actual painting to OpenGL
  }
}
// obsolete soon

template<class TPix, class TWgt, class TCor>
void rsImagePainter<TPix, TWgt, TCor>::drawDottedSpline2(TCor *a, TCor *b, TPix c1, TPix c2)
{
  // obtain arc-length s as (sampled) function of parameter t:
  static const int N = 17; // make this user-adjustable (setDensityCompensationPrecision)
  r.resize(N);
  s.resize(N);
  rsArrayTools::fillWithRangeLinear(&r[0], N, TCor(0), TCor(1));
  cubicSplineArcLength2D(a, b, &r[0], &s[0], N);

#ifdef RS_DEBUG_PLOTTING
  GNUPlotter::plot(N, &r[0], &s[0]);
#endif

  // Compute sequence of t-values that result in equidistant dots:
  TCor splineLength = s[N-1]; // last value in s is total length: s(t=1)
  TCor density = 0.125;       // preliminary - make parameter
  int numSplineDots = rsMax(1, rsRoundToInt(splineLength * density));
  u.resize(numSplineDots);
  t.resize(numSplineDots);
  TCor scaler = splineLength / TCor(numSplineDots-1.0);
  for(int i = 0; i < numSplineDots; i++)
    u[i] = i*scaler;
  resampleNonUniformLinear(&s[0], &r[0], N, &u[0], &t[0], numSplineDots);

#ifdef RS_DEBUG_PLOTTING
  GNUPlotter::plot(numSplineDots, &u[0], &t[0]);
#endif

  // Draw the spline with the computed t-values:
  drawDottedSpline2(a, b, c1, c2, &t[0], numSplineDots);

  // todo: try also a quartic spline that normalizes the integral under the functions x(t) and y(t)
  // to the same value that a linear interpolant would have. This should probably reduce overshoot.
  // Maybe also try to fix the 3rd derivatives at the joints to 0 using a quintic interpolating 
  // spline - figure out, how this affects the dot-density - maybe a quintic could can go without
  // desnity compensation (but probably not)

  // factor out a class that computes the dot-coordinates and colors (rsDottingSomething) because
  // that will be needed for the OpenGL port anyway
}
// obsolete soon
*/

// some helper functions used in Wu algorithm (maybe try to get rid of them):
template<class T> inline int        ipart(T x) { return (int) x;         }
template<class T> inline T          fpart(T x) { return x - ipart(x);    }
template<class T> inline T         rfpart(T x) { return 1 - fpart(x);    }
template<class T> inline int   roundToInt(T x) { return ipart(x + 0.5f); }
//template<class T> inline void swap(T& x, T& y) { T t = x; x = y; y = t;  }
//template<class T> inline float   min(T x, T y) { return x < y ? x : y;   }
 
// Wu line drawing algorithm translated from
// https://en.wikipedia.org/wiki/Xiaolin_Wu's_line_algorithm with a few obvious optimizations:
template<class TPix, class TWgt, class TCor>
void rsImagePainter<TPix, TWgt, TCor>::drawLineWu(TCor x0, TCor y0, TCor x1, TCor y1, TPix color)
{
  bool steep = std::abs(y1 - y0) > std::abs(x1 - x0);

  if(steep){
    rsSwap(x0, y0);
    rsSwap(x1, y1); }
  if(x0 > x1){
    rsSwap(x0, x1);
    rsSwap(y0, y1); }

  TCor dx = x1 - x0;
  TCor dy = y1 - y0;
  TCor gradient = dy / dx;
  if(dx == 0.0)
    gradient = 1.0;

  // handle first endpoint:
  int  xend  = roundToInt(x0);                     
  TCor yend  = y0 + gradient * (xend - x0);
  TCor xgap  = rfpart(x0 + 0.5f);
  int  xpxl1 = xend;                  // will be used in the main loop
  int  ypxl1 = ipart(yend);
  TCor fp    = fpart(yend);           // == yend-ypxl1
  if(steep){
    plot(ypxl1,   xpxl1, (1-fp) * xgap * color);
    plot(ypxl1+1, xpxl1,    fp  * xgap * color); } 
  else {
    plot(xpxl1, ypxl1,   (1-fp) * xgap * color);
    plot(xpxl1, ypxl1+1,    fp  * xgap * color); }
  TCor intery = yend + gradient;      // first y-intersection for the main loop

  // handle second endpoint:  
  xend      = roundToInt(x1);
  yend      = y1 + gradient * (xend - x1);
  xgap      = fpart(x1 + 0.5f);
  int xpxl2 = xend;                    // will be used in the main loop
  int ypxl2 = ipart(yend);
  fp        = fpart(yend);             // == yend-ypxl2
  if(steep){
    plot(ypxl2,   xpxl2, (1-fp) * xgap * color);
    plot(ypxl2+1, xpxl2,    fp  * xgap * color); }
  else {
    plot(xpxl2, ypxl2,   (1-fp) * xgap * color);
    plot(xpxl2, ypxl2+1,    fp  * xgap * color); }

  // main loop:
  int ip;
  if(steep){
    for(int x = xpxl1+1; x <= xpxl2-1; x++){
      ip = ipart(intery);
      fp = intery-ip;
      plot(ip,   x, (1-fp) * color);
      plot(ip+1, x,    fp  * color);
      intery += gradient; }}
  else{
    for(int x = xpxl1+1; x <= xpxl2-1; x++){
      ip = ipart(intery);
      fp = intery-ip;
      plot(x, ip,  (1-fp) * color);
      plot(x, ip+1,   fp  * color);
      intery += gradient; }}
}

