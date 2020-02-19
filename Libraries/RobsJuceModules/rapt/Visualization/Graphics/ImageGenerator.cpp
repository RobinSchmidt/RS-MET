
template<class T>
T squaredDistance(T x1, T y1, T x2, T y2)
{
  T dx = x2 - x1;
  T dy = y2 - y1;
  return dx*dx + dy*dy;
}

template<class T>
T distance(T x1, T y1, T x2, T y2)
{
  return sqrt(squaredDistance(x1, y1, x2, y2));
}
// move to somewhere else


template<class TPix, class TVal> 
TVal rsImageGenerator<TPix, TVal>::spiralRidge(TVal x, TVal y, TVal a, TVal p, TVal sign, 
  int profile, TVal exponent)
{
  // sanity check inputs:
  rsAssert(sign == +1 || sign == -1,     "sign must be +-1");
  rsAssert(profile >= 0 && profile <= 3, "invalid profile");

  // compute raw height:
  TVal r = sqrt(x*x + y*y);
  if(r == 0.0) return 0.0;           // avoid log-of-zero
  TVal t  = log(r) / a;              // parameter t for point on the spiral with radius r
  TVal xs = r * cos(sign * t + p);   // x on the spiral for the given t
  TVal ys = r * sin(sign * t + p);   // y on the spiral for the given t
  TVal d  = distance(xs, ys, x, y);  // distance of input point to point on the spiral
  TVal h  = pow(0.5*d/r, exponent);  // height

  // apply shaping of the height profile:
  if(profile == 2) return h;                        // 2: rectified sine (comes out by raw formula)
  if(profile == 3) return 1-h;                      // 3: inverted rectified sine
  h = asin(h) / (0.5*PI);                           // convert to triangular
  if(profile == 0) return h;                        // 0: triangular
  if(profile == 1) return 0.5*(sin(PI*(h-0.5))+1);  // 1: sinusoidal
  return 0;                                         // unknown profile
}
// optimize: the sqrt in the distance computation can be avoided: compute the distance-squared and 
// then use 0.5*exponent in the subsequent pow call
//
// The algo computes the distance of (x,y) to a point on the spiral that has the same radius as 
// (x,y). It happens that the height-profile (as function of radius for a given angle) comes out as
// a rectified sine shape (when the radius is used a x-axis and the x-axis is logarithmically 
// scaled)
// what if sign is not +-1? what if we use different factors for x- and y: 
//   xs = r * cos(wx * t + px); ys = r * sin(wy * t + py);
// ..in this case, the profile computations will very likely become invalid because the raw profile 
// is not a rectified sine anymore - yes - using, for example cos(2*..),sin(3*..) gives a nice 
// effect - one should increase the shrink-factor accordingly because it get denser otherwise
// for the original profile and the profile converted to a full sine, it makes visually not 
// qualitative difference, when we invert all color channels - for the triangular profile, it does.
// maybe to create audio-signals, we could use dx = (xs-x)/r; dy = (ys-y)/r; as left and right 
// channel signal - but what should the input be? we don't have x,y pixel coordinates as inputs but
// time instants
// 
// when we use d / r, the birghtness of the white ridges is independent for the distance to the 
// center - using a power with exponent < 1, we get a darkening effect towrd the center - but mybe 
// such an effect can be applied as post-processing: 
// circularDarkening(img, x, y, amount)
//   img(i,j) /= pow(r, amount)
//
// try tL = atan2(y,x) + 2*k*pi where k = floor(t0/(2*pi)), tR = tL + 2*pi, compute (xL,yL),(xR,yR)
// by the parametric spiral equations, compute distances dL,dR and use minimum
//
// -these are not the actual distances to the nearest points on the spiral but rather the distances 
//  to two concentric circles that approximate the spiral at the given angle - but they can be used 
//  as an initial estimate for computing the actual distance via netwon iteration - maybe this 
//  refinement can be made optional, controlled by a boolean parameter
//
// see:
// https://en.wikipedia.org/wiki/Logarithmic_spiral
// https://en.wikipedia.org/wiki/Archimedean_spiral