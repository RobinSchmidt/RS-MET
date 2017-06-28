#include "ParticleBouncer.h"


ParticleBouncer::ParticleBouncer()
{
 reset();
}

ParticleBouncer::~ParticleBouncer()
{

}

void ParticleBouncer::reset()
{
  xc = x0; 
  yc = y0;
}

// audio processing:

double ParticleBouncer::getLineEllipseIntersectionParameter(double x, double y, double dx, 
  double dy, double a2r, double b2r)
{
  // Coeffs for the quadratic equation a*x^2 + b*x + c = 0:
  double a, b, c;
  a =    dx*dx*a2r + dy*dy*b2r;
  b = 2*(x*dx*a2r  + y*dy*b2r);
  c =    x*x*a2r   + y*y*b2r - 1;

  // Solutions: t_1,2 = (-b +- sqrt(b^2-4*a*c)) / (2*a):
  double d = b*b - 4*a*c;         // discriminant
  d = sqrt(max(d, 0));            // max, just in case numerical error leads to negative d
  double t1 = (-b + d) / (2*a);   // the solution we are interested in
  //double t2 = (-b - d) / (2*a); // the other solution which is irrelevant for us

  return t1;

  // The formula was derived by setting :
  // (1) xi = x + t*dx = a * cos(2*PI*t)   
  // (2) yi = y + t*dy = b * sin(2*PI*t)  t in 0..1
  // divide (1) by a, (2) by b, square both equations and add them. That gives a right-hand-side
  // of sin^2(..) + cos^2(..) = 1 and leads to a quadratic equation for t. We are interested in a
  // t between 0 and 1. xi, yi are the x,y coordinates of the intersection point given by
  // xi = x + ti*dx, yi = y + ti*dy where ti is the value of t where the intersection occurs.
}

void ParticleBouncer::getSampleFrame(double &x, double &y)
{
  // assign outputs:
  x = xc;
  y = yc;

  // ...Update current particle position...

  // Compute new (tentative) x,y coordinates:
  double xn, yn;
  xn = xc + dx;
  yn = yc + dy;

  // Check, if new coordinates are inside enclosure by plugging them into the implicit ellipse 
  // equation x^2/a^2 + y^2/b^2 - 1 = 0. If it's not 0, let's call the right hand side d. This 
  // is proportional to the signed distance between the particle and the perimeter of the 
  // ellipse (what's the proportionality constant?):

  double a2r = 1 / (a*a), b2r = 1 / (b*b); // maybe make members
  double d   = xn*xn*a2r + yn*yn*b2r - 1;
  if(d > 0.0)  
  {
    // xn,yn is outside ellipse - we need to reflect...

    // Compute intersection point between current line segment and ellipse:
    double ti, xi, yi;
    ti = getLineEllipseIntersectionParameter(xc, yc, dx, dy, a2r, b2r);
    xi = xc + ti*dx;
    yi = yc + ti*dy;

    // for debug - check that xi,yi is indeed on the ellipse:
    //double err = xi*xi*a2r + yi*yi*b2r - 1; // should be 0 up to roundoff

    int dummy = 0;

    // Compute tangent to the ellipse at intersection point xi, yi:
    // ...
  }


  // update current coordinates:
  xc = xn;
  yc = yn;
}

/*

General Idea:
Let a particle bounce around inside an enclosure in 2D space and use its x- and y-coordinates as
signal (for left and right channel, or mid and side or whatever). The enclosure can be defined 
geometrically for example as an ellipse or a polygon. Maybe there can be also other shapes inside
the enclosure that act as obstacles for the particle. The user hase some means to configure the 
shape of the enclosure and obstacles and can give an initial position and velocity for the 
particle (if the initial position happens to be inside some obstacle, just ignore the collision
with the obstacle boundary - we imagine these boundaries a semi-permeable: things can get out but 
not in, except, of course, for the enclosure boundary).

Implementation:
At each sample, our virtual particle has a position given by (x,y) and a velocity vector (vx, vy).
To update the particle, we move it a bit along a line segment, i.e. we compute: x += vx, y += vy.
However, when the particle hits the boundary of the enclosure (or one of the obstacles), it gets
bounced back according to the reflection law: the angle of incidence should equal the angle of
reflection. The absolute value of the velocity (i.e. the speed) should stay the same.


Enclosure/Obstacle Shapes:


Relevant Equations:

-Line:
 -implicit equation: A*x + B*y + C = 0
 -2-point equation: (y-y1)/(x-x1) = (y2-y1)/(x2-x1)
 -point/direction equation: y-y1 = m*(x-x1) where m = (y2-y1)/(x2-x1)
 -parametric equation:
 -explicit equation: y = m*x + b
 -angle between two lines:

-Ellipse
 -implicit equation: x^2/a^2 + y^2/b^2 = 1
 -parametric equation: x(t) = a * cos(t), y(t) = b * sin(t), t in 0..2*PI
 -explicit equation: y = +- (b/a) * sqrt(a^2 - x^2)
 -area: A = PI*a*b
 -tangent at x0,y0: x*x0/a^2 + y*y0/b^2 = 1

-Polygon:
 ...

Notes:
-I think, the "frequency" of the resulting waveform should be proportional to the particle speed 
 and inversely proportional to the total area inside which the particle is free to bounce around.
 This fact can be used to keep the frequency constant when the shapes change.
-Maybe we can apply some geometric transformations (rotation, stretch, shear, etc.) on the shapes 
 and/or the particle position/velocity before applying the formulas (for line, ellipse, etc.).
 -for this, it may be useful to know that the determinant of the transformation matrix gives the
  (signed) area change for any shape


*/