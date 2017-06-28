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

void ParticleBouncer::getSampleFrame(double &x, double &y)
{
  // preliminary - later, we need to check for intersections with the enclosure and compute 
  // reflections:

  // update current particle position:
  xc += dx;
  yc += dy;

  // assign outputs:
  x = xc;
  y = yc;
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