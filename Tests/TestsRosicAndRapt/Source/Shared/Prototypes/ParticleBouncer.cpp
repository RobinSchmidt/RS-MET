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


}

/*

General Idea:
Let a particle bounce around inside an enclosure in 2D space and use its x- and y-coordinates as
signal (for left and right channel, or mid and side or whatever). The enclosure can be defined 
geometrically for example as an ellipse or a polygon. Maybe there can be also other shapes inside
the enclosure that act as obstacles for the particle.

Implementation:
At each sample, our virtual particle has a position given by (x,y) and a velocity vector (vx, vy).
To update the particle, we move it a bit along a line segment, i.e. we compute: x += vx, y += vy.
However, when the particle hits the boundary of the enclosure (or one of the obstacles), it gets
bounced back according to the reflection law: the angle of incidence should equal the angle of
reflection. The absolute value of the velocity (i.e. the speed) should stay the same.


Enclosure/Obstacle Shapes:


Relevant Equations:

-Line:
 -implicit equation:
 -parameteric equation:
 -angle between two lines:

-Ellipse
 -implicit equation: x^2/a^2 + y^2/b^2 = 1
 -parametric equation: x(t) = a * cos(t), y(t) = b * sin(t), t in 0..2*PI
 -area: A = PI*a*b
 -tangent at x0,y0: x*x0/a^2 + y*y0/b^2 = 1

-Polygon:
 ...

Notes:
-I think, the "frequency" of the resulting waveform should be proportional to the particle speed 
 and inversely proportional to the total area inside which the particle is free to bounce around.
 This fact can be used to keep the frequency constant when the shapes change.


*/