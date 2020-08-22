#ifndef romos_Interpolation_h
#define romos_Interpolation_h

//namespace romos
//{

//-------------------------------------------------------------------------------------------------

// The following functions are for computing an intermediate data value that is in between the 
// known data samples (interpolation)

/** Interpolates between the points (xL, yL), (xR, yR) that represent x- and y-coordinates to the 
left and right to the point x at which we want to obtain an interpolated data value. Assuming 
yR > yL, the curve drawn between the two points is an exponential growth function when "shape" is 
greater than 0, an (up-down flipped) exponential decay function when shape is less than 0 and a 
straight line when shape equals zero. Analogously, when yR < yL, the shape is an exponential decay 
for shape < 0 and a flipped exponential growth for shape > 0. The absolute value of "shape" 
determines how many exponential time constants occur between the endpoints of the curve. Values of
shape < 0 draw curves that resemble the behaviour of analog (RC-type) envelope generators.
References:
 -Moore: Elements of Computer Music, page 184 */
double interpolateExponentially(double x, double xL, double yL, double xR, double yR, double shape);

//}

#endif
