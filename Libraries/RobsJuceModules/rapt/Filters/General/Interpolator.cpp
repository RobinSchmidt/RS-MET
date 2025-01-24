
// construction/destruction:
/*
rsInterpolator::rsInterpolator()
{
  previousOutput      = 0;
  interpolationMethod = LINEAR;
}

rsInterpolator::~rsInterpolator()
{

}

// parameter settings:

void rsInterpolator::setInterpolationMethod(int newMethod)
{
  interpolationMethod = newMethod;
}

// others:

void rsInterpolator::reset()
{
  previousOutput = 0;
}
*/

/*

Ideas:

- Implement an interpolation formula based on sines and cosines. Let's assume the usual normalized 
  setting where we want to interpolate between x0 = 0 and x1 = 1 with prescribed values 
  y0,y1,y0',y1',y0'',y1'',... as we do in Hermite interpolation. Now make the ansatz:

    w = pi*x,  y(x) = a0 + a1*cos(w) + b1*sin(w) + a2*cos(2*w) + b2*sin(2*w) + ...

  i.e. we interpolate with sines and cosines that have an integer number of half-cycles between
  x = 0 and x = 1. Maybe we should include a linear term, too - but I'm not sure about that. Maybe
  we should use one cosine more than we use sines to make the number of constraints match the 
  number of coefficients. Or maybe impose an additional constraint such as the integral between 0 
  and 1 being whatever linear interpolation would give. I have once done this with cubic Hermite 
  interpolation to "enhance" it into a quartic interpolation scheme. The code should be somewhere 
  in the experiments. Or maybe the interpolant should exactly pass through the midpoint, i.e. 
  through (x,y) = (1/2, (y1+y2)/2). I don't know, if that improves it, though. But if we need an
  additional constraint, it shouldn't be too hard to come up with one that makes the function 
  somehow look "nicer". Maybe when we use sines and cosines of these frequencies, we can make sure
  that any aliasing frequencies will end up on 0 Hz and/or fs/2 and therefore be inaudible? But 
  no - actual aliasing will only occur, when the readout increment is != 1. Otherwise, we will just
  see a frequency dependent damping (or maybe boosting?). Plot the frequency responses of such 
  interpolators for various values of the fractional part of the position.


*/

