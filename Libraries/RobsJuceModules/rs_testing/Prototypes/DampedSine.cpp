



/*

ToDo:
-implement addition and multiplication operators
-implement functions for derivative and integral, if possible, i.e. if the result is again a sum of 
 damped sines (i think so, but i'm not sure). The definite integral is just a number, so it can be
 computed anyway. I think, the integral look like a constant minus a decay (RC loading curve) with a 
 sine-wiggle on top...but no - it's monotonic...at least in the case of the energy integral, which
 is a the integral of a damped sine multiplied by itsefl...we'll see...
-implement efficient computation of an array of values by means of recursion
 -maybe this function should accumulate into the array, i.e. use +=, because in practice, we may 
  want to work with sums of damped sines

rsDampedSineSum:
-Maybe keep the array of components sorted in some way - maybe sort according to frequency first 
 (ascending), then according to total energy of the envelope (descending), then according to 
 amplitude (descending), then according to phase (ascending).
-When the array is known to be sorted, we can also meaningfully imlement the == operator.
...but maybe that is not useful - we'll see...


Math for the multiplication:
  
  f(x) = a1 * exp(-d1*x) * sin(w1*x + p1)  *  a2 * exp(-d2*x) * sin(w2*x + p2)
       = a1*a2 * exp(-(d1+d2)*x) * sin(w1*x + p1) * sin(w2*x + p2)
 
Then, for the product of sines, we use the general result:

  sin(a*x+b)*sin(c*x+d) = ( cos((a-c)*x+(b-d)) - cos((a+c)*x+(b+d)) ) / 2

which can be verified via wolfram alpha. For the resulting frequency and phase, there is some 
freedom: we can always choose a positive or negative frequency and then compensate by adding or 
subtracting pi from the phase. When we flip the sign of the frequency, we need to compensate by a 
phase-shift of 180�. ToDo: maybe we should always enforce positive frequencies in the result?



*/