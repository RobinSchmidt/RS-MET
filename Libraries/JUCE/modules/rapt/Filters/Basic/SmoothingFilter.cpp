template<class TSig, class TPar>
rsSmoothingFilter<TSig, TPar>::rsSmoothingFilter()
{
  y1.resize(1);
  reset();
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::setTimeConstantAndSampleRate(TPar timeConstant, TPar sampleRate)
{
  decay = sampleRate * timeConstant;
  updateCoeff();
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::setOrder(int newOrder)
{
  order = rsMax(1, newOrder);

  y1.resize(order);

  reset();
  // todo: if newOrder > oldOrder, init only the vector values in y1 above oldOrder-1 to 0
  // not all of them

  updateCoeff();
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::reset()
{
  for(int i = 0; i < order; i++)
    y1[i] = 0;
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::updateCoeff()
{
  //coeff = exp(-order/decay); // amounts to divide the time-constant by the order
  coeff = exp(-1/decay); // no scaling for test purposes
}


/*

ToDo:
-make an envelope follower based on this filter for use in dynamics processors
-the env-follower should be availbale as modulation source in chainer

Trying to work out a formula for where the step response goes through 0.5. The impulse response
of the 1st order filter with unit time constant is given by:
h1(t) = e^(-t) for t >= 0 and 0 for t < 0
the impulse response of the 2nd order filter is given by the convolution of h1 with itself:
h2(t) = conv(h1(t), h1(t)) = integrate e^(-(t-T))*e^(-t) dT from T=0 to t 
this can be plugged to wolfram alpha and gives:
h2(t) = e^(-2*t) (e^t-1)
the 3rd oder impulse response is given by the convolution of that h2(t) with h1(t):
h3(t) = conv(h1(t), h2(t)) = e^(-(t-T)) * e^(-2*t) (e^t-1) dT from T=0 to t 
again, wolfram alpha can solve this as:
h3(t) = e^(-3 t) (-1 + e^t)^2
continuing this way, we get:
h4(t) = conv(h1(t), h3(t)) = e^(-(t-T)) * e^(-3 t) (-1 + e^t)^2 dT from T=0 to t 
h4(t) = e^(-4 t) (-1 + e^t)^3
the pattern that emerges is:
hN(t) = e^(-N t) (-1 + e^t)^(N-1)
OK, we have a nice formula for the impulse response - now for the step response, we must convolve
the impulse response with the unit step. this just amounts to integrating the impulse response
from 0 to t (i.e. the step-response is the running sum of the impulse response). Doing the
integral for N=4 (integrate e^(-4 t) (e^t - 1)^3 dt from t=0 to t) gives:
s4(t) = 1/4 e^(-4 t) (e^t - 1)^4
and in general:
sN(t) = 1/N e^(-N t) (e^t - 1)^N
...however, when plotting these, they only approach 1/N instead of 1, so this 1/N factor seems 
wrong (why?). should it not be just: 
sN(t) ?= e^(-N t) (e^t - 1)^N
assuming that is right, we set that expression equal to our target value a = 1/2 and solve for t:
t = -log(1 - pow(a, 1/N))
...that should be the time instant, where step response goes through a...check that



FUCK - NO this is all wrong - it has to be:

h1(t) = e^(-t) for t >= 0 and 0 for t < 0
h2(t) = conv(h1(t), h1(t)) = integrate e^(-(t-T))*e^(-T) dT from T=0 to t 
h2(t) = t * e^-t
h3(t) = integrate e^(-(t-T)) * T*e^(-T) dT from T=0 to t 
h3(t) = t^2/2 * e^-t
h4(t) = integrate e^(-(t-T)) * (T^2/2)*e^(-T) dT from T=0 to t
h4(t) = t^3/6 * e^-t
h5(t) = integrate e^(-(t-T)) * (T^3/6)*e^(-T) dT from T=0 to t
h5(t) = t^4/24 * e^-t
h6(t) = integrate e^(-(t-T))*(T^4/24)*e^(-T) dT from T=0 to t
h6(t) = t^5/120 * e^-t
the general expression for the impulse response is then:

hN(t) = e^-t * t^(N-1) / (N-1)!

step responses are the running integral (sum) of the impulse responses:

sN(t) = integrate hN(T) dT from T=0 to t
s1(t) = integrate             e^(-T) dT from T=0 to t = 
s2(t) = integrate T         * e^(-T) dT from T=0 to t
s3(t) = integrate (T^2/2)   * e^(-T) dT from T=0 to t
s4(t) = integrate (T^3/6)   * e^(-T) dT from T=0 to t
s5(t) = integrate (T^4/24)  * e^(-T) dT from T=0 to t
s6(t) = integrate (T^5/120) * e^(-T) dT from T=0 to t

sN(t) = 1 - e^-t * ( sum_{k=0}^{N-1} (t^k)/(k!) )
      = 1 - e^-t * (1 + t + t^2/2 + t^3/6 + t^4/24 + t^5/120 + ... + t^(N-1)/(N-1)! )

...damn! that expression is harder to handle - we need to set it equal to our target value 
a = 1/2 and solve for t. I think, it's not explicitly soluble for t, so maybe we need a root 
finder. I think, to do this, we need to drag over the polynomial and root-finder class

*/