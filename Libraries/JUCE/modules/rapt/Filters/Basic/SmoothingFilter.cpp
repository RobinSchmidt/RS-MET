template<class TSig, class TPar>
const TPar rsSmoothingFilter<TSig, TPar>::maxAsym = 2.0;

template<class TSig, class TPar>
bool rsSmoothingFilter<TSig, TPar>::tableIsFilled = false;

template<class TSig, class TPar>
rsMatrix<TPar> rsSmoothingFilter<TSig, TPar>::tauScalers; 



template<class TSig, class TPar>
rsSmoothingFilter<TSig, TPar>::rsSmoothingFilter()
{
  if(!tableIsFilled)
    createTauScalerTable();
  y1.resize(1);
  coeffs.resize(1);
  reset();
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::setTimeConstantAndSampleRate(TPar timeConstant, TPar sampleRate)
{
  decay = sampleRate * timeConstant;
  updateCoeffs();
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::setOrder(int newOrder)
{
  order = rsMax(1, newOrder);

  y1.resize(order);
  coeffs.resize(order);

  reset();
  // todo: if newOrder > oldOrder, init only the vector values in y1 above oldOrder-1 to 0
  // not all of them

  updateCoeffs();
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::setShape(int newShape)
{
  shape = newShape;
  updateCoeffs();
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::setShapeParameter(TPar newParam)
{
  shapeParam = newParam;
  updateCoeffs();
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::reset()
{
  for(int i = 0; i < order; i++)
    y1[i] = 0;
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::updateCoeffs()
{
  TPar tmp;

  if(shape == FAST_ATTACK)
  {
    for(int i = 0; i < order; i++)
    {
      tmp  = decay / (TPar) pow(i+1, shapeParam); // scaled decay time-constant
      //if(i > 0)
      //  tmp *= (shapeParam+1);
      coeffs[i] = exp(-order/tmp);      // tau[n] = tau[0] / n^p // p == shapeParam
    }
    // maybe try, if it responds different to modulations of the time-constants are in reverse
    // order (from short to long instead of long to short)
  }
  else
  {
    // all filter stages use the same time-constant
    tmp = exp(-order/decay); // amounts to divide the time-constant by the order
    for(int i = 0; i < order; i++)
      coeffs[i] = tmp;
  }

  //// test:
  //TPar scaler = 1;
  //if(order > 1)
  //  scaler = 1.0 / (order-1);
  //TPar scaledDecay = scaler * decay;
  //coeff = exp(-1/scaledDecay);

  //coeff = exp(-1/decay); // no scaling for test purposes
}

template<class TSig, class TPar>
void rsSmoothingFilter<TSig, TPar>::createTauScalerTable()
{
  y1.resize(maxOrder);
  coeffs.resize(maxOrder);
  reset();
  tauScalers.setSize(maxOrder, numAsyms);
  tauScalers.setAllValues(1);

  // something to do...

  y1.resize(maxOrder);
  coeffs.resize(maxOrder);
  setOrder(1);
  tableIsFilled = true;
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
s1(t) = integrate             e^(-T) dT from T=0 to t
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

wolfram alpha can produce numerical results for this equation:

1 - e^-t * (1 + t + t^2/2! + t^3/3! + t^4/4! + t^5/5! + t^6/6! + t^7/7!) = 1/2

maybe when using the result of this computation in the analog domain, a step-invariant digital
filter approximation should be used - it's explained here:
http://homes.esat.kuleuven.be/~maapc/Sofia/slides_chapter8.pdf
because it doesn't really match the 0.5-crossing of the digital filter

...but maybe we should just create a table by reading off the time-instants where the curves cross
1/2 from the step responses

or maybe we could normalize it in such a way that the inflection point (zero of derivative of 
hN(t)) is always at the same time instant. The derivative of is:
hN'(t) = ((N-1)*t^(N-2) * e^-t - t^(N-1) * e^-t) / (N-1)!  ...verify that formula
       = e^-t * t^(N-2) * (N-1-t) / (N-1)!

Or try a lowpass chain where each stage has a different cutoff frequency. The individual stages 
have impulse responses:
h1(t) = e^(-a_1*t), h2(t) = e^(-a_2*t), h3(t) = e^(-a_3*t), ..., hN(t) = e^(-a_N*t)
the accumulated impulse reponses after each stage n, say gn, are given by: 
g1(t) = h1(t), g2(t) = conv(h2(t), g1(t)), g3(t) = conv(h3(t) g2(t)), ....
so
g2(t) = integrate e^(-a_2 (t-T)) e^(-a_1 T) dT for T=0 to t
      = (e^(-a_1*t) - e^(-a_2*t)) / (a_2 - a_1)
g3(t) = integrate e^(-a_3 (t-T)) (e^(-a_1*T) - e^(-a_2*T)) / (a_2 - a_1) dT for T=0 to t
...hmm..the expressions seems to get more complicated at an exponential rate (each is twice
as complicated as the one before)...impractical to work with

...maybe try the bandwidth scaling formula?

...maybe instead of using the convolution integrals for the analog counterpart, we may use the
convolution sums for the actual digital version? If that works out, it will be even more precise.
We take the difference equation:
y[n] = x[n] + c * (y[n-1] - x[n]) = (1-c) * x[n] + c * y[n-1]
The impulse response is:
h[n] = (1-c) * c^n
We assume not necessarily equal coeffs for each stage, so we index the coeffients c like 
c1, c2, c3, ..., cN. The impulse responses of the individual stages are then given by:
h1[n] = (1-c1) * c1^n, h2[n] = (1-c2) * c2^n, ..., hN[n] = (1-cN) * cN^n
The accumulated impulse response of each m-th stage is given by the convolution sum of the 
accumulated impulse response of the (m-1)th stage and the impulse reponse of the m-th stage:
g1[n] = h1[n]
g2[n] = sum_k=0^n g1[k] * h2[n-k] 
g3[n] = sum_k=0^n g2[k] * h3[n-k] 
etc. The step response of the N-th stage is then given by the running sum of its impulse response:
sN[n] = sum_k=0^n gN[k]
If we can work out a nice closed form expression for sN[n] which is soluble for n, we are
done.

pass to wolfram alpha:
sum (1-c) c^k for k=0 to n (just for test)
sum (1-c_1) c_1^k (1-c_2) c_2^(n-k) for k=0 to n

g1[n] = (1-c_1) c_1^n
g2[n] = sum_(k=0)^n (1-c_1) c_1^k (1-c_2) c_2^(n-k) 
      = ((c_1-1) (c_2-1) (c_1^(n+1) - c_2^(n+1)))/(c_1-c_2)
g3[n] = sum_(k=0)^n ((c_1-1) (c_2-1) (c_1^(k+1) - c_2^(k+1)))/(c_1-c_2)  (1-c_3) c_3^(n-k)
      = ...no solution

*/