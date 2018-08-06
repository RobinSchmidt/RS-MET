template<class T>
void rsTriSawOscillator<T>::updateTriSawCoeffs()
{
  // range variables - make settable members later:
  //T min = 0;
  //T max = 1;
  // hmm...but actually we need thme to be 0 and 1 because the waveshaping expects its input in
  // this range - if the range should really be different we may scale/shift as last step...

  //// coeffs:
  //a0 = min;
  //a1 = (max-a0) / h;
  //b0 = (max-h*min) / (1-h);
  //b1 = min-b0;

  // coeffs:
  a0 = 0;  // can actually be optimized away
  a1 = T(1) / h;
  b0 = T(1) / (1-h);
  b1 = -b0;

  // todo: we need to catch cases when h = 0 or 1-h = 0
  // if(h < eps)
  //   a1 = 0;
  // if(1-h < eps)
  //   b1 = 0;
  // or something
}

// for shapes, see:
// rational:    https://www.desmos.com/calculator/ql1hh1byy5 (a*x+x)/(2*a*x-a+1)
// exponential: https://www.desmos.com/calculator/vtymebgkrr (1-exp(a*x))/(1-exp(a))
// sinh(a*x)/sinh(a), tanh(a*x)/tanh(a)

// maybe try something based on:
// (a*exp(x) + b*exp(-x) + p) / (c*exp(x) + d*exp(-x) + q)
// tanh: a=1,b=-1,p=0,c=1,d=0,q=0
// sinh:
//  exp: 
// logistic: 1/(1+exp(-x))

// weirdo: tanh(a*x)/(1+exp(-x))

// maybe have a look at these:
// https://en.wikipedia.org/wiki/Gompertz_function
// https://en.wikipedia.org/wiki/Generalised_logistic_function

// or try: f(x) = (a*x + b) / (c*x + d) and to keep equations simple, consider the range -1..+1 and
// require f(-1) = -1, f(+1) = +1, f(0) = p, f'(0) = s where p and s are adjustable parameters, p
// corresponds roughly to "a" in the rational formula (can we make it, such that when s=0, it 
// reduces to the old formula?)
// ...but playing with the parameters, it doesn't seem to be capabable of a sigmoid shape
// anyway: f(1)=1 -> a+b = c+d, f(-1) = -1 -> c-d = b-a, f(0) = p -> d*p = b
// f'(x) = (a*d-b*c) / (c*x+d)^2, f'(0) = (a*d-b*c)/d^2
// f'(0) = s ->

// maybe try f(x) = (a*x + b*x^3) / (c*x + d) this seems to be able to do sigmoids
// http://www.wolframalpha.com/input/?i=(a*x+%2B+b*x%5E3)+%2F+(c*x+%2B+d)
// f'(x) = (a*d+b*x^2*(2*c*x+3*d)) / (c*x+d)^2
// f'(0) = a*d/d^2 = s
// f(0)  = 0 ...damn - we loose the ability to adjust the y-value at the midpoint - we need both

// how about f(x) = a0 + a1*x + a3*x^3 / (b0 + b1*x)
// should allow to adjust both...maybe also try an a2*x^2 term insead of a3*x^3 ...or maybe both?
// or maybe try a 2nd order denominator, too - maybe a general biquadratic formula?