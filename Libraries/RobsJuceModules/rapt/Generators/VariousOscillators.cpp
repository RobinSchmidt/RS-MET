

template<class T>
void rsTriSawOscillator<T>::updateTriSawCoeffs()
{
  //// range variables - make settable members later:
  //T min = -1;
  //T max = +1;
  //// hmm...but actually we need them to be -1 and 1 because the waveshaping expects its input in
  //// this range - if the range should really be different we may scale/shift as last step...

  //// coeffs:
  //a0 = min;
  //a1 = (max-a0) / h;
  //b0 = (max-h*min) / (1-h);
  //b1 = min-b0;

  // simpler:
  a0 = -1;
  a1 = 2 / h;
  b0 = (1+h)/(1-h);
  b1 = -1-b0;

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

// how about f(x) = (a0 + a1*x + a3*x^3) / (b0 + b1*x)
// should allow to adjust both...maybe also try an a2*x^2 term insead of a3*x^3 ...or maybe both?
// or maybe try a 2nd order denominator, too - maybe a general biquadratic formula?
// looks like it's flexible enough: https://www.desmos.com/calculator/4dxeptvth3
// question would be, how convenient it is to control
// http://www.wolframalpha.com/input/?i=(a+%2B+b+x+%2B+c+x%5E2)+%2F+(p+%2B+q+x+%2B+r+x%5E2)

/*
for the biquadratic formula, giving this to sage:

var("x a b c p q r s t")
f(x)  = (a + b*x + c*x^2) / (p + q*x + r*x^2)
fp(x) = diff(f(x), x)
e1 = f(-1) == -1
e2 = f( 1) ==  1
e3 = f( 0) ==  t
e4 = fp(0) ==  s
solve([e1,e2,e3,e4],[a,b,c,p,q,r]) 

simplified:

var("a b c p q r s t")
e1 = a+b+c == p+q+r
e2 = a-b+c == q-p-r
e3 = t == a/p
e4 = s == (b*p-a*q)/p^2
solve([e1,e2,e3,e4],[a,b,c,p,q,r]) 

..yep, gives the same result

// first solution has 2 params and looks like this:
https://www.desmos.com/calculator/tv3w6uquxr
// not useful
// write rational map as ((T+1)x)/(2*T*x-T+1): https://www.desmos.com/calculator/lybdjuzzln
// where T is the "tension" in -1..+1 - we see that when the sigmodity S is zero, we need to have
// a=0, b=(T+1), c=0, p=1-T, q=2*T, r=0
//...ahh...but we need to re-derive it so the function is the range -1..+1, not 0..1



produces a solution that contains new variables r1,r2,r3,r4 - what are they? roots? 
but of which equation? probably a quartic because there are 4 roots? here is a simpler example in 
which just r1 occurs

var("x a b c d s")
f(x)  = a + b*x + c*x^3 + d*x^5
fp(x) = diff(f(x), x)
e1 = f(-1) == -1
e2 = f( 1) ==  1
e3 = f( 0) ==  0
e4 = fp(0) ==  s
solve([e1,e2,e3,e4],[a,b,c,d]) 

oh - no - here (just below the center)
http://doc.sagemath.org/html/en/reference/calculus/sage/symbolic/relation.html
it says:
If there is a parameter in the answer, that will show up as a new variable. In the following 
example, r1 is a real free variable (because of the r):

sage: forget()
sage: x, y = var('x,y')
sage: solve([x+y == 3, 2*x+2*y == 6],x,y)
[[x == -r1 + 3, y == r1]]

Especially with trigonometric functions, the dummy variable may be implicitly an integer 
(hence the z):

or maybe compose y = g(x) = (1-t)*x + t*x^3 with h(y) = (a*y + b)/(c*y + d)
the first function gives the s-shape, the 2nd the tension
use h(y) = (y+p)/(p*y+1) where -1 < p < +1 is the parameter

var("a b c d x p")
f(x) = (a*x + b) / (c*x + d)
e1 = f(1)  ==  1
e2 = f(-1) == -1
e3 = f(0)  == p
solve([e1,e2,e3],[a,b,c])

maybe use clip(x / (1-s), -1, +1) where s is a "squarishness" parameter. when 1, the denominator is
0, so the output just alternates between -1 and +1 according to the numerator, so that would give a 
square wave (or general pulse-wave, if h is != 0.5) - adds another degree of flexibility

*/