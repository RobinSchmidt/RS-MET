


/*

Ideas:

- Make classes similar to rsNormalizedSigmoids for functions of other qualitative kinds:

  - rsFadeFunctions: for fade-in/out and crossfade functions. See Notes/FadeInOut.txt

  - rsCycleFunctions: functions that have one cycle of an oscillation. See the experiments
    singleSineCycleWobbles(), multiSineCycleWobbles, etc.

  - rsCompressionFunctions: functions that look a bit like sigmoids but are unbounded. For example
    cube-root, asinh, odd-roots, inverses of x + x3, x + x^5, ...

  - rsExpansionFunctions: functions like sinh, x^3, x + x^3, ...

  - rsExplosionFunctions: anti sigmoids, functions with vertical asymptotes at -1 amd +1 and a zero 
    at 0. Examples: -x / (x^2 - 1), atanh(x), tan(x * pi/2), ...


How to create certain fetaures:

 - Zeros at any position z with order m can be created by introducing a factor (x-z)^m.

 - Poles of finite order m at any x-position p can be created by dividing by (x-p)^m

 - A pole of infinite order at x = 0 can be created by exp(1/x) +- exp(-1/x) where + vs - decides 
   if it's a symmetric or antisymmetric function.

 - I think, we can impose any desired asymptotic function g(x) to any function f(x) by means of
   mixing functions m_f(x), m_g(x) as  (m_f(x)*f(x) + m_g(x)*g(x)) / (m_f(x) + m_g(x)). When m_g
   increases asymptotically faster than m_f, the result will asymptotically look like g. The 
   division normalizes to weights for f and g to sum to 1.


*/