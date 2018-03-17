#ifndef Ideas_h
#define Ideas_h


/*
-continuous/probabilistic logic functions that take in real numbers a,b (probabilities):
 -not(a)    = 1-a
 -and(a, b) = a*b;
 -or( a, b) = 1 - (1-a)*(1-b) = a + b - a*b
  here, (1-a)*(1-b) is the probability that a and b are both false
  nand(a, b) = 1-a*b
 -xor(a, b) = or(and(a,not(b)),and(b,not(a))) = a + b - a*b(3-a-b+a*b)
 -maybe these functions can be extended to take a 3rd argument c for the correlation between 
  a and b
  -and(a,b,c) should reduce to the form above when c=0: and(a,b,0) = a*b
   when c = +1: and(a=1,b=0,+1) = and(a=0,b=1,+1) = 0
   when c = -1: and(a=1,b=1,-1) = and(a=0,b=0,-1) = 0
   ...maybe we can find more conditions that must hold and derive a general formula for and(a,b,c)
  -when a formula for "and" is established, a formula for "or" can be derived from 
   or(a,b) = not(and(not(a),not(b)))
  -maybe it should also output two new correlation values for (a,result), (b,result) as side 
   results
 -maybe a class for such variables can be written, using operators "*" for "and" and "+" for "or"
  -can we then also derive meaningful formulas for "/" and "-"?
 -...maybe using p,q instead of a,b is a nicer convention

-3 valued logic: true,false,unknown
 T and U = U, T or U = T, F and U = F, F or U = U, U and U = U, U or U = U
 the rest stays the same (F and T = F etc.)
 -can (?) be modelled by probabilistic logic by using 0.5 as "unknown" - maybe by clamping unknown
  results to 0.5

-use that probabilistic logic in a cellular automaton - the new activation of each cell is a 
 logical combination of its and neighbour activations
 -read out the celluar automaton one cell per sample, after each cycle compute new cell activations
 -have a global "fuzziness" parameter that crossfades all cell-activtions between hard 0/1 and 0.5
 -maybe lowpass filter cell contents during readout
 -maybe use a 2D automaton with various readout rules (scan lines, diagonal, hilbert-curve, etc.)
 -needs on-the-fly sample-rate conversion - read out the automaton always at its "natural" rate and
  then convert to target sample-rate
 -maybe neighbour correlations can be computed on the fly and be fed into the and/or operations 
  (maybe scaled by a factor)
 -the cells could be initialized with random values, maybe with an optional "balancing" of 
  0s and 1s: just count the 1s, and if there are n too many, change n randomly chosen 1-cells to 0
  (or the other way around)
*/

#endif