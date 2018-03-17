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
  -and(a,b,c) should reduce to the form above when c=0, 
   when: c=1: and(a=1,b=1,1) = 0, and(a=1,b=0,1) = 0
  -when a formula for "and" is established, a formual for "or" can be derivedn from 
   or(a,b) = not(not(a),not(b))
  -maybe it should also output two new correlation values for (a,result), (b,result) as side 
   results

-3 valued logic: true,false,unknown
 T and U = U, T or U = T, F and U = F, F or U = U, U and U = U, U or U = U
 the rest stays the same (F and T = F etc.)
 -can be modelled by probabilistic logic by using 0.5 as "unknown" - maybe by clamping unknown
  results to 0.5
*/

#endif