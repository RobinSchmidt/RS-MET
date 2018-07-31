#pragma once

// Functions for continuous logic, where the truth values are seen as probability of an proposition 
// to be true.

/** Probability that two independent events a and b occur both where pa is the probability for a to
occur and pb is the probability for b to occur. */
inline double rsAnd(double pa, double pb) { return pa*pb; }

/** Probability that event a occurs or event b occurs or both occur, assuming a and b are
independent. */
inline double rsOr(double pa, double pb) { return pa + pb - pa*pb; }

/** A class for representing continuous probabilities which are real numbers in the range [0..1]. 
The class allows for various computations
*/

/*
class Probability
{

public:





protected:


  double p; // the actual probability value


};
*/

