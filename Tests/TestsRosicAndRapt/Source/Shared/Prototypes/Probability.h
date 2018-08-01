#pragma once

// maybe rename file to Logic.h


//=================================================================================================

/** A class for representing a value in three valued logic. The truth value of a proposition can be 
either false (0,F), true (1,T) or unknown (2,U). */

class TriLogic
{

  /** Negation. not(U) = U. 
  in:   T  F  U
  out:  F  T  U  */
  static char not(char a);

  /** Logical and. and(T, U) = U, and(F, U) = F, and(U, U) = U. 
     T F U
  T  T F U
  F  F F F
  U  U F U   */
  static char and(char a, char b);

  /** Logical or. or(T, U) = T, or(F, U) = U, or(U, U) = U 
     T F U
  T  T T T
  F  T F U
  U  T U U   */
  static char or(char a, char b);


  inline void setTrue() { value = 1; }

  inline void setFalse() { value = 0; }

  inline void setUnknown() { value = 2; }

  inline bool isTrue() { return value == 1; }

  inline bool isFalse() { return value == 0; }

  inline bool isUnknown() { return value == 2; }

  inline bool isPossible() { return isTrue() || isUnknown(); }

  inline bool isQuestionable() { return isFalse() || isUnknown(); }


protected:

  char value = 0;  // 0: false, 1: true, 2: unknown

};
// todo: maybe explore other binary functions besides and/or and their algebraic properties. Maybe
// addition and multiplication modulo 3 could be useful and/or have interesting algebraic 
// properties (algebraic structure should be isomorphic to (Z_3,+,*) )


// Functions for continuous logic, where the truth values are seen as probability of an proposition 
// to be true.



inline double rsNot(double p) { return 1 - p; }

/** Probability that two independent events a and b occur both where pa is the probability for a to
occur and pb is the probability for b to occur. */
inline double rsAnd(double pa, double pb) { return pa*pb; }

/** Probability that event a occurs or event b occurs or both occur, assuming a and b are
independent. */
inline double rsOr(double pa, double pb) { return pa + pb - pa*pb; }

/** == or(and(a,not(b)),and(b,not(a))) */
inline double rsXor(double pa, double pb) { return pa + pb - pa*pb*(3-pa-pb+pa*pb); }

//=================================================================================================

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

