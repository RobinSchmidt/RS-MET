#pragma once

/** Given a filter-design specification in terms of polynomial coefficients, this function analyzes
whether or not the filter satisfies complementariness conditions. These conditions mean that if we 
obtain a "complementary" signal by subtracting the filter's output from the input, the resulting 
frequency response of the difference-filter is a mirror image of the frequency response of the 
given filter. */
bool isComplementary(const RAPT::rsFilterSpecificationBA<double>& specBA);
// todo: maybe rename the condition...or maybe not
// Actually, that condition of being a mirror image only applies to lowpass and highpass filters. 
// For bandpass filters, the complementary filter will be a bandreject filter and in this case, the
// condition for complementarity may be a bit more complicated to specify. Maybe this function 
// should be renamed to reflect that. Maybe isComplementaryLowpass - we usually think of the 
// lowpass as the actual filter and the highpass is obtained by subtraction although it could be
// done the other way around as well.

/** Makes some plots, etc. */
bool analyzeComplementaryFilter(const RAPT::rsFilterSpecificationBA<double>& specBA);


// The actual prototype designs that satisfy the mirror-image conditions:

// Probably useful:
RAPT::rsFilterSpecificationBA<double> complementaryLowpass1p1z();
RAPT::rsFilterSpecificationBA<double> complementaryLowpass2p2z();

// Probably useless:
RAPT::rsFilterSpecificationBA<double> complementaryLowpass2p3z();
RAPT::rsFilterSpecificationBA<double> complementaryLowpass4p4z1t(); // 1t: one tweakable
RAPT::rsFilterSpecificationBA<double> complementaryLowpass4p4z();   // 2t
RAPT::rsFilterSpecificationBA<double> complementaryLowpass4p5z();
//RAPT::rsFilterSpecificationBA<double> complementaryLowpass3p3z();


/** Used internally by zLowpassToLowpass, zLowpassToHighpass. */
void zMapFirstOrder(RAPT::rsFilterSpecificationZPK<double>& zpk, double g, double c, 
  std::complex<double> zNorm);


/** Transforms the given prototype lowpass filter specification with cutoff wp (in radians) into 
another lowpass specification with some target cutoff wt (also in radians) using the 
Constantinides transform formula. */
RAPT::rsFilterSpecificationBA<double> zLowpassToLowpass(
  const RAPT::rsFilterSpecificationBA<double>& spec, double wp, double wt);
// Maybe this function should go into class rsPoleZeroMapper - the functionality fits there well. But
// there, it should operate on raw arrays for production use. We can the keep this function here as
// convenience function that callsthe other one

RAPT::rsFilterSpecificationBA<double> zLowpassToHighpass(
  const RAPT::rsFilterSpecificationBA<double>& spec, double wp, double wt);

RAPT::rsFilterSpecificationBA<double> zLowpassToBandpass(
  const RAPT::rsFilterSpecificationBA<double>& spec, double wp, double wl, double wu);

RAPT::rsFilterSpecificationBA<double> zLowpassToBandreject(
  const RAPT::rsFilterSpecificationBA<double>& spec, double wp, double wl, double wu);



// ToDo: 
// -Maybe wrap into a class mainly to have a scope under which we can do a
//    using FilterSpecBA = RAPT::rsFilterSpecificationBA<double>
//  so we don't have to write it out that often. Maybe a class rsComplementaryFilterDesigner
// -Obtain frequency warped versions via Constantinides formulas. The functions above design
//  prototypes with their cutoff frequency being a quarter of the sample rate (halfband filters).
//  we need the Constantinides transforms to eventually make the split frequency tweakable by the
//  user.
// -For the useful 2p2z design, try using the LP->BP transform. I hope that this will turn this 
//  LP/HP pair into an equally useful BP/BR pair (each of which with 4 poles and zeros).
// -Add complementaryFilter2p4z, 4p4z maybe try to place 3 zeros at z = -1
// -Try another design approach: instead of manually placing poles and zeros in the z-plane by 
//  trial and error, start with a (digital) allpass and obtain HP/LP (or BP/BR) signals by forming
//  sum and difference between input and allpass output (and dividing by two). This works for the 
//  1st order 1p1z case. Maybe it can be generalized. The advantage is that we have fewer degrees 
//  of freedom to tweak manually. We just set the poles, zeros are then obtained as their 
//  reciprocals. Alternatively, we may set the a-coeffs and the b-coeffs are obtained by reversal
//  (and perhaps normalization).
// -Try to express the above designs in terms of the allpass approach. That may or may not be 
//  possible...not sure yet...but if it is, it may be advantageous to have the allpass based 
//  version available because it may make for a good implementation structure, namely one that 
//  needs less coeffs. Only the a-coeffs need to be computed and stored. The b-coeffs are just the
//  reversed array. ...that actually already prescribes the algorithm for obtaining the allpass 
//  version, if it works at all: just take the a-array of lowpass or highpass (they are the same)
//  and reverse it. All filters: LP, HP and AP have the same poles, only the zeros differ.
// -Maybe try to do the same things in the s-domain. Or: take our z-domain designs and transfer 
//  them to the s-plane by the inverse bilinear transform. Then, when we have an analog prototype, 
//  we can use "prescribed Nyquist gain" techniques to move it back into the z-plane, thereby 
//  avoiding the cramping artifacts which the Constantinides transforms will surely produce.
