#pragma once

/** Given a filter-design specification in terms of polynomial coefficients, this function analyzes
whether or not the filter satisfies complementariness conditions. These conditions mean that if we 
obtain a "complementary" signal by subtracting the filter's output from the input, the resulting 
frequency response of the difference-filter is a mirror image of the frequency response of the 
given filter. */
bool isComplementary(const RAPT::rsFilterSpecificationBA<double>& specBA);
// todo: maybe rename the condition...or maybe not

/** Makes some plots, etc. */
bool analyzeComplementaryFilter(const RAPT::rsFilterSpecificationBA<double>& specBA);


// the actual prototype designs that satisfy the mirror-imag e conditions:
RAPT::rsFilterSpecificationBA<double> complementaryLowpass1p1z();
RAPT::rsFilterSpecificationBA<double> complementaryLowpass2p2z();
RAPT::rsFilterSpecificationBA<double> complementaryLowpass2p3z();

RAPT::rsFilterSpecificationBA<double> complementaryLowpass2p4z();


//RAPT::rsFilterSpecificationBA<double> complementaryLowpass3p3z();

// add complementaryFilter2p4z, 4p4z maybe try to place 3 zeros at z = -1

// todo: obtain frequency warped versions via Constantinides formulas (the functions above design
// prototypes with their cutoff frequency being a quarter of the sample rate (halfband filters)
