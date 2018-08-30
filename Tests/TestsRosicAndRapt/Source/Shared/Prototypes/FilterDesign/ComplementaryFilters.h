#pragma once


RAPT::rsFilterSpecificationBA<double> complementaryFilter(
  const RAPT::rsFilterSpecificationBA<double>& baSpec);

bool testSplitConditions(const RAPT::rsFilterSpecificationBA<double>& lpfBA);
// rename to testComplementarity


template<class T>
void plotMagnitudesBA(int numFreqs, T lowFreq, T highFreq, bool logFreqAxis, bool decibels,
  const std::vector<RAPT::rsFilterSpecificationBA<T>>& filterSpecs);
// make a simpler function

/** Given a filter-design specification in terms of polynomial coefficients, this function analyzes
whether or not the filter satisfies complementariness conditions. */
void analyzeComplementaryFilter(const RAPT::rsFilterSpecificationBA<double>& specBA);


RAPT::rsFilterSpecificationBA<double> complementaryLowpass2p3z();
// rename to complementaryFilter2p3z

// add complementaryFilter1p1z, 2p4z

