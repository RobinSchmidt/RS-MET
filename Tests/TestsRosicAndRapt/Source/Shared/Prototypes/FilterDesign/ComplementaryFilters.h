#pragma once


RAPT::rsFilterSpecificationBA<double> complementaryFilter(
  const RAPT::rsFilterSpecificationBA<double>& baSpec);

bool testSplitConditions(const RAPT::rsFilterSpecificationBA<double>& lpfBA);


template<class T>
void plotMagnitudesBA(int numFreqs, T lowFreq, T highFreq, bool logFreqAxis, bool decibels,
  const std::vector<RAPT::rsFilterSpecificationBA<T>>& filterSpecs);
// make a simpler function


RAPT::rsFilterSpecificationBA<double> splitterPrototype_2_3_new();
// rename to complementaryFilter2p3z

// add complementaryFilter1p1z, 2p4z

