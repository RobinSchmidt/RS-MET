#ifndef RAPT_FILTEREXPERIMENTS_H
#define RAPT_FILTEREXPERIMENTS_H

void bandSplittingTwoWay();
void bandSplittingMultiWay();
void bandSplittingTreeAlgo();
void bandSplitFreqResponses();
void complementaryFiltersIIR();
void firstOrderFilters();
void ladderResonanceManipulation();
// insert movingAverage here

void nonUniformMovingAverage();
void nonUniformOnePole1();
void nonUniformOnePole2();
void nonUniformBiquad();

void smoothingFilterOrders();
void smoothingFilterTransitionTimes();

void prototypeDesign();    // still using the old code
void poleZeroPrototype();  // new implementation


#endif