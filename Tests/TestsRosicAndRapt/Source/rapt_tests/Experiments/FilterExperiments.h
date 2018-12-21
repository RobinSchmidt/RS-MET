#ifndef RAPT_FILTEREXPERIMENTS_H
#define RAPT_FILTEREXPERIMENTS_H

#include "../../Shared/Shared.h"

void bandSplittingTwoWay();
void bandSplittingMultiWay();
void bandSplittingTreeAlgo();
void bandSplitFreqResponses();
void complementaryFiltersIIR();
void firstOrderFilters();
void ladderResonanceManipulation();
// insert movingAverage here
void nonUniformMovingAverage();
void smoothingFilterOrders();
void smoothingFilterTransitionTimes();

void prototypeDesign();    // still using the old code
void poleZeroPrototype();  // new implementation


#endif