//#include "romos_WorkArea.h"
//using namespace romos;

double  WorkArea::tmpInFrame [maxNumPins];
double  WorkArea::tmpOutFrame[maxNumPins];

double  WorkArea::tmpInFramePoly [maxNumVoices*maxNumPins];
double  WorkArea::tmpOutFramePoly[maxNumVoices*maxNumPins];

double* WorkArea::inVoiceFramePointer[maxNumPins];

//int     WorkArea::inFrameStrides[maxNumPins];


/*
WorkArea romos::workArea;  // definition of the global object

WorkArea::WorkArea()
{
  dummySourceModule = ModuleFactory::createModule(ModuleTypeRegistry::CONSTANT, "DummySourceModule"

}

WorkArea::~WorkArea()
{

}
*/

