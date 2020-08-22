#ifndef romos_WorkArea_h
#define romos_WorkArea_h

#define MAX(x,y) (((x)>(y))?(x):(y))  // move somewhere else - or replace by std::max or rsMax


/** This class is used internally inside the DSP code to store some temporary data. We don't want 
to allocate memory dynamically in the DSP code, that's why we use a global work area for this. */

//class ConstantModule;

struct WorkArea
{
  WorkArea();
  ~WorkArea();

  static const int maxNumVoices = 32;

  static const int maxNumPins   = MAX(64, ProcessingStatus::maxBufferSize);
  // necessary because in ContainerModule's processBlockFrameWise functions, the voice-index 
  // temporarily takes over the role of the frame index - there is only one frame for which the 
  // process function is called, but we may need a lot of pins, so we do this role-changing to 
  // avoid unnecessary memory overheads ...comment this better

  static double tmpInFrame[maxNumPins];
  static double tmpOutFrame[maxNumPins];

  static double tmpInFramePoly[maxNumVoices*maxNumPins];
  static double tmpOutFramePoly[maxNumVoices*maxNumPins];

  //static double tmpInFramePoly [maxNumVoices][maxNumPins]; 
  //static double tmpOutFramePoly[maxNumVoices][maxNumPins]; 
  // we should use flat arrays and maxNumPins should be >= processingStatus.bufferSize - we want 
  // the distance between successive voices to be equal to the allocated blocksize

  static double *inVoiceFramePointer[maxNumPins];

  //static int    inFrameStrides[maxNumPins];  // maybe rename to inFrameSizes - obsolete

  // maybe include facilities to let the maxNumInputPins grow, if necessary

  //romos::ConstantModule *dummySourceModule;
};

//extern WorkArea workArea;



#endif
