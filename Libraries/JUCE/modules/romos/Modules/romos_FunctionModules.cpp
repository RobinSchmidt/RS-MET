#include "romos_FunctionModules.h"

namespace romos
{

void ClipperModule::initialize()
{
  initInputPins(3, "In", "Min", "Max");
  initOutputPins(1, "Out");
}
INLINE void ClipperModule::process(Module *module, double *in1, double *in2, double *in3, double *out, int voiceIndex)
{
  *out = clip(*in1, *in2, *in3);
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_3(ClipperModule);


void SinCosModule::initialize()
{
  initInputPins(1, "In");
  initOutputPins(2, "Sin", "Cos");
}
INLINE void SinCosModule::process(Module *module, double *in, double *out, int voiceIndex)
{
  rosic::sinCos(2.0*PI*(*in), out, out+1);
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_1(SinCosModule);


void TriSawModule::initialize()
{
  initInputPins(6, "In", "Asym", "AtBn", "AtSg", "DcBn", "DcSg");
  initOutputPins(1, "");
}
INLINE void TriSawModule::process(Module *module, double *In, double *Asym, double *AtBn, 
  double *AtSg, double *DcBn, double *DcSg, double *out, int voiceIndex)
{
  double h  = 0.5 * (*Asym+1);
  double a0 = -1;
  double a1 = 2 / h;
  double b0 = (1+h)/(1-h);
  double b1 = -1 - b0;
  double p  = fmod(*In, 1);
  typedef RAPT::rsTriSawOscillator<double> TrSw;
  if(p < h)
    *out = TrSw::shape(a0+a1*p,  *AtBn, -0.5 * (*AtSg));    // upward section
  else 
    *out = TrSw::shape(b0+b1*p, -(*DcBn), -0.5 * (*DcSg));  // downward section
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_6(TriSawModule);



}
