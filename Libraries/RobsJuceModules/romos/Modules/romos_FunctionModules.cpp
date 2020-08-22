//#include "romos_FunctionModules.h"

//-------------------------------------------------------------------------------------------------

void ClipperModule::initialize()
{
  initInputPins({ "In", "Min", "Max" });
  initOutputPins({ "Out" });
  inputPins[1].setDefaultValue(-1.0); 
  inputPins[2].setDefaultValue(+1.0); 
}
INLINE void ClipperModule::process(Module *module, double *in1, double *in2, double *in3, double *out, int voiceIndex)
{
  *out = RAPT::rsClip(*in1, *in2, *in3);
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_3(ClipperModule);

//-------------------------------------------------------------------------------------------------

void SinCosModule::initialize()
{
  initInputPins({ "In" });
  initOutputPins({ "Sin", "Cos" });
}
INLINE void SinCosModule::process(Module *module, double *in, double *out, int voiceIndex)
{
  RAPT::rsSinCos(2.0*PI*(*in), out, out+1);
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_1(SinCosModule);

//-------------------------------------------------------------------------------------------------

void TriSawModule::initialize()
{
  initInputPins({ "In", "Asym", "AtBn", "AtSg", "DcBn", "DcSg" });
  initOutputPins({ "" });
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

//-------------------------------------------------------------------------------------------------

void SaturatorModule::initialize()
{
  //initInputPins(3, "In", "Width", "Center");
  //initOutputPins(1, "Out");
  initInputPins({ "In", "Width", "Center" });
  initOutputPins({ "Out" });
  inputPins[1].setDefaultValue(2); // Width is 2 by default
}
INLINE void SaturatorModule::process(Module *module, double *In, double *Width, double *Center,
  double *out, int voiceIndex)
{
  SaturatorModule* sat = static_cast<SaturatorModule*> (module);
  double scaleX =  2 / *Width;
  double shiftX = -1 - (scaleX * (*Center - 0.5 * *Width));
  double scaleY =  1 / scaleX; 
  double shiftY = -shiftX * scaleY;
  *out = shiftY + scaleY * sat->sigmoid(scaleX * *In + shiftX);

  //*out = sat->sigmoid(*In);  // preliminary
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_3(SaturatorModule);
