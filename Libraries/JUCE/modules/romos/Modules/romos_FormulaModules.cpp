namespace romos
{
    
//-------------------------------------------------------------------------------------------------

void FormulaModule1In1Out::initialize()
{
  initInputPins({ "In" });
  initOutputPins({ "Out" });
}
INLINE void FormulaModule1In1Out::process(Module *module, double *in, double *out, int voiceIndex)
{
  //*out = 0.0; 

  *out = tanh(2 * (*in) * (*in)); // preliminary: y = tanh(2*x^2)
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_1(FormulaModule1In1Out);


}