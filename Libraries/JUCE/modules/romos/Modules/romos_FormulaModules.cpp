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
  *out = 0.0; // preliminary
}
CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_1(FormulaModule1In1Out);

}

}
