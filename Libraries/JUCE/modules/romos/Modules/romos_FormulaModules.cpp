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
  //*out = *in;
  //return; // for debugging memelak


  FormulaModule1In1Out *evaluatorModule = static_cast<FormulaModule1In1Out*> (module);
  rosic::ExpressionEvaluator* evaluator = evaluatorModule->evaluators[voiceIndex];

  evaluator->assignVariable("In", *in); 
  // maybe use getVariableAddress - the doc says, it's more efficient than assignVariable

  //// for debug:
  //double tmp1, tmp2;
  //tmp1 = tanh(2 * (*in) * (*in));
  //tmp2 = evaluator->evaluateExpression(); // this should become our output
  ////RAPT::rsAssert(tmp1 == tmp2);

  *out = evaluator->evaluateExpression();
}

void FormulaModule1In1Out::resetVoiceState(int voiceIndex)
{
  AtomicModule::resetVoiceState(voiceIndex);
  // ...more to do...
}

void FormulaModule1In1Out::allocateMemory()
{
  AtomicModule::allocateMemory();
  evaluators.resize(getNumVoices());
  for(int i = 0; i < getNumVoices(); i++)
    evaluators[i] = new rosic::ExpressionEvaluator;

  // create buffers for formula variables
}

void FormulaModule1In1Out::freeMemory()
{
  AtomicModule::freeMemory();
  for(int i = 0; i < getNumVoices(); i++)
    delete evaluators[i];
  evaluators.clear();

  // free buffers for formula variables
}

void FormulaModule1In1Out::setFormula(const std::string& newFormula)
{
  formula = newFormula;
  updateEvaluatorFormulas();
}

void FormulaModule1In1Out::updateEvaluatorFormulas()
{
  for(int i = 0; i < evaluators.size(); i++)
    evaluators[i]->setExpressionString(formula.c_str());
}


CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_1(FormulaModule1In1Out);


}