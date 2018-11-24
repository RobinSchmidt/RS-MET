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
  FormulaModule1In1Out *evaluatorModule = static_cast<FormulaModule1In1Out*> (module);
  if(evaluatorModule->inVariables[voiceIndex] == nullptr) { 
    *out = 0; 
    return; 
  }

  rosic::ExpressionEvaluator* evaluator = evaluatorModule->evaluators[voiceIndex];
  //evaluator->assignVariable("In", *in); 
  *(evaluatorModule->inVariables[voiceIndex]) = *in;  // probably more efficient than assignVariable
  *out = evaluator->evaluateExpression();
}

void FormulaModule1In1Out::resetVoiceState(int voiceIndex)
{
  AtomicModule::resetVoiceState(voiceIndex);
  // ...more to do?
}

void FormulaModule1In1Out::allocateMemory()
{
  AtomicModule::allocateMemory();
  inVariables.resize(getNumVoices());
  evaluators.resize(getNumVoices());
  for(int i = 0; i < getNumVoices(); i++)
    evaluators[i] = new rosic::ExpressionEvaluator;
  updateInputVariables();
}

void FormulaModule1In1Out::freeMemory()
{
  AtomicModule::freeMemory();
  for(int i = 0; i < getNumVoices(); i++)
    delete evaluators[i];
  evaluators.clear();
  inVariables.clear();
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
  updateInputVariables();
}

void FormulaModule1In1Out::updateInputVariables()
{
  RAPT::rsAssert(inVariables.size() == evaluators.size());
  for(int i = 0; i < evaluators.size(); i++)
    inVariables[i] = evaluators[i]->getVariableAddress("In");
}


CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_1(FormulaModule1In1Out);


}