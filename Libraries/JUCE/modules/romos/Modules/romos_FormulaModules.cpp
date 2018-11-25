namespace romos
{
    
//-------------------------------------------------------------------------------------------------

void FormulaModule1In1Out::initialize()
{
  initInputPins({ "x" });
  initOutputPins({ "y" });
  setFormula("y=x");
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

std::map<std::string, std::string> FormulaModule1In1Out::getState()
{
  std::map<std::string, std::string> state;
  state.emplace("Formula", formula);
  return state;
}

// move to RAPT:
template<class Key, class Value>
bool rsContains(const std::map<Key, Value>& map, const Key& key)
{
  auto iterator = map.find(key);
  if( iterator == map.end() )
    return false;
  return true;
}

void FormulaModule1In1Out::setState(const std::map<std::string, std::string>& state)
{
  //std::string tmp = state.find("Formula");
  //std::string tmp = state[state.find("Formula")];

  if(rsContains(state, std::string("Formula"))) 
  {
    std::string tmp = state.at(std::string("Formula"));
    RAPT::rsAssert(isFormulaValid(tmp));
    setFormula(tmp);
    // hmm...setFormula sets the formula only when it is valid - maybe we should set it
    // when it is invalid anyway. that would be more friendly behavior, if the user enters invalid
    // formulas in an xml state. It would be weird, if the module would keep the old formula 
    // without notifying the user - maybe the best thing would be to pop up an error message
  }
  else
  {
    // return false
  }


  //auto iterator = state.find("Formula");
  //if( iterator == state.end() )
  //  return;  // maybe we should return false, i.e. make that function bool

  //std::string tmp = state[iterator];
  //RAPT::rsAssert(isFormulaValid(tmp));
  //setFormula(tmp);
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

bool FormulaModule1In1Out::isFormulaValid(const std::string& formulaToTest)
{
  return trialEvaluator.setExpressionString(formulaToTest.c_str());
}

bool FormulaModule1In1Out::setFormula(const std::string& newFormula)
{
  if(isFormulaValid(newFormula)) {
    formula = newFormula;
    updateEvaluatorFormulas();
    return true;
  }
  return false;
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
    inVariables[i] = evaluators[i]->getVariableAddress("x");
}

CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_1(FormulaModule1In1Out);

/* Maybe, with a little trick, the formula module can be made to have memory:
old = 0;  // "declaration" of memory variable "old"
out = old;
old = In;
out;
...try it...If it works, we can build modules with memory (for example filters) from the 
formula module. ..but no - that doesn't work because the "declaration" would reset it in each
call, so out would alsway be zero. We somehow need to be able to declare variables without
assigning them to a value - check the ExprEval doc, if we can create variables. */


}