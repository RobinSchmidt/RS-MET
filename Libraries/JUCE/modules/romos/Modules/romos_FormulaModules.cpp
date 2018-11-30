namespace romos
{

//-------------------------------------------------------------------------------------------------

void FormulaModule_1_1::initialize()
{
  initInputPins({ "x" });
  initOutputPins({ "y" });
  setFormula("y=x");
}

INLINE void FormulaModule_1_1::process(Module *module, double *in, double *out, int voiceIndex)
{
  FormulaModule_1_1 *evaluatorModule = static_cast<FormulaModule_1_1*> (module);
  if(evaluatorModule->inVariables[voiceIndex] == nullptr) { 
    *out = 0; 
    return;  }

  rosic::ExpressionEvaluator* evaluator = evaluatorModule->evaluators[voiceIndex];
  //evaluator->assignVariable("x", *in); 
  *(evaluatorModule->inVariables[voiceIndex]) = *in;  // more efficient than assignVariable
  *out = evaluator->evaluateExpression();
}

void FormulaModule_1_1::resetVoiceState(int voiceIndex)
{
  AtomicModule::resetVoiceState(voiceIndex);
  // ...more to do?
}

std::map<std::string, std::string> FormulaModule_1_1::getState() const
{
  std::map<std::string, std::string> state = Module::getState();
  state.emplace("Formula", formula);
  return state;
}

bool FormulaModule_1_1::setState(const std::map<std::string, std::string>& state)
{
  bool result = Module::setState(state);
  if(RAPT::rsContains(state, std::string("Formula"))) 
  {
    std::string tmp = state.at(std::string("Formula"));
    RAPT::rsAssert(isFormulaValid(tmp));
    if(isFormulaValid(tmp)) {
      setFormula(tmp);
      return result;
    }
    else
      return false;
    // hmm...setFormula sets the formula only when it is valid - maybe we should set it
    // when it is invalid anyway. that would be more friendly behavior, if the user enters invalid
    // formulas in an xml state. It would be weird, if the module would keep the old formula 
    // without notifying the user - maybe the best thing would be to pop up an error message
  }
  else
    return false;
}

void FormulaModule_1_1::allocateMemory()
{
  AtomicModule::allocateMemory();
  inVariables.resize(getNumVoices());
  evaluators.resize(getNumVoices());
  for(int i = 0; i < getNumVoices(); i++)
    evaluators[i] = new rosic::ExpressionEvaluator;
  updateInputVariables();
}

void FormulaModule_1_1::freeMemory()
{
  AtomicModule::freeMemory();
  for(int i = 0; i < getNumVoices(); i++)
    delete evaluators[i];
  evaluators.clear();
  inVariables.clear();
}

bool FormulaModule_1_1::isFormulaValid(const std::string& formulaToTest)
{
  return trialEvaluator.setExpressionString(formulaToTest.c_str());
}

bool FormulaModule_1_1::setFormula(const std::string& newFormula)
{
  if(isFormulaValid(newFormula)) {
    formula = newFormula;
    updateEvaluatorFormulas();
    return true;
  }
  return false;
}

void FormulaModule_1_1::updateEvaluatorFormulas()
{
  for(int i = 0; i < evaluators.size(); i++)
  {
    evaluators[i]->initVariableList();  // experimental
    evaluators[i]->setExpressionString(formula.c_str());
  }
  updateInputVariables();
}

void FormulaModule_1_1::updateInputVariables()
{
  RAPT::rsAssert(inVariables.size() == evaluators.size());
  for(int i = 0; i < evaluators.size(); i++)
    inVariables[i] = evaluators[i]->getVariableAddress("x");
}

CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_1(FormulaModule_1_1);

/* Maybe, with a little trick, the formula module can be made to have memory:
old = 0;  // "declaration" of memory variable "old"
out = old;
old = In;
out;
...try it...If it works, we can build modules with memory (for example filters) from the 
formula module. ..but no - that doesn't work because the "declaration" would reset it in each
call, so out would alsway be zero. We somehow need to be able to declare variables without
assigning them to a value - check the ExprEval doc, if we can create variables. Maybe we can
introdcue memory variables, by just doing some assignments *after* all the outputs have been 
calculated. */

//-------------------------------------------------------------------------------------------------

double FormulaModule_N_1::dummyInput = 0.0;

void FormulaModule_N_1::initialize()
{
  initInputPins({ "x" });
  initOutputPins({ "y" });
  setFormula("y=x");
  setInputVariables("x"); // does this make the call to initInputPins superfluous? check this out
}

INLINE void FormulaModule_N_1::process(Module *module, double *in, double *out, int voiceIndex)
{
  FormulaModule_N_1 *formulaModule = static_cast<FormulaModule_N_1*> (module);
  rosic::ExpressionEvaluator* evaluator = formulaModule->evaluators[voiceIndex];

  for(unsigned int i = 0; i < formulaModule->numInputs; i++) 
    *(formulaModule->inVariablesN[voiceIndex][i]) = in[i]; // inject inputs
  // what, if the expression contains variables that don't appear in our inputVariablesN? maybe, we
  // should set them to zero

  *out = evaluator->evaluateExpression();

  // in a multi-output module later, we need to just call evaluator->evaluateExpression(); 
  // and then collect the outputs outputs here - the subclass can actually call the baseclass 
  // method and then collect
}

bool FormulaModule_N_1::setFormula(const std::string& newFormula)
{
  bool result = FormulaModule_1_1::setFormula(newFormula);
  updateInputVariables(); 
  return result;
}

std::map<std::string, std::string> FormulaModule_N_1::getState() const
{
  std::map<std::string, std::string> state = FormulaModule_1_1::getState();
  state.emplace("Inputs", inputVariableString);
  return state;
}

bool FormulaModule_N_1::setState(const std::map<std::string, std::string>& state)
{
  bool result = FormulaModule_1_1::setState(state);
  if(RAPT::rsContains(state, std::string("Inputs"))) {
    std::string inputStr = state.at(std::string("Inputs"));
    result &= setInputVariables(inputStr);
  }
  else {
    setInputVariables("x");
    result = false;
  }
  RAPT::rsAssert(result);
  return result;
}


// move to rosic:

// use function from RAPT::rsArray - but needs adaption of parameter types (constness)
inline int findIndexOf(const char* buffer, char elementToFind, int length)
{
  for(int i = 0; i < length; i++) {
    if( buffer[i] == elementToFind )
      return i;
  }
  return -1;
}
// http://www.cplusplus.com/reference/string/string/substr/
std::vector<std::string> tokenize(const std::string& str, const char splitChar)
{
  std::vector<std::string> result;
  int start = 0;
  while(start < str.size()) {

    int delta = findIndexOf(&str[start], splitChar, (int)str.size()-start);
    // use http://www.cplusplus.com/reference/string/string/find/


    if(delta == -1)
      break;
    std::string token = str.substr(start, delta);
    result.push_back(token);
    start += delta+1; // +1 for the splitChar itself
  }
  result.push_back(str.substr(start, str.size()-start)); // add tail
  return result;
}
void removeChar(std::string& str, const char chr)
{
  std::string::iterator end_pos = std::remove(str.begin(), str.end(), chr);
  str.erase(end_pos, str.end());
  // from https://stackoverflow.com/questions/83439/remove-spaces-from-stdstring-in-c
}


bool FormulaModule_N_1::setInputVariables(const std::string& newInputs)
{
  inputVariableString = newInputs;
  std::vector<std::string> strArr = tokenize(newInputs, ','); // collect variable names in array

  for(size_t i = 0; i < strArr.size(); i++) {
    removeChar(strArr[i], ' '); // remove whitespaces
    // todo: remove anything after a colon to allow comments - but maybe keep the strings after
    // the colon in a separate array to be used for the infoline whe user hovers over pins
  }

  setInputVariables(strArr);
  return true; // preliminary
}

void FormulaModule_N_1::setInputVariables(const std::vector<std::string>& newInVars)
{
  // todo:
  // -check, if variable names are valid
  // -have a boolean return value to report error, when names are invalid

  size_t oldSize = audioInputNames.size(); // == inputPins.size()
  size_t newSize = newInVars.size();
  audioInputNames.resize(newSize);
  inputPins.resize(newSize);
  for(size_t i = 0; i < newSize; i++)
    audioInputNames[i] = newInVars[i];
  numInputs = (int) inputPins.size();  // try to get rid of numInputs - i think, it's redundant
  updateInputVariables();
}

void FormulaModule_N_1::allocateMemory()
{
  AtomicModule::allocateMemory();
  inVariablesN.resize(getNumVoices());     // sole difference to baseclass - maybe refactor
  evaluators.resize(getNumVoices());
  for(int i = 0; i < getNumVoices(); i++)
    evaluators[i] = new rosic::ExpressionEvaluator;
  updateInputVariables();
}

void FormulaModule_N_1::freeMemory()
{
  AtomicModule::freeMemory();
  for(int i = 0; i < getNumVoices(); i++)
    delete evaluators[i];
  evaluators.clear();
  inVariablesN.clear(); // sole difference to baseclass - maybe refactor
  // maybe we should set all input pointers to point to the dummyInput - maybe this should be
  // done in the baseclass as well, like
  // invalidateInputVariablePointers();
}

void FormulaModule_N_1::updateInputVariables()
{
  RAPT::rsAssert(inVariablesN.size() == evaluators.size());
  for(size_t i = 0; i < evaluators.size(); i++) {        // loop over the voices
    inVariablesN[i].resize(audioInputNames.size());
    for(size_t j = 0; j < inVariablesN[i].size(); j++) { // loop over the input variables
      std::string varName = audioInputNames[j].asStdString();
      double* varPtr = evaluators[i]->getVariableAddress(varName.c_str());
      if(varPtr != nullptr)
        inVariablesN[i][j] = varPtr;
      else
        inVariablesN[i][j] = &dummyInput;  // points to a zero valued memory location
    }
  }
  int dummy = 0;
}

CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_N(FormulaModule_N_1);



}