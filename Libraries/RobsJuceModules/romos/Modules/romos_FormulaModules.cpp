
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
  resetVariables(voiceIndex);
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

void FormulaModule_1_1::resetVariables(int voiceIndex)
{
  evaluators[voiceIndex]->resetVariables();
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
  for(size_t i = 0; i < evaluators.size(); i++)
  {
    evaluators[i]->initVariableList();  // without it, the evaluator may still remember old input
                                        // variables from previous formulas - maybe that can be
                                        // exploited for introducing memory variables later
    evaluators[i]->setExpressionString(formula.c_str());
  }
  updateInputVariables();
}

void FormulaModule_1_1::updateInputVariables()
{
  RAPT::rsAssert(inVariables.size() == evaluators.size());
  for(size_t i = 0; i < evaluators.size(); i++)
    inVariables[i] = evaluators[i]->getVariableAddress("x");
}

CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_1(FormulaModule_1_1);

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
  *out = evaluator->evaluateExpression();
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

  std::vector<std::pair<AudioConnection, std::string>> 
    inputConnections = getInputVariableConnections(); // remember connectivity

  size_t oldSize = audioInputNames.size(); // == inputPins.size()
  size_t newSize = newInVars.size();
  audioInputNames.resize(newSize);
  inputPins.resize(newSize);
  for(size_t i = 0; i < newSize; i++)
    audioInputNames[i] = newInVars[i];
  numInputs = (int) inputPins.size();  // try to get rid of numInputs - i think, it's redundant
  updateInputVariables();

  restoreInputVariableConnections(inputConnections); // restore connectivity
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
}

std::vector<std::pair<AudioConnection, std::string>> 
FormulaModule_N_1::getInputVariableConnections()
{
  std::vector<std::pair<AudioConnection, std::string>> pairs;
  std::vector<AudioConnection> cons = getIncomingAudioConnections();
  for(size_t i = 0; i < cons.size(); i++) {
    AudioConnection connection = cons[i];
    int inPinIndex = connection.getTargetInputIndex();
    std::string inPinName = audioInputNames[i].asStdString();
    pairs.push_back(std::pair<AudioConnection, std::string>(connection, inPinName));
  }
  return pairs;
}

void FormulaModule_N_1::restoreInputVariableConnections(
  const std::vector<std::pair<AudioConnection, std::string>>& connections)
{
  // disconnect all inputs and create all connections anew - that's easiest to implement, but
  // maybe it can optimized later
  disconnectAllInputPins();
  for(size_t i = 0; i < connections.size(); i++) {
    AudioConnection con = connections[i].first;
    std::string inPinName = connections[i].second;
    size_t inPinIndex = RAPT::rsFind(audioInputNames, rosic::rsString(inPinName));
    if(inPinIndex < audioInputNames.size())
      connectInputPinTo((int)inPinIndex, con.getSourceModule(), con.getSourceOutputIndex());
  }
  // maybe it should return whether or not the connections were actually modified
}

CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_N(FormulaModule_N_1);


//-------------------------------------------------------------------------------------------------

void FormulaModule_N_M::initialize()
{
  initInputPins({ "x" });
  initOutputPins({ "y" });
  setFormula("y=x");
  setInputVariables("x"); // does this make the call to initInputPins superfluous? check this out
  setOutputVariables("y");
}  

INLINE void FormulaModule_N_M::process(Module *module, double *in, double *out, int voiceIndex)
{
  FormulaModule_N_1::process(module, in, out, voiceIndex);
  FormulaModule_N_M *formulaModule = static_cast<FormulaModule_N_M*> (module);
  for(unsigned int i = 0; i < formulaModule->outFrameStride; i++) 
    out[i] = *(formulaModule->outVariablesM[voiceIndex][i]); // collect outputs
}

void FormulaModule_N_M::resetVoiceState(int voiceIndex)
{
  FormulaModule_N_1::resetVoiceState(voiceIndex);
  for(unsigned int i = 0; i < outFrameStride; i++) 
    *(outVariablesM[voiceIndex][i]) = 0.0; // reset outputs
}

bool FormulaModule_N_M::setFormula(const std::string& newFormula)
{
  bool result = FormulaModule_N_1::setFormula(newFormula);
  updateOutputVariables(); 
  return result;
}

std::map<std::string, std::string> FormulaModule_N_M::getState() const
{
  std::map<std::string, std::string> state = FormulaModule_N_1::getState();
  state.emplace("Outputs", outputVariableString);
  return state;
}

bool FormulaModule_N_M::setState(const std::map<std::string, std::string>& state)
{
  bool result = FormulaModule_N_1::setState(state);
  if(RAPT::rsContains(state, std::string("Outputs"))) {
    std::string outputStr = state.at(std::string("Outputs"));
    result &= setOutputVariables(outputStr);
  }
  else {
    setOutputVariables("y");
    result = false;
  }
  RAPT::rsAssert(result);
  return result;
}

void FormulaModule_N_M::allocateMemory()
{
  FormulaModule_N_1::allocateMemory();
  outVariablesM.resize(getNumVoices());
  updateOutputVariables();
}

void FormulaModule_N_M::freeMemory()
{
  FormulaModule_N_1::freeMemory();
  outVariablesM.clear();
}

bool FormulaModule_N_M::setOutputVariables(const std::string& newOutputs)
{
  outputVariableString = newOutputs;
  std::vector<std::string> strArr = tokenize(newOutputs, ',');
  for(size_t i = 0; i < strArr.size(); i++)
    removeChar(strArr[i], ' ');
  setOutputVariables(strArr);
  return true; // preliminary
}

void FormulaModule_N_M::setOutputVariables(const std::vector<std::string>& newOutVars)
{
  std::vector<std::pair<AudioConnection, std::string>>
    outputConnections = getOutputVariableConnections();

  size_t oldSize = audioOutputNames.size();
  size_t newSize = newOutVars.size();
  audioOutputNames.resize(newSize);
  for(size_t i = 0; i < newSize; i++)
    audioOutputNames[i] = newOutVars[i];

  if(newSize != outFrameStride) {
    outFrameStride = (int) newSize;
    allocateAudioOutputs(); // or maybe call allocateMemory?
  }
  updateOutputVariables();

  restoreOutputVariableConnections(outputConnections);
}

void FormulaModule_N_M::updateOutputVariables()
{
  RAPT::rsAssert(outVariablesM.size() == evaluators.size());
  for(size_t i = 0; i < evaluators.size(); i++) {         // loop over the voices
    outVariablesM[i].resize(audioOutputNames.size());
    for(size_t j = 0; j < outVariablesM[i].size(); j++) { // loop over the output variables
      std::string varName = audioOutputNames[j].asStdString();
      double* varPtr = evaluators[i]->getVariableAddress(varName.c_str());
      if(varPtr != nullptr)
        outVariablesM[i][j] = varPtr;
      else
        outVariablesM[i][j] = &dummyInput;  // points to a zero valued memory location
    }
  }
}

std::vector<std::pair<AudioConnection, std::string>> 
FormulaModule_N_M::getOutputVariableConnections()
{
  std::vector<std::pair<AudioConnection, std::string>> pairs;
  std::vector<AudioConnection> cons = getOutgoingAudioConnections();
  for(size_t i = 0; i < cons.size(); i++) {
    AudioConnection connection = cons[i];
    int outPinIndex = connection.getSourceOutputIndex();
    std::string outPinName = audioOutputNames[i].asStdString();
    pairs.push_back(std::pair<AudioConnection, std::string>(connection, outPinName));
  }
  return pairs;
}

void FormulaModule_N_M::restoreOutputVariableConnections(
  const std::vector<std::pair<AudioConnection, std::string>>& connections)
{
  disconnectAllOutputPins();
  for(size_t i = 0; i < connections.size(); i++) {
    AudioConnection con = connections[i].first;
    std::string outPinName = connections[i].second;
    size_t outPinIndex = RAPT::rsFind(audioOutputNames, rosic::rsString(outPinName));
    if(outPinIndex < audioOutputNames.size())
      con.getTargetModule()->connectInputPinTo(con.getTargetInputIndex(), this, (int)outPinIndex);
  }
}

CREATE_AND_ASSIGN_PROCESSING_FUNCTIONS_N(FormulaModule_N_M);


/*
Ideas:
Make an ordinary differential equation solver based on an expression evaluator. The ODE is 
expressed as y' = f(x,y) where y, y' are vectors and f is a vector valued function (maybe we can 
absorb the independent variable x into the vector y as well - we'll see)
-inputs: parameters to the ODE, maybe time in case of non-autonomous systems (we need an absolute 
 time module then, that outputs the time passed since note-start in samples or seconds)
-outputs: elements of the vector y (maybe not all of them)
-the formula should then give expressions for all the element functions of f
-we need to establish a convention for the notation for this - maybe if the output variables are
 named x,y,z, the derivatives could be named xp, yp, zp (p for prime) or dx, dy, dz as shothand
 for dx/dt, dy/dt, dz/dt, or xd, yd, zd (d for "dot" or "derivative")
-on the GUI teh user selects the algorithm (Euler, Runge-Kutta, etc.) and maybe the set size...
 but maybe the stepsize should be an input to make it dependent on sample-rate, frequency, etc.
-maybe the solver needs an internal stepsize independent from the sampl-rate/frequency and a 
 resampler
-for the implementation either use multiple inheritance (from ExpressionEvaluator and 
 std::function) for the definition of the function f and have an ODESolver object for each 
 expression evaluator that points to the std::function part or have a subclass of std::function
 that maintains its expression evaluator as member...maybe that's better because it's more 
 generally useful
*/