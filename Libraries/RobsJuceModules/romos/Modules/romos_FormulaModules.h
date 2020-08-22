#ifndef romos_FormulaModules_h
#define romos_FormulaModules_h

//-------------------------------------------------------------------------------------------------

/** A module that a applies a user defined formula y = f(x) with one input and one output. */

class FormulaModule_1_1 : public AtomicModule // rename to FormulaModule_1_1
{
  CREATE_COMMON_DECLARATIONS_1(FormulaModule_1_1);

public:

  virtual void resetVoiceState(int voiceIndex) override;
  virtual bool isFormulaValid(const std::string& formulaToTest);
  virtual bool setFormula(const std::string& newFormula);
  virtual std::string getFormula() { return formula; }

  virtual std::map<std::string, std::string> getState() const override;
  virtual bool setState(const std::map<std::string, std::string>& state) override;


protected:

  virtual void resetVariables(int voiceIndex);
  //void setInitialVariableValues()
  virtual void allocateMemory();
  virtual void freeMemory();
  virtual void updateEvaluatorFormulas();
  virtual void updateInputVariables();

  rosic::ExpressionEvaluator trialEvaluator; // to validate formulas maybe rename to formulaValidator
  std::vector<rosic::ExpressionEvaluator*> evaluators; // the evaluators used for actual dsp
  std::string formula;
  std::vector<double*> inVariables; // pointers to the input variables in the expression evaluator object
};
class FormulaModule_1_1TypeInfo : public ModuleTypeInfo
{
public:
  FormulaModule_1_1TypeInfo() {
    shortName    = "Formula";
    fullName     = "Formula_1_1";
    description  = "A custom formula with one input and one output";
    category     = "Functions";
    createModule =  []()->Module* { return new FormulaModule_1_1; };
  }
};

//-------------------------------------------------------------------------------------------------

/** A formula module with multiple inputs and one output. */

class FormulaModule_N_1 : public FormulaModule_1_1
{
  CREATE_COMMON_DECLARATIONS_N(FormulaModule_N_1);

public:

  // overrides:
  virtual bool setFormula(const std::string& newFormula) override;
  //virtual void resetVoiceState(int voiceIndex) override;
  virtual std::map<std::string, std::string> getState() const override;
  virtual bool setState(const std::map<std::string, std::string>& state) override;


  /** Sets up the input variables from a string that encodes what input variables the fomula 
  should have using a simple syntax. If the formula is, for example, y = a*x + b, we would have the 
  3 input variables a, b, x. These could be set up, by using the string: "a; b; x" - i.e. variable
  names a re separated by a semicolon

  ...todo: allow comments/descriptions after a colon, like  "x: input value; a: slope; b: offset",
  and state clearly how white-space is handled, what valid variable names are, etc. (i.e. start 
  with letter (or underscore or maybe not), contain letters and numbers, etc.) */
  virtual bool setInputVariables(const std::string& newInputs);

  /** Returns the string that defines the input variables. */
  virtual std::string getInputVariables() { return inputVariableString; }


protected:

  virtual void setInputVariables(const std::vector<std::string>& newInputVariables);

  /** Returns a array of pairs of an incoming audio-connection and an associated input variable 
  name. This is needed to restore the desired connections after a call to updateInputVariables. */
  virtual std::vector<std::pair<AudioConnection, std::string>>
    getInputVariableConnections();

  virtual void restoreInputVariableConnections(
    const std::vector<std::pair<AudioConnection, std::string>>& connections);


  // overrides:
  virtual void allocateMemory() override;
  virtual void freeMemory() override;
  virtual void updateInputVariables() override;

  std::string inputVariableString;   // the string that specifies the input variables


  //std::vector<double**> inVariablesN;
  std::vector<std::vector<double*>> inVariablesN;
   // 1st index: voice, 2nd index: variable
   // makes inherited inVariables variable obsolete (it's not used in this subclass)

  static double dummyInput; // = 0, that's where we point to, when we don't find a variable in the
                            // expression evaluator objects - maybe move to Module baseclass

  friend bool testFormulaModules(); // unit test, needs access to protected members
};
class FormulaModule_N_1TypeInfo : public ModuleTypeInfo
{
public:
  FormulaModule_N_1TypeInfo() {
    shortName    = "Formula";
    fullName     = "Formula_N_1";
    description  = "A custom formula with multiple inputs and one output";
    category     = "Functions";
    createModule =  []()->Module* { return new FormulaModule_N_1; };
  }
};

//-------------------------------------------------------------------------------------------------

/** A formula module with multiple inputs and multiple outputs. */

class FormulaModule_N_M : public FormulaModule_N_1
{
  CREATE_COMMON_DECLARATIONS_N(FormulaModule_N_M);


public:

  virtual void resetVoiceState(int voiceIndex) override;
  virtual bool setFormula(const std::string& newFormula) override;
  virtual std::map<std::string, std::string> getState() const override;
  virtual bool setState(const std::map<std::string, std::string>& state) override;


  virtual bool setOutputVariables(const std::string& newOutputs);

  /** Returns the string that defines the output variables. */
  virtual std::string getOutputVariables() { return outputVariableString; }


protected:

  virtual void setOutputVariables(const std::vector<std::string>& newOutputVariables);

  virtual std::vector<std::pair<AudioConnection, std::string>>
    getOutputVariableConnections();

  virtual void restoreOutputVariableConnections(
    const std::vector<std::pair<AudioConnection, std::string>>& connections);


  virtual void allocateMemory() override;
  virtual void freeMemory() override;
  virtual void updateOutputVariables();

  std::string outputVariableString;

  std::vector<std::vector<double*>> outVariablesM;


  friend bool testFormulaModules(); // unit test, needs access to protected members
};
class FormulaModule_N_MTypeInfo : public ModuleTypeInfo
{
public:
  FormulaModule_N_MTypeInfo() {
    shortName    = "Formula";
    fullName     = "Formula"; // "Formula_N_M";
    description  = "A custom formula with multiple inputs and multiple outputs";
    category     = "Functions";
    createModule =  []()->Module* { return new FormulaModule_N_M; };
  }
};

#endif 
