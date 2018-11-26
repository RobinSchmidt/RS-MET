#ifndef romos_FormulaModules_h
#define romos_FormulaModules_h

namespace romos
{

//-------------------------------------------------------------------------------------------------

/** A module that a applies a user defined formula y = f(x) with one input and one output. */
class FormulaModule1In1Out : public AtomicModule // rename to FormulaModule_1_1
{
  CREATE_COMMON_DECLARATIONS_1(FormulaModule1In1Out);

public:

  virtual bool isFormulaValid(const std::string& formulaToTest);
  virtual bool setFormula(const std::string& newFormula);
  virtual void resetVoiceState(int voiceIndex);
  virtual std::string getFormula() { return formula; }

  virtual std::map<std::string, std::string> getState() const override;
  virtual bool setState(const std::map<std::string, std::string>& state) override;

protected:

  virtual void allocateMemory();
  virtual void freeMemory();
  virtual void updateEvaluatorFormulas();
  virtual void updateInputVariables();

  rosic::ExpressionEvaluator trialEvaluator; // to validate formulas maybe rename to formulaValidator
  std::vector<rosic::ExpressionEvaluator*> evaluators; // the evaluators used for actual dsp
  std::string formula;
  std::vector<double*> inVariables; // pointers to the input variables in the expression evaluator object
};
class FormulaModule1In1OutTypeInfo : public ModuleTypeInfo
{
public:
  FormulaModule1In1OutTypeInfo() {
    shortName    = "Formula";
    fullName     = "Formula1In1Out"; // maybe use _1_1
    description  = "A custom formula with one input and one output";
    category     = "Functions";
    createModule =  []()->Module* { return new FormulaModule1In1Out; };
  }
};

//-------------------------------------------------------------------------------------------------

/** A formula module with multiple inputs and one output. */
class FormulaModule_N_1 : public FormulaModule1In1Out
{
  CREATE_COMMON_DECLARATIONS_N(FormulaModule_N_1);

public:

  // overrides:
  virtual bool setFormula(const std::string& newFormula);
  virtual void resetVoiceState(int voiceIndex);
  virtual std::map<std::string, std::string> getState() const override;
  virtual bool setState(const std::map<std::string, std::string>& state) override;

protected:

  // overrides:
  virtual void allocateMemory();
  virtual void freeMemory();
  virtual void updateInputVariables();

  std::vector<double**> inVariables; // 1st index: voice, 2nd index: variable, masks inherited
                                     // variable of same name
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




}

#endif 
