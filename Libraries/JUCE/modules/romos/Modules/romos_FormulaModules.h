#ifndef romos_FormulaModules_h
#define romos_FormulaModules_h

namespace romos
{

//-------------------------------------------------------------------------------------------------

/** A module that a applies a user defined formula y = f(x) with one input and one output. */
class FormulaModule1In1Out : public AtomicModule
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
    fullName     = "Formula1In1Out";
    description  = "A custom formula with one input and one output";
    category     = "Functions";
    createModule =  []()->Module* { return new FormulaModule1In1Out; };
    //hasHeader = false;
  }
};

//-------------------------------------------------------------------------------------------------
// a formula module with multiple ins and outs:


}

#endif 
