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

  virtual void setFormula(const std::string& newFormula);
  virtual void resetVoiceState(int voiceIndex);

protected:
  virtual void allocateMemory();
  virtual void freeMemory();
  virtual void updateEvaluatorFormulas();

  std::vector<rosic::ExpressionEvaluator*> evaluators;
  std::string formula;

  double *variables; // not yet used and may not be needed - we'll see
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

}

#endif 
