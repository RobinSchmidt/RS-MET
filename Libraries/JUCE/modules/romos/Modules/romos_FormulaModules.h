#ifndef romos_FormulaModules_h
#define romos_FormulaModules_h

namespace romos
{

//-------------------------------------------------------------------------------------------------

/** A module that a applies a user defined formula y = f(x) with one input and one output. */
class FormulaModule1In1Out : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_1(FormulaModule1In1Out);
};
class FormulaModule1In1OutTypeInfo : public ModuleTypeInfo
{
public:
  FormulaModule1In1OutTypeInfo() {
    shortName    = "Formula_1_1";
    fullName     = "Formula, 1 In, 1 Out";
    description  = "A custom formula with one input and one output";
    category     = "Functions";
    createModule =  []()->Module* { return new FormulaModule1In1Out; };
    //hasHeader = false;
  }
};

}

#endif 
