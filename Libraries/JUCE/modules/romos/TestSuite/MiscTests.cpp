

bool romos::testFormulaModules()
{
  bool result = true;

  typedef romos::Module MDL;
  ContainerModule* cm = (ContainerModule*)moduleFactory.createModule("Container",
    "FormulaTest", 0, 0, true);
  // actually, we don't really need this container...yet...but later we may want to make tests for
  // the updating of the connectivity when the input variables change...

  /*
  MDL* constant3    = cm->addChildModule("Constant", "3", 2, 4, false, false);
  MDL* constant5    = cm->addChildModule("Constant", "5", 2, 6, false, false);
  MDL* constant7    = cm->addChildModule("Constant", "7", 2, 8, false, false);
  MDL* audioOutput1 = cm->addChildModule("AudioOutput", "Out1", 20, 6, false, false);
  */

  FormulaModule_N_1* formula_N_1 = (FormulaModule_N_1*)
    cm->addChildModule("Formula_N_1", "Formula_N_1", 13, 6, false, false);

  cm->sortChildModuleArray();

  double ins[4] = { 2, 3, 5, 7};
  double outs[4];

  // 3 inputs, all used in formula:
  formula_N_1->setFormula("y = a*x + b");
  formula_N_1->setInputVariables("x,a,b");
  FormulaModule_N_1::process(formula_N_1, ins, outs, 0);
  result &= outs[0] == 3*2 + 5;

  // 3 inputs, 3rd not used:
  formula_N_1->setFormula("y = a*x");
  formula_N_1->setInputVariables("x,a,b");
  FormulaModule_N_1::process(formula_N_1, ins, outs, 0);
  result &= outs[0] == 3*2;

  // for c, there's no input - should be treated as 0:
  formula_N_1->setFormula("y = a*x + b + c"); 
  formula_N_1->setInputVariables("x,a,b");
  FormulaModule_N_1::process(formula_N_1, ins, outs, 0);
  result &= outs[0] == 3*2 + 5;

  // now there is a 4th input for c = 7 - 4 inputs, all used:
  formula_N_1->setFormula("y = a*x + b + c");
  formula_N_1->setInputVariables("x,a,b,c");  
  FormulaModule_N_1::process(formula_N_1, ins, outs, 0);
  result &= outs[0] == 3*2 + 5 + 7;

  // 4 inputs, all used:
  formula_N_1->setFormula("y = a + b + c + d");
  formula_N_1->setInputVariables("a,b,c,d");
  FormulaModule_N_1::process(formula_N_1, ins, outs, 0);
  result &= outs[0] == 2 + 3 + 5 + 7;

  // 4 inputs a,b,c,d - d is not used in the formula, but an x is used which doesn't exist as 
  // input (so is treated as x = 0) - so we should get b + c
  formula_N_1->setFormula("y = a*x + b + c");
  formula_N_1->setInputVariables("a,b,c,d");
  FormulaModule_N_1::process(formula_N_1, ins, outs, 0);
  result &= outs[0] == 3 + 5; //...but we get c + d - but that's coincidence - the expression                            
  // evaluator actaully still knows "x" which has a value 2, so we get: 2*2 + 3 + 5 = 12 - 
  // we should set all variables to zero that appear in the formula but don't have an input
  // associated with them...but we should also make sure, that the constants (e, pi, etc.) are
  // not touched

  // maybe we need to uncomment and use ExpressionEvaluator::initVariableList
  // ...Ok - seems to work




  return result;




  /*
  cm->addAudioConnection(formula_N_1, 0, audioOutput1, 0);
  cm->addAudioConnection(constant3,   0, formula_N_1,  0);  // x = 3
  cm->addAudioConnection(constant5,   0, formula_N_1,  1);  // a = 5
  cm->addAudioConnection(constant7,   0, formula_N_1,  2);  // b = 7
  */


}