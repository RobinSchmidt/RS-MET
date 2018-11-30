

bool romos::testFormulaModules()
{
  bool result = true;

  typedef romos::Module MDL;

  ContainerModule* cm = (ContainerModule*)moduleFactory.createModule("Container",
    "FormulaTest", 0, 0, true);

  MDL* constant3    = cm->addChildModule("Constant", "3", 2, 4, false, false);
  MDL* constant5    = cm->addChildModule("Constant", "5", 2, 6, false, false);
  MDL* constant7    = cm->addChildModule("Constant", "7", 2, 8, false, false);
  MDL* audioOutput1 = cm->addChildModule("AudioOutput", "Out1", 20, 6, false, false);



  FormulaModule_N_1* formula_N_1 = (FormulaModule_N_1*)
    cm->addChildModule("Formula_N_1", "Formula_N_1", 13, 6, false, false);
  cm->sortChildModuleArray();

  formula_N_1->setFormula("y = a*x + b");
  formula_N_1->setInputVariables("x,a,b");

  double ins[4] = { 2, 3, 5, 7};
  double outs[4];

  FormulaModule_N_1::process(formula_N_1, ins, outs, 0);

  result &= outs[0] == 3*2 + 5;




  /*
  cm->addAudioConnection(formula_N_1, 0, audioOutput1, 0);
  cm->addAudioConnection(constant3,   0, formula_N_1,  0);  // x = 3
  cm->addAudioConnection(constant5,   0, formula_N_1,  1);  // a = 5
  cm->addAudioConnection(constant7,   0, formula_N_1,  2);  // b = 7
  */

  return result;
}