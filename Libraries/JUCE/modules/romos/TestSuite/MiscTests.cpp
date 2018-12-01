

bool romos::testFormulaModules()
{
  bool result = true;

  typedef romos::Module MDL;
  ContainerModule* cm = (ContainerModule*)moduleFactory.createModule("Container",
    "FormulaTest", 0, 0, true);
  // actually, we don't really need this container...yet...but later we may want to make tests for
  // the updating of the connectivity when the input variables change...

  MDL* constant2 = cm->addChildModule("Constant", "2", 2, 4, false, false);
  MDL* constant3 = cm->addChildModule("Constant", "3", 2, 4, false, false);
  MDL* constant5 = cm->addChildModule("Constant", "5", 2, 6, false, false);
  MDL* constant7 = cm->addChildModule("Constant", "7", 2, 8, false, false);
  MDL* audioOutput1 = cm->addChildModule("AudioOutput", "Out1", 20, 6, false, false);

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
  result &= outs[0] == 3 + 5;

  // check, how the module connectivity behaves when inputs are added or removed...

  formula_N_1->setFormula("y = a + b + c + d");
  formula_N_1->setInputVariables("a,b,c,d");
  
  cm->addAudioConnection(formula_N_1, 0, audioOutput1, 0);
  cm->addAudioConnection(constant2, 0, formula_N_1,  0);  // a = 2
  cm->addAudioConnection(constant3, 0, formula_N_1,  1);  // b = 3
  cm->addAudioConnection(constant5, 0, formula_N_1,  2);  // c = 5
  cm->addAudioConnection(constant7, 0, formula_N_1,  3);  // d = 7

  std::vector<AudioConnection> inCons;
  inCons = formula_N_1->getIncomingAudioConnections();
  RAPT::rsAssert(inCons.size() == 4);
  result &= inCons.size() == 4;
  result &= inCons[0].getSourceModule() == constant2;
  result &= inCons[1].getSourceModule() == constant3;
  result &= inCons[2].getSourceModule() == constant5;
  result &= inCons[3].getSourceModule() == constant7;

  // remove 2nd input variable ("b") - this should result in the connection from 3 to b 
  // (pinIndex 1) to disappear, instead 5 should now be connected to pin 1 and 7 to pin 2:
  formula_N_1->setInputVariables("a,c,d");
  inCons = formula_N_1->getIncomingAudioConnections();
  RAPT::rsAssert(inCons.size() == 3);
  result &= inCons.size() == 3;
  result &= inCons[0].getSourceModule() == constant2;
  result &= inCons[1].getSourceModule() == constant5;
  result &= inCons[2].getSourceModule() == constant7;
  // this test still fails (as expected) - todo: implement pin re-assignment
  // ..oh - and we have memory leaks!

  //formula_N_1->setInputVariables("d,a,c"); // permute inputs, see if connections get permuted accordingly



  return result;
}