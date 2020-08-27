#include "romos_TestModuleBuilder.h"
//using namespace rsTestRomos;

namespace rsTestRomos
{

romos::Module* TestModuleBuilder::createWrappedAdder(const rosic::rsString &name, int x, int y, 
  bool polyphonic)
{
  romos::ContainerModule* m = (romos::ContainerModule*) romos::moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);

  // change argument order to: typeID, name, x, y, poly, sort - pass as last argument not the polyphonic
  // parameter that is passed in here, but either true or false, depending on what is actually desired
  romos::Module *audioInput1  = m->addChildModule("AudioInput",  "In1",  2,  2, false, false);
  romos::Module *audioInput2  = m->addChildModule("AudioInput",  "In2",  2,  5, false, false);
  romos::Module *adder1       = m->addChildModule("Adder",       "+",    9,  3, false, false);
  romos::Module *audioOutput1 = m->addChildModule("AudioOutput", "Out", 13,  3, false, false);

  m->addAudioConnection(audioInput1,  0, adder1,       0);
  m->addAudioConnection(audioInput2,  0, adder1,       1);
  m->addAudioConnection(adder1,       0, audioOutput1, 0);

  return m;  
}

romos::Module* TestModuleBuilder::createGain(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  romos::ContainerModule* m = (romos::ContainerModule*) romos::moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);

  romos::Module *audioInput1  = m->addChildModule("AudioInput",  "In",   1,  2, false, false);
  romos::Module *audioInput2  = m->addChildModule("AudioInput",  "G",    1,  5, false, false);
  romos::Module *multiplier1  = m->addChildModule("Multiplier",  "*",    7,  2, false, false);
  romos::Module *audioOutput1 = m->addChildModule("AudioOutput", "Out", 12,  2, false, false);

  m->addAudioConnection(audioInput1,  0, multiplier1,  0);
  m->addAudioConnection(audioInput2,  0, multiplier1,  1);
  m->addAudioConnection(multiplier1,  0, audioOutput1, 0);

  return m;
}

romos::Module* TestModuleBuilder::createSumDiff(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* m = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);

  romos::Module *audioInput1  = m->addChildModule("AudioInput",  "In1",   1,  2, false, false);
  romos::Module *audioInput2  = m->addChildModule("AudioInput",  "In2",   1,  7, false, false);
  romos::Module *adder1       = m->addChildModule("Adder",       "+",     9,  2, false, false);
  romos::Module *subtractor1  = m->addChildModule("Subtractor",  "-",     9,  6, false, false);
  romos::Module *audioOutput1 = m->addChildModule("AudioOutput", "Sum",  14,  2, false, false);
  romos::Module *audioOutput2 = m->addChildModule("AudioOutput", "Diff", 14,  6, false, false);

  m->addAudioConnection(audioInput1,  0, adder1,       0);
  m->addAudioConnection(audioInput2,  0, adder1,       1);
  m->addAudioConnection(audioInput1,  0, subtractor1,  0);
  m->addAudioConnection(audioInput2,  0, subtractor1,  1);
  m->addAudioConnection(adder1,       0, audioOutput1, 0);
  m->addAudioConnection(subtractor1,  0, audioOutput2, 0);

  return m;
}

romos::Module* TestModuleBuilder::createWrappedSumDiff(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* m = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);

  romos::Module *audioInput1  = m->addChildModule("AudioInput",  "In1",   1,  3, false, false);
  romos::Module *audioInput2  = m->addChildModule("AudioInput",  "In2",   1,  6, false, false);
  romos::Module *audioOutput1 = m->addChildModule("AudioOutput", "Out1", 16,  3, false, false);
  romos::Module *audioOutput2 = m->addChildModule("AudioOutput", "Out2", 16,  6, false, false);

  romos::Module *sumDiff = m->addChildModule(createSumDiff("SumDiff",  8,  4, false));

  m->sortChildModuleArray();
  
  m->addAudioConnection(audioInput1,  0, sumDiff,      0);
  m->addAudioConnection(audioInput2,  0, sumDiff,      1);
  m->addAudioConnection(sumDiff,      0, audioOutput1, 0);
  m->addAudioConnection(sumDiff,      1, audioOutput2, 0);

  return m;
}

romos::Module* TestModuleBuilder::createSumDiffProd(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* m = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);

  romos::Module *audioInput1  = m->addChildModule("AudioInput",  "In1",   1,  2, false, false);
  romos::Module *audioInput2  = m->addChildModule("AudioInput",  "In2",   1,  9, false, false);
  romos::Module *adder1       = m->addChildModule("Adder",       "+",     9,  2, false, false);
  romos::Module *subtractor1  = m->addChildModule("Subtractor",  "-",     9,  5, false, false);
  romos::Module *multiplier1  = m->addChildModule("Multiplier",  "*",     9,  8, false, false);
  romos::Module *audioOutput1 = m->addChildModule("AudioOutput", "Sum",  14,  2, false, false);
  romos::Module *audioOutput2 = m->addChildModule("AudioOutput", "Diff", 14,  5, false, false);
  romos::Module *audioOutput3 = m->addChildModule("AudioOutput", "Prod", 14,  8, false, false);

  m->sortChildModuleArray();

  m->addAudioConnection(audioInput1,  0, adder1,       0);
  m->addAudioConnection(audioInput2,  0, adder1,       1);
  m->addAudioConnection(audioInput1,  0, subtractor1,  0);
  m->addAudioConnection(audioInput2,  0, subtractor1,  1);
  m->addAudioConnection(audioInput1,  0, multiplier1,  0);
  m->addAudioConnection(audioInput2,  0, multiplier1,  1);
  m->addAudioConnection(adder1,       0, audioOutput1, 0);
  m->addAudioConnection(subtractor1,  0, audioOutput2, 0);
  m->addAudioConnection(multiplier1,  0, audioOutput3, 0);

  return m;
}

romos::Module* TestModuleBuilder::createWrappedSumDiffProd(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* m = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);

  romos::Module *audioInput1  = m->addChildModule("AudioInput",  "In1",   1,  2, false, false);
  romos::Module *audioInput2  = m->addChildModule("AudioInput",  "In2",   1,  7, false, false);
  romos::Module *audioOutput1 = m->addChildModule("AudioOutput", "Sum",  20,  2, false, false);
  romos::Module *audioOutput2 = m->addChildModule("AudioOutput", "Diff", 20,  5, false, false);
  romos::Module *audioOutput3 = m->addChildModule("AudioOutput", "Prod", 20,  8, false, false);

  romos::Module *sumDiffProd = m->addChildModule(createSumDiffProd("SumDiffProd",  8,  4, false));

  m->sortChildModuleArray();

  m->addAudioConnection(audioInput1,  0, sumDiffProd,  0);
  m->addAudioConnection(audioInput2,  0, sumDiffProd,  1);
  m->addAudioConnection(sumDiffProd,  0, audioOutput1, 0);
  m->addAudioConnection(sumDiffProd,  1, audioOutput2, 0);
  m->addAudioConnection(sumDiffProd,  2, audioOutput3, 0);

  return m;
}

romos::Module* TestModuleBuilder::createWrappedAdderN(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* m = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);

  romos::Module *audioInput1  = m->addChildModule("AudioInput",  "In1",    1,  2, false, false);
  romos::Module *adderN1      = m->addChildModule("AdderN",      "Sum",   10,  2, false, false);
  romos::Module *audioOutput1 = m->addChildModule("AudioOutput", "Out1",  20,  2, false, false);

  m->sortChildModuleArray();

  static const int numValuesToAdd = 10;
  for(int i = 0; i < numValuesToAdd; i++)
    m->addAudioConnection(audioInput1,  0, adderN1, i);
  m->addAudioConnection(adderN1, 0, audioOutput1,  0);

  return m;
}

romos::Module* TestModuleBuilder::createDifferences(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* m = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);

  romos::Module *audioInput1  = m->addChildModule("AudioInput",  "In1",   1,  2, false, false);
  romos::Module *audioInput2  = m->addChildModule("AudioInput",  "In2",   1,  7, false, false);
  romos::Module *audioInput3  = m->addChildModule("AudioInput",  "In3",   1, 12, false, false);
  romos::Module *subtractor1  = m->addChildModule("Subtractor",  "-",    12,  2, false, false);
  romos::Module *subtractor2  = m->addChildModule("Subtractor",  "-",    12,  5, false, false);
  romos::Module *subtractor3  = m->addChildModule("Subtractor",  "-",    12,  8, false, false);
  romos::Module *subtractor4  = m->addChildModule("Subtractor",  "-",    12, 11, false, false);
  romos::Module *audioOutput1 = m->addChildModule("AudioOutput", "Out1", 19,  2, false, false);
  romos::Module *audioOutput2 = m->addChildModule("AudioOutput", "Out2", 19,  5, false, false);
  romos::Module *audioOutput3 = m->addChildModule("AudioOutput", "Out3", 19,  8, false, false);
  romos::Module *audioOutput4 = m->addChildModule("AudioOutput", "Out4", 19, 11, false, false);

  m->sortChildModuleArray();

  m->addAudioConnection(audioInput1,  0, subtractor1,  0);
  m->addAudioConnection(audioInput2,  0, subtractor1,  1);
  m->addAudioConnection(audioInput2,  0, subtractor2,  0);
  m->addAudioConnection(audioInput1,  0, subtractor2,  1);
  m->addAudioConnection(audioInput3,  0, subtractor3,  0);
  m->addAudioConnection(audioInput2,  0, subtractor3,  1);
  m->addAudioConnection(audioInput2,  0, subtractor4,  0);
  m->addAudioConnection(audioInput3,  0, subtractor4,  1);
  m->addAudioConnection(subtractor1,  0, audioOutput1, 0);
  m->addAudioConnection(subtractor2,  0, audioOutput2, 0);
  m->addAudioConnection(subtractor3,  0, audioOutput3, 0);
  m->addAudioConnection(subtractor4,  0, audioOutput4, 0);

  return m;
}

romos::Module* TestModuleBuilder::createSummedDiffs(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* m = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);

  romos::Module *audioInput1  = m->addChildModule("AudioInput",  "In1",   1,  3, false, false);
  romos::Module *audioInput2  = m->addChildModule("AudioInput",  "In2",   1, 12, false, false);
  romos::Module *audioInput3  = m->addChildModule("AudioInput",  "In3",   1, 21, false, false);
  romos::Module *subtractor1  = m->addChildModule("Subtractor",  "-",    12,  3, false, false);
  romos::Module *subtractor2  = m->addChildModule("Subtractor",  "-",    12,  8, false, false);
  romos::Module *subtractor3  = m->addChildModule("Subtractor",  "-",    12, 15, false, false);
  romos::Module *subtractor4  = m->addChildModule("Subtractor",  "-",    12, 20, false, false);
  romos::Module *adderN1      = m->addChildModule("AdderN",      "+",    20,  3, false, false);
  romos::Module *adderN2      = m->addChildModule("AdderN",      "+",    20, 18, false, false);
  romos::Module *adderN3      = m->addChildModule("AdderN",      "+",    25,  8, false, false);
  romos::Module *adderN4      = m->addChildModule("AdderN",      "+",    25, 14, false, false);
  romos::Module *audioOutput1 = m->addChildModule("AudioOutput", "Out1", 30,  3, false, false);
  romos::Module *audioOutput2 = m->addChildModule("AudioOutput", "Out2", 30,  8, false, false);
  romos::Module *audioOutput3 = m->addChildModule("AudioOutput", "Out3", 30, 14, false, false);
  romos::Module *audioOutput4 = m->addChildModule("AudioOutput", "Out4", 30, 18, false, false);

  m->sortChildModuleArray();

  m->addAudioConnection(audioInput1,  0, subtractor1,  0);
  m->addAudioConnection(audioInput2,  0, subtractor1,  1);
  m->addAudioConnection(audioInput2,  0, subtractor2,  0);
  m->addAudioConnection(audioInput1,  0, subtractor2,  1);
  m->addAudioConnection(audioInput3,  0, subtractor3,  0);
  m->addAudioConnection(audioInput2,  0, subtractor3,  1);
  m->addAudioConnection(audioInput2,  0, subtractor4,  0);
  m->addAudioConnection(audioInput3,  0, subtractor4,  1);
  m->addAudioConnection(subtractor1,  0, adderN1,      0);
  m->addAudioConnection(subtractor2,  0, adderN1,      1);
  m->addAudioConnection(subtractor4,  0, adderN1,      2);
  m->addAudioConnection(subtractor1,  0, adderN2,      0);
  m->addAudioConnection(subtractor3,  0, adderN2,      1);
  m->addAudioConnection(subtractor4,  0, adderN2,      2);
  m->addAudioConnection(subtractor2,  0, adderN3,      0);
  m->addAudioConnection(subtractor4,  0, adderN3,      1);
  m->addAudioConnection(subtractor1,  0, adderN4,      0);
  m->addAudioConnection(subtractor3,  0, adderN4,      1);
  m->addAudioConnection(adderN1,      0, audioOutput1, 0);
  m->addAudioConnection(adderN3,      0, audioOutput2, 0);
  m->addAudioConnection(adderN4,      0, audioOutput3, 0);
  m->addAudioConnection(adderN2,      0, audioOutput4, 0);

  return m;
}

romos::Module* TestModuleBuilder::createMovingAverage(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* m = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);

  romos::Module *audioInput1  = m->addChildModule("AudioInput",  "In",   1,  2, false, false);
  romos::Module *audioInput2  = m->addChildModule("AudioInput",  "b0",   1,  5, false, false);
  romos::Module *audioInput3  = m->addChildModule("AudioInput",  "b1",   1,  8, false, false);
  romos::Module *unitDelay1   = m->addChildModule("UnitDelay",   "D",    8,  6, false, false);
  romos::Module *multiplier1  = m->addChildModule("Multiplier",  "*",   11,  2, false, false);
  romos::Module *multiplier2  = m->addChildModule("Multiplier",  "*",   11,  7, false, false);
  romos::Module *adder1       = m->addChildModule("Adder",       "+",   15,  4, false, false);
  romos::Module *audioOutput1 = m->addChildModule("AudioOutput", "Out", 19,  4, false, false);

  m->addAudioConnection(audioInput1,  0, unitDelay1,   0);
  m->addAudioConnection(audioInput2,  0, multiplier1,  1);
  m->addAudioConnection(audioInput1,  0, multiplier1,  0);
  m->addAudioConnection(audioInput3,  0, multiplier2,  1);
  m->addAudioConnection(unitDelay1,   0, multiplier2,  0);
  m->addAudioConnection(multiplier1,  0, adder1,       0);
  m->addAudioConnection(multiplier2,  0, adder1,       1);
  m->addAudioConnection(adder1,       0, audioOutput1, 0);

  return m;
}

romos::Module* TestModuleBuilder::createDelayedConnection( const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* m = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *audioInput1  = m->addChildModule("AudioInput",  "In",   2,  2, false, false);
  romos::Module *identity1    = m->addChildModule("Identity",    "",    10,  2, false, false);
  romos::Module *identity2    = m->addChildModule("Identity",    "",     8,  5, false, false);
  romos::Module *audioOutput1 = m->addChildModule("AudioOutput", "Out", 20,  5, false, false);

  m->sortChildModuleArray();

  m->addAudioConnection(audioInput1,  0, identity1,     0);
  m->addAudioConnection(identity1,    0, identity2,     0);
  m->addAudioConnection(identity2,    0, audioOutput1,  0);

  return m;
}

romos::Module* TestModuleBuilder::createLeakyIntegrator(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *audioInput1  = module->addChildModule("AudioInput",  "In",   2,  8, false, false);
  romos::Module *audioInput2  = module->addChildModule("AudioInput",  "c",    2, 11, false, false);
  romos::Module *constant1    = module->addChildModule("Constant",    "1",    2,  5, false, false);
  romos::Module *subtractor1  = module->addChildModule("Subtractor",  "-",    8,  6, false, false);
  romos::Module *multiplier1  = module->addChildModule("Multiplier",  "*",   12,  5, false, false);
  romos::Module *multiplier2  = module->addChildModule("Multiplier",  "*",   12,  9, false, false);
  romos::Module *adder1       = module->addChildModule("Adder",       "+",   18,  7, false, false);
  romos::Module *identity1    = module->addChildModule("Identity",    "",    22,  2, false, false);
  romos::Module *audioOutput1 = module->addChildModule("AudioOutput", "Out", 23,  7, false, false);

  // using the function without the  doesn't work yet (crashes, because in one test, a 
  // module is searched for via its old type-id)

  module->sortChildModuleArray();

  module->addAudioConnection(constant1,    0, subtractor1,  0);
  module->addAudioConnection(audioInput2,  0, subtractor1,  1);
  module->addAudioConnection(subtractor1,  0, multiplier1,  1);
  module->addAudioConnection(identity1,    0, multiplier1,  0);
  module->addAudioConnection(audioInput2,  0, multiplier2,  0);
  module->addAudioConnection(audioInput1,  0, multiplier2,  1);
  module->addAudioConnection(multiplier1,  0, adder1,       0);
  module->addAudioConnection(multiplier2,  0, adder1,       1);
  module->addAudioConnection(adder1,       0, identity1,    0);
  module->addAudioConnection(adder1,       0, audioOutput1, 0);

  return module;
}

romos::Module* TestModuleBuilder::createTestFilter1(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *audioInput1  = module->addChildModule(("AudioInput"),  "In",    2,  2, false, false);
  romos::Module *audioInput2  = module->addChildModule(("AudioInput"),  "b0",    2,  5, false, false);
  romos::Module *audioInput3  = module->addChildModule(("AudioInput"),  "b1",    2,  8, false, false);
  romos::Module *audioInput4  = module->addChildModule(("AudioInput"),  "c",     2, 12, false, false);
  romos::Module *adder1       = module->addChildModule(("Adder"),       "+",    27,  4, false, false);
  romos::Module *subtractor1  = module->addChildModule(("Subtractor"),  "-",    27,  7, false, false);
  romos::Module *multiplier1  = module->addChildModule(("Multiplier"),  "*",    27, 10, false, false);
  romos::Module *audioOutput1 = module->addChildModule(("AudioOutput"), "Sum",  31,  4, false, false);
  romos::Module *audioOutput2 = module->addChildModule(("AudioOutput"), "Diff", 31,  7, false, false);
  romos::Module *audioOutput3 = module->addChildModule(("AudioOutput"), "Prod", 31, 10, false, false);

  romos::Module *movingAverage   = module->addChildModule(createMovingAverage(  "MovingAverage",   10,  4, false));
  romos::Module *leakyIntegrator = module->addChildModule(createLeakyIntegrator("LeakyIntegrator", 10, 11, false));

  module->sortChildModuleArray();

  module->addAudioConnection(audioInput1,     0, movingAverage,   0);
  module->addAudioConnection(audioInput2,     0, movingAverage,   1);
  module->addAudioConnection(audioInput3,     0, movingAverage,   2);
  module->addAudioConnection(audioInput1,     0, leakyIntegrator, 0);
  module->addAudioConnection(audioInput4,     0, leakyIntegrator, 1);
  module->addAudioConnection(movingAverage,   0, adder1,          0);
  module->addAudioConnection(leakyIntegrator, 0, adder1,          1);
  module->addAudioConnection(movingAverage,   0, subtractor1,     0);
  module->addAudioConnection(leakyIntegrator, 0, subtractor1,     1);
  module->addAudioConnection(movingAverage,   0, multiplier1,     0);
  module->addAudioConnection(leakyIntegrator, 0, multiplier1,     1);
  module->addAudioConnection(adder1,          0, audioOutput1,    0);
  module->addAudioConnection(subtractor1,     0, audioOutput2,    0);
  module->addAudioConnection(multiplier1,     0, audioOutput3,    0);

  return module;
}

romos::Module* TestModuleBuilder::createBiquadMacro(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *audioInput1  = module->addChildModule(("AudioInput"),  "In",   1,  3, false, false);
  romos::Module *audioInput2  = module->addChildModule(("AudioInput"),  "b0",   1,  5, false, false);
  romos::Module *audioInput3  = module->addChildModule(("AudioInput"),  "b1",   1,  9, false, false);
  romos::Module *audioInput4  = module->addChildModule(("AudioInput"),  "b2",   1, 13, false, false);
  romos::Module *audioInput5  = module->addChildModule(("AudioInput"),  "a1",   1, 17, false, false);
  romos::Module *audioInput6  = module->addChildModule(("AudioInput"),  "a2",   1, 21, false, false);
  romos::Module *unitDelay1   = module->addChildModule(("UnitDelay"),   "D",    7,  7, false, false);
  romos::Module *unitDelay2   = module->addChildModule(("UnitDelay"),   "D",   10, 11, false, false);
  romos::Module *unaryMinus1  = module->addChildModule(("UnaryMinus"),  "-",   11, 17, false, false);
  romos::Module *unitDelay3   = module->addChildModule(("UnitDelay"),   "D",   11, 19, false, false);
  romos::Module *unaryMinus2  = module->addChildModule(("UnaryMinus"),  "-",   11, 21, false, false);
  romos::Module *multiplier1  = module->addChildModule(("Multiplier"),  "*",   15,  3, false, false);
  romos::Module *multiplier2  = module->addChildModule(("Multiplier"),  "*",   15,  7, false, false);
  romos::Module *multiplier3  = module->addChildModule(("Multiplier"),  "*",   15, 11, false, false);
  romos::Module *multiplier4  = module->addChildModule(("Multiplier"),  "*",   15, 15, false, false);
  romos::Module *multiplier5  = module->addChildModule(("Multiplier"),  "*",   15, 19, false, false);
  romos::Module *adderN1      = module->addChildModule(("AdderN"),      "+",   23,  3, false, false);
  romos::Module *identity1    = module->addChildModule(("Identity"),    "",    29,  7, false, false);
  romos::Module *audioOutput1 = module->addChildModule(("AudioOutput"), "Out", 30,  3, false, false);

  module->sortChildModuleArray();

  module->addAudioConnection(audioInput1,  0, unitDelay1,   0);
  module->addAudioConnection(unitDelay1,   0, unitDelay2,   0);
  module->addAudioConnection(audioInput5,  0, unaryMinus1,  0);
  module->addAudioConnection(identity1,    0, unitDelay3,   0);
  module->addAudioConnection(audioInput6,  0, unaryMinus2,  0);
  module->addAudioConnection(audioInput1,  0, multiplier1,  0);
  module->addAudioConnection(audioInput2,  0, multiplier1,  1);
  module->addAudioConnection(unitDelay1,   0, multiplier2,  0);
  module->addAudioConnection(audioInput3,  0, multiplier2,  1);
  module->addAudioConnection(unitDelay2,   0, multiplier3,  0);
  module->addAudioConnection(audioInput4,  0, multiplier3,  1);
  module->addAudioConnection(unaryMinus1,  0, multiplier4,  1);
  module->addAudioConnection(identity1,    0, multiplier4,  0);
  module->addAudioConnection(unaryMinus2,  0, multiplier5,  1);
  module->addAudioConnection(unitDelay3,   0, multiplier5,  0);
  module->addAudioConnection(multiplier1,  0, adderN1,      0);
  module->addAudioConnection(multiplier2,  0, adderN1,      1);
  module->addAudioConnection(multiplier3,  0, adderN1,      2);
  module->addAudioConnection(multiplier4,  0, adderN1,      3);
  module->addAudioConnection(multiplier5,  0, adderN1,      4);
  module->addAudioConnection(adderN1,      0, identity1,    0);
  module->addAudioConnection(adderN1,      0, audioOutput1, 0);

  return module;
}


/*
romos::Module* TestModuleBuilder::createAddedConstants(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule *module = (ContainerModule*) ModuleFactory::createModule(ModuleTypeRegistry::CONTAINER, name, x, y, polyphonic);

  romos::Module *audioOutput1 = module->addChildModule(("AudioOutput"), "Out1", 19, 10, false, false);

  int xStart     = 1;
  int yStart     = 2;
  int xStep      = 4;
  int yStep      = 3;
  int numColumns = 4;
  int numRows    = 6; 
  for(int i=1; i<=numRows; i++)
  {
    for(int j=1; j<=numColumns; j++)
    {
      //ConstantModule *constantModule = new ConstantModule(NULL, polyphonic);
      Module *constantModule = ModuleFactory::createModule(ModuleTypeRegistry::CONSTANT);
      constantModule->setModuleName( String(i) + String(j) );
      constantModule->setPositionXY(xStart + (j-1)*xStep, yStart + (i-1)*yStep);
      module->addChildModule(constantModule);
      module->addAudioConnection(constantModule, 0, audioOutput1, 0);
    }
  }
  
  return module;
}
*/

romos::Module* TestModuleBuilder::createPinSortingInner(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *audioInput1  = module->addChildModule(("AudioInput"),  "1",    1,  2, false, false);
  romos::Module *audioInput2  = module->addChildModule(("AudioInput"),  "2",    1,  6, false, false);
  romos::Module *audioInput3  = module->addChildModule(("AudioInput"),  "3",    1, 11, false, false);
  romos::Module *subtractor1  = module->addChildModule(("Subtractor"),  "-",    7,  2, false, false);
  romos::Module *subtractor2  = module->addChildModule(("Subtractor"),  "-",    7,  6, false, false);
  romos::Module *subtractor3  = module->addChildModule(("Subtractor"),  "-",    7, 10, false, false);
  romos::Module *audioOutput1 = module->addChildModule(("AudioOutput"), "1m2", 11,  2, false, false);
  romos::Module *audioOutput2 = module->addChildModule(("AudioOutput"), "1m3", 11,  6, false, false);
  romos::Module *audioOutput3 = module->addChildModule(("AudioOutput"), "2m3", 11, 10, false, false);

  module->addAudioConnection(audioInput1,  0, subtractor1,  0);
  module->addAudioConnection(audioInput2,  0, subtractor1,  1);
  module->addAudioConnection(audioInput1,  0, subtractor2,  0);
  module->addAudioConnection(audioInput3,  0, subtractor2,  1);
  module->addAudioConnection(audioInput3,  0, subtractor3,  1);
  module->addAudioConnection(audioInput2,  0, subtractor3,  0);
  module->addAudioConnection(subtractor1,  0, audioOutput1, 0);
  module->addAudioConnection(subtractor2,  0, audioOutput2, 0);
  module->addAudioConnection(subtractor3,  0, audioOutput3, 0);

  return module;
}

romos::Module* TestModuleBuilder::createPinSortTest(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *audioInput1  = module->addChildModule(("AudioInput"),  "1",    1,  2, false, false);
  romos::Module *audioInput2  = module->addChildModule(("AudioInput"),  "2",    1,  6, false, false);
  romos::Module *audioInput3  = module->addChildModule(("AudioInput"),  "3",    1, 10, false, false);
  romos::Module *audioOutput1 = module->addChildModule(("AudioOutput"), "1m2", 20,  2, false, false);
  romos::Module *audioOutput2 = module->addChildModule(("AudioOutput"), "1m3", 20,  6, false, false);
  romos::Module *audioOutput3 = module->addChildModule(("AudioOutput"), "2m3", 20, 10, false, false);

  romos::Module *pinSortingInner = module->addChildModule(createPinSortingInner("PinSortingInner",  7,  5, false));

  module->sortChildModuleArray();

  module->addAudioConnection(audioInput1,     0, pinSortingInner, 0);
  module->addAudioConnection(audioInput2,     0, pinSortingInner, 1);
  module->addAudioConnection(audioInput3,     0, pinSortingInner, 2);
  module->addAudioConnection(pinSortingInner, 0, audioOutput1,    0);
  module->addAudioConnection(pinSortingInner, 1, audioOutput2,    0);
  module->addAudioConnection(pinSortingInner, 2, audioOutput3,    0);

  return module;
}

romos::Module* TestModuleBuilder::createDifferencer(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *audioInput1  = module->addChildModule(("AudioInput"),       "In1",               1,  4, false, false);
  romos::Module *constant1    = module->addChildModule(("Constant"),         "-1",                3,  8, false, false);
  romos::Module *constant2    = module->addChildModule(("Constant"),         "1",                 4,  6, false, false);
  romos::Module *filter1p1z1  = module->addChildModule(("FirstOrderFilter"), "FirstOrderFilter",  8,  4, false, false);
  romos::Module *audioOutput1 = module->addChildModule(("AudioOutput"),      "Out2",             17,  4, false, false);

  module->addAudioConnection(constant2,    0, filter1p1z1,  1);
  module->addAudioConnection(constant1,    0, filter1p1z1,  2);
  module->addAudioConnection(audioInput1,  0, filter1p1z1,  0);
  module->addAudioConnection(filter1p1z1,  0, audioOutput1, 0);

  return module;
}

romos::Module* TestModuleBuilder::createImpulse(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *noteGate1    = module->addChildModule(("NoteGate"),    "NoteGate",  1,  4, false, false);
  romos::Module *constant1    = module->addChildModule(("Constant"),    "0",        16,  7, false, false);
  romos::Module *constant2    = module->addChildModule(("Constant"),    "1",        16, 10, false, false);
  romos::Module *clipper1     = module->addChildModule(("Clipper"),     "Clipper",  21,  4, false, false);
  romos::Module *audioOutput1 = module->addChildModule(("AudioOutput"), "Out",      28,  4, false, false);

  romos::Module *differencer = module->addChildModule(createDifferencer("Differencer", 10,  4, false));

  module->sortChildModuleArray();

  module->addAudioConnection(noteGate1,    0, differencer,  0);
  module->addAudioConnection(differencer,  0, clipper1,     0);
  module->addAudioConnection(constant1,    0, clipper1,     1);
  module->addAudioConnection(constant2,    0, clipper1,     2);
  module->addAudioConnection(clipper1,     0, audioOutput1, 0);

  return module;
}

romos::Module* TestModuleBuilder::createBlip(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *noteFrequency1  = module->addChildModule(("NoteFrequency"),          "NoteFrequency",  1, 11, false, false);
  romos::Module *constant1       = module->addChildModule(("Constant"),               "20",             1, 13, false, false);
  romos::Module *biquadDesigner1 = module->addChildModule(("BiquadDesigner"),         "Designer",      14, 10, false, false);
  romos::Module *biquad1         = module->addChildModule(("Biquad"),                 "Biquad",        24,  9, false, false);
  romos::Module *audioOutput1    = module->addChildModule(("AudioOutput"),            "Out1",          33,  9, false, false);

  romos::Module *impulse = module->addChildModule(createImpulse("Impulse", 15,  4, false));

  module->sortChildModuleArray();

  ((romos::ModuleWithParameters*) biquadDesigner1)->setParameter("Mode", "Bandpass, const. skirt, BLT");

  module->addAudioConnection(constant1,       0, biquadDesigner1, 1);
  module->addAudioConnection(noteFrequency1,  0, biquadDesigner1, 0);
  module->addAudioConnection(impulse,         0, biquad1,         0);
  module->addAudioConnection(biquadDesigner1, 0, biquad1,         1);
  module->addAudioConnection(biquadDesigner1, 1, biquad1,         2);
  module->addAudioConnection(biquadDesigner1, 2, biquad1,         3);
  module->addAudioConnection(biquadDesigner1, 3, biquad1,         4);
  module->addAudioConnection(biquadDesigner1, 4, biquad1,         5);
  module->addAudioConnection(biquad1,         0, audioOutput1,    0);

  return module;
}

romos::Module* TestModuleBuilder::createMonoToPoly(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *constant1    = module->addChildModule(("Constant"),    "1",    1,  2, false, false);
  romos::Module *unaryMinus1  = module->addChildModule(("UnaryMinus"),  "-",    5,  2, false, false);
  romos::Module *audioOutput1 = module->addChildModule(("AudioOutput"), "Out",  9,  2, false, false);

  constant1->setPolyphonic(  false);
  unaryMinus1->setPolyphonic(true );

  module->addAudioConnection(constant1,    0, unaryMinus1,  0);
  module->addAudioConnection(unaryMinus1,  0, audioOutput1, 0);

  return module;
}

romos::Module* TestModuleBuilder::createVoiceCombiner(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *constant1      = module->addChildModule(("Constant"),      "1",    1,  2, false, false);
  romos::Module *voiceCombiner1 = module->addChildModule(("VoiceCombiner"), "}",    5,  2, false, false);
  romos::Module *audioOutput1   = module->addChildModule(("AudioOutput"),   "Out",  9,  2, false, false);

  constant1->setPolyphonic(     true );
  voiceCombiner1->setPolyphonic(false);

  module->addAudioConnection(constant1,       0, voiceCombiner1,  0);
  module->addAudioConnection(voiceCombiner1,  0, audioOutput1,    0);

  return module;
}

romos::Module* TestModuleBuilder::createIn3Out5(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *audioInput1  = module->addChildModule(("AudioInput"),  "In1",   1,  2, false, false);
  romos::Module *audioInput2  = module->addChildModule(("AudioInput"),  "In2",   1,  8, false, false);
  romos::Module *audioInput3  = module->addChildModule(("AudioInput"),  "In3",   1, 15, false, false);
  romos::Module *subtractor1  = module->addChildModule(("Subtractor"),  "-",    13,  2, false, false);
  romos::Module *subtractor2  = module->addChildModule(("Subtractor"),  "-",    13,  8, false, false);
  romos::Module *subtractor3  = module->addChildModule(("Subtractor"),  "-",    13, 14, false, false);
  romos::Module *subtractor4  = module->addChildModule(("Subtractor"),  "-",    19,  5, false, false);
  romos::Module *subtractor5  = module->addChildModule(("Subtractor"),  "-",    19, 11, false, false);
  romos::Module *audioOutput1 = module->addChildModule(("AudioOutput"), "Out1", 25,  2, false, false);
  romos::Module *audioOutput2 = module->addChildModule(("AudioOutput"), "Out2", 25,  5, false, false);
  romos::Module *audioOutput3 = module->addChildModule(("AudioOutput"), "Out3", 25,  8, false, false);
  romos::Module *audioOutput4 = module->addChildModule(("AudioOutput"), "Out4", 25, 11, false, false);
  romos::Module *audioOutput5 = module->addChildModule(("AudioOutput"), "Out5", 25, 14, false, false);

  module->sortChildModuleArray();

  module->addAudioConnection(audioInput1,  0, subtractor1,  0);
  module->addAudioConnection(audioInput2,  0, subtractor1,  1);
  module->addAudioConnection(audioInput1,  0, subtractor2,  0);
  module->addAudioConnection(audioInput3,  0, subtractor2,  1);
  module->addAudioConnection(audioInput3,  0, subtractor3,  1);
  module->addAudioConnection(audioInput2,  0, subtractor3,  0);
  module->addAudioConnection(subtractor1,  0, subtractor4,  0);
  module->addAudioConnection(subtractor2,  0, subtractor4,  1);
  module->addAudioConnection(subtractor2,  0, subtractor5,  0);
  module->addAudioConnection(subtractor3,  0, subtractor5,  1);
  module->addAudioConnection(subtractor1,  0, audioOutput1, 0);
  module->addAudioConnection(subtractor4,  0, audioOutput2, 0);
  module->addAudioConnection(subtractor2,  0, audioOutput3, 0);
  module->addAudioConnection(subtractor5,  0, audioOutput4, 0);
  module->addAudioConnection(subtractor3,  0, audioOutput5, 0);

  return module;
}

romos::Module* TestModuleBuilder::createGateAndKill(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *noteGate1    = module->addChildModule(("NoteGate"),    "NoteGate",     2,  2, true,  false);
  romos::Module *constant1    = module->addChildModule(("Constant"),    "0.001",        5,  8, true,  false);
  romos::Module *constant2    = module->addChildModule(("Constant"),    "0.001",        5, 11, true,  false);
  romos::Module *voiceKiller1 = module->addChildModule(("VoiceKiller"), "VoiceKiller", 15,  7, true,  false);
  romos::Module *audioOutput1 = module->addChildModule(("AudioOutput"), "Out1",        24,  2, true,  false);

  module->sortChildModuleArray();

  module->addAudioConnection(constant1,    0, voiceKiller1, 1);
  module->addAudioConnection(constant2,    0, voiceKiller1, 2);
  module->addAudioConnection(noteGate1,    0, voiceKiller1, 0);
  module->addAudioConnection(noteGate1,    0, audioOutput1, 0);

  return module;
}


romos::Module* TestModuleBuilder::createTriggerAndKill(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *noteOnTrigger1 = module->addChildModule(("NoteOnTrigger"), "NoteOnTrigger",  1,  2, false, false);
  romos::Module *constant1      = module->addChildModule(("Constant"),      "0.001",          9,  6, false, false);
  romos::Module *constant2      = module->addChildModule(("Constant"),      "0.001",          9,  8, false, false);
  romos::Module *voiceKiller1   = module->addChildModule(("VoiceKiller"),   "VoiceKiller",   18,  6, false, false);
  romos::Module *audioOutput1   = module->addChildModule(("AudioOutput"),   "Out1",          29,  2, false, false);

  module->addAudioConnection(noteOnTrigger1, 0, voiceKiller1,   0);
  module->addAudioConnection(constant1,      0, voiceKiller1,   1);
  module->addAudioConnection(constant2,      0, voiceKiller1,   2);
  module->addAudioConnection(noteOnTrigger1, 0, audioOutput1,   0);

  return module;
}

romos::Module* TestModuleBuilder::createGatedNoteFrequency(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *noteGate1      = module->addChildModule(("NoteGate"),      "NoteGate",       1,  2, false, false);
  romos::Module *noteFrequency1 = module->addChildModule(("NoteFrequency"), "NoteFrequency",  1,  5, false, false);
  romos::Module *multiplier1    = module->addChildModule(("Multiplier"),    "*",             16,  2, false, false);
  romos::Module *audioOutput1   = module->addChildModule(("AudioOutput"),   "Out",           22,  2, false, false);

  module->addAudioConnection(noteGate1,      0, multiplier1,    0);
  module->addAudioConnection(noteFrequency1, 0, multiplier1,    1);
  module->addAudioConnection(multiplier1,    0, audioOutput1,   0);

  return module;
}





romos::Module* TestModuleBuilder::createVoiceKill(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *audioInput1  = module->addChildModule(("AudioInput"),  "In1",          2,  2, true,  false);
  romos::Module *constant1    = module->addChildModule(("Constant"),    "0.001",        2,  5, true,  false);
  romos::Module *constant2    = module->addChildModule(("Constant"),    "0.01",         2,  8, true,  false);
  romos::Module *voiceKiller1 = module->addChildModule(("VoiceKiller"), "VoiceKiller", 12,  4, true,  false);

  module->sortChildModuleArray();

  module->addAudioConnection(audioInput1,  0, voiceKiller1, 0);
  module->addAudioConnection(constant1,    0, voiceKiller1, 1);
  module->addAudioConnection(constant2,    0, voiceKiller1, 2);

  return module;
}

romos::Module* TestModuleBuilder::createIn1Out2(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *audioInput1  = module->addChildModule(("AudioInput"),  "In1",   2,  4, true,  false);
  romos::Module *audioOutput1 = module->addChildModule(("AudioOutput"), "Out1", 11,  2, true,  false);
  romos::Module *audioOutput2 = module->addChildModule(("AudioOutput"), "Out2", 11,  6, true,  false);

  module->sortChildModuleArray();

  module->addAudioConnection(audioInput1,  0, audioOutput1, 0);
  module->addAudioConnection(audioInput1,  0, audioOutput2, 0);

  return module;
}

romos::Module* TestModuleBuilder::createPolyBlipStereo(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *voiceCombiner1 = module->addChildModule(("VoiceCombiner"), "}",       22,  3, false, false);
  romos::Module *voiceCombiner2 = module->addChildModule(("VoiceCombiner"), "}",       22,  6, false, false);
  romos::Module *audioOutput1   = module->addChildModule(("AudioOutput"),   "Out1",    27,  3, false, false);
  romos::Module *audioOutput2   = module->addChildModule(("AudioOutput"),   "Out2",    27,  6, false, false);

  romos::Module *testBlip  = module->addChildModule(createBlip(     "TestBlip",   3,  4, true ));
  romos::Module *in1Out2   = module->addChildModule(createIn1Out2(  "In1Out2",   13,  4, true ));
  romos::Module *voiceKill = module->addChildModule(createVoiceKill("VoiceKill", 14, 12, true ));

  ((ContainerModule*) testBlip)->setPolyphonicRecursively(true);

  module->sortChildModuleArray();

  module->addAudioConnection(testBlip,       0, in1Out2,        0);
  module->addAudioConnection(testBlip,       0, voiceKill,      0);
  module->addAudioConnection(in1Out2,        0, voiceCombiner1, 0);
  module->addAudioConnection(in1Out2,        1, voiceCombiner2, 0);
  module->addAudioConnection(voiceCombiner1, 0, audioOutput1,   0);
  module->addAudioConnection(voiceCombiner2, 0, audioOutput2,   0);

  return module;
}








romos::Module* TestModuleBuilder::createNoteFilter(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *audioInput1     = module->addChildModule(("AudioInput"),     "In1",            16,  4, true,  false);
  romos::Module *audioInput2     = module->addChildModule(("AudioInput"),     "In2",            16, 18, true,  false);
  romos::Module *noteFrequency1  = module->addChildModule(("NoteFrequency"),  "NoteFreq",        2, 10, true,  false);
  romos::Module *constant1       = module->addChildModule(("Constant"),       "50",              2, 12, true,  false);
  romos::Module *biquadDesigner1 = module->addChildModule(("BiquadDesigner"), "BiquadDesigner", 11, 10, true,  false);
  romos::Module *biquad1         = module->addChildModule(("Biquad"),         "Biquad",         30,  4, true,  false);
  romos::Module *biquad2         = module->addChildModule(("Biquad"),         "Biquad",         30, 14, true,  false);
  romos::Module *audioOutput1    = module->addChildModule(("AudioOutput"),    "Out1",           40,  4, true,  false);
  romos::Module *audioOutput2    = module->addChildModule(("AudioOutput"),    "Out2",           40, 14, true,  false);

  module->sortChildModuleArray();

  ((romos::ModuleWithParameters*) biquadDesigner1)->setParameter("Mode", "Bandpass, const. skirt, BLT");

  module->addAudioConnection(noteFrequency1,  0, biquadDesigner1, 0);
  module->addAudioConnection(constant1,       0, biquadDesigner1, 1);
  module->addAudioConnection(audioInput1,     0, biquad1,         0);
  module->addAudioConnection(biquadDesigner1, 0, biquad1,         1);
  module->addAudioConnection(biquadDesigner1, 1, biquad1,         2);
  module->addAudioConnection(biquadDesigner1, 2, biquad1,         3);
  module->addAudioConnection(biquadDesigner1, 3, biquad1,         4);
  module->addAudioConnection(biquadDesigner1, 4, biquad1,         5);
  module->addAudioConnection(audioInput2,     0, biquad2,         0);
  module->addAudioConnection(biquadDesigner1, 0, biquad2,         1);
  module->addAudioConnection(biquadDesigner1, 1, biquad2,         2);
  module->addAudioConnection(biquadDesigner1, 2, biquad2,         3);
  module->addAudioConnection(biquadDesigner1, 3, biquad2,         4);
  module->addAudioConnection(biquadDesigner1, 4, biquad2,         5);
  module->addAudioConnection(biquad1,         0, audioOutput1,    0);
  module->addAudioConnection(biquad2,         0, audioOutput2,    0);

  return module;
}

romos::Module* TestModuleBuilder::createGatedNoise(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *constant1    = module->addChildModule(("Constant"),    "0",           2,  4, true,  false);
  romos::Module *constant2    = module->addChildModule(("Constant"),    "1",           2, 14, true,  false);
  romos::Module *whiteNoise1  = module->addChildModule(("WhiteNoise"),  "WhiteNoise",  6,  4, true,  false);
  romos::Module *noteGate1    = module->addChildModule(("NoteGate"),    "NoteGate",    6,  8, true,  false);
  romos::Module *whiteNoise2  = module->addChildModule(("WhiteNoise"),  "WhiteNoise",  6, 14, true,  false);
  romos::Module *multiplier1  = module->addChildModule(("Multiplier"),  "*",          18,  4, true,  false);
  romos::Module *multiplier2  = module->addChildModule(("Multiplier"),  "*",          18, 13, true,  false);
  romos::Module *audioOutput1 = module->addChildModule(("AudioOutput"), "Out1",       23,  4, true,  false);
  romos::Module *audioOutput2 = module->addChildModule(("AudioOutput"), "Out2",       23, 13, true,  false);

  module->sortChildModuleArray();

  ((romos::ModuleWithParameters*) whiteNoise2)->setParameter("Seed", "1.0");

  module->addAudioConnection(constant1,   0, whiteNoise1,  0);
  module->addAudioConnection(constant2,   0, whiteNoise2,  0);
  module->addAudioConnection(whiteNoise1, 0, multiplier1,  0);
  module->addAudioConnection(noteGate1,   0, multiplier1,  1);
  module->addAudioConnection(noteGate1,   0, multiplier2,  0);
  module->addAudioConnection(whiteNoise2, 0, multiplier2,  1);
  module->addAudioConnection(multiplier1, 0, audioOutput1, 0);
  module->addAudioConnection(multiplier2, 0, audioOutput2, 0);

  return module;
}

romos::Module* TestModuleBuilder::createNoiseFlute(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *voiceCombiner1 = module->addChildModule(("VoiceCombiner"), "}",          25,  4, false, false);
  romos::Module *constant1      = module->addChildModule(("Constant"),      "0.125",      25,  7, false, false);
  romos::Module *voiceCombiner2 = module->addChildModule(("VoiceCombiner"), "}",          25, 10, false, false);
  romos::Module *multiplier1    = module->addChildModule(("Multiplier"),    "*",          31,  4, false, false);
  romos::Module *multiplier2    = module->addChildModule(("Multiplier"),    "*",          31,  9, false, false);
  romos::Module *audioOutput1   = module->addChildModule(("AudioOutput"),   "Out1",       36,  4, false, false);
  romos::Module *audioOutput2   = module->addChildModule(("AudioOutput"),   "Out2",       36,  9, false, false);

  romos::Module *gatedNoise = module->addChildModule(createGatedNoise("GatedNoise",  1,  4, true ));
  romos::Module *noteFilter = module->addChildModule(createNoteFilter("NoteFilter", 13,  4, true ));
  romos::Module *voiceKill  = module->addChildModule(createVoiceKill( "VoiceKill",  25, 16, true ));

  module->sortChildModuleArray();

  module->addAudioConnection(gatedNoise,     0, noteFilter,     0);
  module->addAudioConnection(gatedNoise,     1, noteFilter,     1);
  module->addAudioConnection(noteFilter,     0, voiceCombiner1, 0);
  module->addAudioConnection(noteFilter,     1, voiceCombiner2, 0);
  module->addAudioConnection(noteFilter,     1, voiceKill,      0);
  module->addAudioConnection(voiceCombiner1, 0, multiplier1,    0);
  module->addAudioConnection(constant1,      0, multiplier1,    1);
  module->addAudioConnection(constant1,      0, multiplier2,    0);
  module->addAudioConnection(voiceCombiner2, 0, multiplier2,    1);
  module->addAudioConnection(multiplier1,    0, audioOutput1,   0);
  module->addAudioConnection(multiplier2,    0, audioOutput2,   0);

  return module;
}








/*
romos::Module* TestModuleBuilder::createVoiceKill(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule *module = (ContainerModule*) ModuleFactory::createModule(ModuleTypeRegistry::CONTAINER, name, x, y, polyphonic);
  
  romos::Module *audioInput1  = module->addChildModule(("AudioInput"),  "In1",          2,  2, true,  false);
  romos::Module *constant1    = module->addChildModule(("Constant"),    "0.001",        2,  5, true,  false);
  romos::Module *constant2    = module->addChildModule(("Constant"),    "0.01",         2,  8, true,  false);
  romos::Module *voiceKiller1 = module->addChildModule(("VoiceKiller"), "VoiceKiller", 12,  4, true,  false);

  module->sortChildModuleArray();

  module->addAudioConnection(audioInput1,  0, voiceKiller1, 0);
  module->addAudioConnection(constant1,    0, voiceKiller1, 1);
  module->addAudioConnection(constant2,    0, voiceKiller1, 2);

  return module;
}

romos::Module* TestModuleBuilder::create1In2Out(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule *module = (ContainerModule*) ModuleFactory::createModule(ModuleTypeRegistry::CONTAINER, name, x, y, polyphonic);
  
  romos::Module *audioInput1  = module->addChildModule(("AudioInput"),  "In1",   2,  4, true,  false);
  romos::Module *audioOutput1 = module->addChildModule(("AudioOutput"), "Out1", 11,  2, true,  false);
  romos::Module *audioOutput2 = module->addChildModule(("AudioOutput"), "Out2", 11,  6, true,  false);

  module->sortChildModuleArray();

  module->addAudioConnection(audioInput1,  0, audioOutput1, 0);
  module->addAudioConnection(audioInput1,  0, audioOutput2, 0);

  return module;
}

romos::Module* TestModuleBuilder::createGatedNoise(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule *module = (ContainerModule*) ModuleFactory::createModule(ModuleTypeRegistry::CONTAINER, name, x, y, polyphonic);
  
  romos::Module *constant1                   = module->addChildModule(("Constant"),                   "0",           2,  4, true,  false);
  romos::Module *noiseGeneratorWhiteUniform1 = module->addChildModule(("WhiteNoise"), "WhiteNoise",  6,  4, true,  false);
  romos::Module *noteGate1                   = module->addChildModule(("NoteGate"),                   "NoteGate",    6,  8, true,  false);
  romos::Module *multiplier1                 = module->addChildModule(("Multiplier"),                 "*",          18,  4, true,  false);
  romos::Module *audioOutput1                = module->addChildModule(("AudioOutput"),                "Out1",       23,  4, true,  false);

  module->sortChildModuleArray();

  module->addAudioConnection(constant1,                   0, noiseGeneratorWhiteUniform1, 0);
  module->addAudioConnection(noiseGeneratorWhiteUniform1, 0, multiplier1,                 0);
  module->addAudioConnection(noteGate1,                   0, multiplier1,                 1);
  module->addAudioConnection(multiplier1,                 0, audioOutput1,                0);

  return module;
}

romos::Module* TestModuleBuilder::createGatedNoiseStereo(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule *module = (ContainerModule*) ModuleFactory::createModule(ModuleTypeRegistry::CONTAINER, name, x, y, polyphonic);
  
  romos::Module *voiceCombiner1 = module->addChildModule(("VoiceCombiner"), "}",         23,  4, false, false);
  romos::Module *constant1      = module->addChildModule(("Constant"),      "0.125",     23,  7, false, false);
  romos::Module *multiplier1    = module->addChildModule(("Multiplier"),    "*",         29,  4, false, false);
  romos::Module *audioOutput1   = module->addChildModule(("AudioOutput"),   "Out1",      35,  4, false, false);

  romos::Module *gatedNoise = module->addChildModule(createGatedNoise("GatedNoise",  3,  4, true ));
  romos::Module *In1Out2    = module->addChildModule(create1In2Out(   "In1Out2",    14,  4, true ));
  romos::Module *voiceKill  = module->addChildModule(createVoiceKill( "VoiceKill",  14, 12, true ));

  module->sortChildModuleArray();

  module->addAudioConnection(gatedNoise,     0, In1Out2,        0);
  module->addAudioConnection(gatedNoise,     0, voiceKill,      0);
  module->addAudioConnection(In1Out2,        1, voiceCombiner1, 0);
  module->addAudioConnection(voiceCombiner1, 0, multiplier1,    0);
  module->addAudioConnection(constant1,      0, multiplier1,    1);
  module->addAudioConnection(multiplier1,    0, audioOutput1,   0);

  return module;
}
*/











// an input connected to a minus connected to two outputs
romos::Module* TestModuleBuilder::createContainerize01(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *audioInput1  = module->addChildModule("AudioInput",  "In1",   2,  3, false, false);
  romos::Module *unaryMinus1  = module->addChildModule("UnaryMinus",  "-",    11,  3, false, false);
  romos::Module *audioOutput1 = module->addChildModule("AudioOutput", "Out1", 19,  3, false, false);
  romos::Module *audioOutput2 = module->addChildModule("AudioOutput", "Out2", 19,  8, false, false);

  module->sortChildModuleArray();

  module->addAudioConnection(audioInput1,  0, unaryMinus1,  0);
  module->addAudioConnection(unaryMinus1,  0, audioOutput1, 0);
  module->addAudioConnection(unaryMinus1,  0, audioOutput2, 0);

  return module;
}

romos::Module* TestModuleBuilder::createContainerize02(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *constant1       = module->addChildModule(("Constant"),       "1",               2,  5, false, false);
  romos::Module *constant2       = module->addChildModule(("Constant"),       "2",               2, 15, false, false);
  romos::Module *constant3       = module->addChildModule(("Constant"),       "3",               7, 10, false, false);
  romos::Module *biquadDesigner1 = module->addChildModule(("BiquadDesigner"), "BiquadDesigner", 12,  7, false, false);
  romos::Module *adder1          = module->addChildModule(("Adder"),          "+",              31, 12, false, false);
  romos::Module *multiplier1     = module->addChildModule(("Multiplier"),     "*",              35, 11, false, false);
  romos::Module *audioOutput1    = module->addChildModule(("AudioOutput"),    "Out1",           39,  7, false, false);
  romos::Module *audioOutput2    = module->addChildModule(("AudioOutput"),    "Out2",           39, 10, false, false);
  romos::Module *audioOutput3    = module->addChildModule(("AudioOutput"),    "Out3",           39, 15, false, false);

  module->sortChildModuleArray();

  module->addAudioConnection(constant1,       0, biquadDesigner1, 0);
  module->addAudioConnection(constant3,       0, biquadDesigner1, 1);
  module->addAudioConnection(constant2,       0, biquadDesigner1, 2);
  module->addAudioConnection(biquadDesigner1, 3, adder1,          0);
  module->addAudioConnection(biquadDesigner1, 2, multiplier1,     0);
  module->addAudioConnection(biquadDesigner1, 0, audioOutput1,    0);
  module->addAudioConnection(biquadDesigner1, 1, audioOutput2,    0);
  module->addAudioConnection(constant2,       0, audioOutput3,    0);

  return module;
}

romos::Module* TestModuleBuilder::createOneTwoThree(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *constant1    = module->addChildModule(("Constant"),    "1",     2,  2, false, false);
  romos::Module *constant2    = module->addChildModule(("Constant"),    "2",     2,  5, false, false);
  romos::Module *constant3    = module->addChildModule(("Constant"),    "3",     2,  8, false, false);
  romos::Module *audioOutput1 = module->addChildModule(("AudioOutput"), "Out1",  7,  2, false, false);
  romos::Module *audioOutput2 = module->addChildModule(("AudioOutput"), "Out2",  7,  5, false, false);
  romos::Module *audioOutput3 = module->addChildModule(("AudioOutput"), "Out3",  7,  8, false, false);

  module->sortChildModuleArray();

  module->addAudioConnection(constant1,    0, audioOutput1, 0);
  module->addAudioConnection(constant2,    0, audioOutput2, 0);
  module->addAudioConnection(constant3,    0, audioOutput3, 0);

  return module;
}

romos::Module* TestModuleBuilder::createOutputDeletion(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *unaryMinus1 = module->addChildModule(("UnaryMinus"), "-", 16,  3, false, false);
  romos::Module *unaryMinus2 = module->addChildModule(("UnaryMinus"), "-", 16,  5, false, false);
  romos::Module *unaryMinus3 = module->addChildModule(("UnaryMinus"), "-", 16,  7, false, false);

  romos::Module *oneTwoThree = module->addChildModule(createOneTwoThree("OneTwoThree",  1,  4, false));

  module->sortChildModuleArray();

  module->addAudioConnection(oneTwoThree, 0, unaryMinus1, 0);
  module->addAudioConnection(oneTwoThree, 1, unaryMinus2, 0);
  module->addAudioConnection(oneTwoThree, 2, unaryMinus3, 0);

  return module;
}






romos::Module* TestModuleBuilder::createAddedConstants(const rosic::rsString &name, int x, int y, bool polyphonic)
{
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);
  
  romos::Module *constant1    = module->addChildModule(("Constant"),    "1",     2,  4, false, false);
  romos::Module *constant2    = module->addChildModule(("Constant"),    "3",     2,  6, false, false);
  romos::Module *constant3    = module->addChildModule(("Constant"),    "5",     2,  8, false, false);
  romos::Module *constant4    = module->addChildModule(("Constant"),    "7",     2, 10, false, false);
  romos::Module *constant5    = module->addChildModule(("Constant"),    "9",     2, 12, false, false);
  romos::Module *constant6    = module->addChildModule(("Constant"),    "11",    2, 14, false, false);
  romos::Module *constant7    = module->addChildModule(("Constant"),    "13",    2, 16, false, false);
  romos::Module *constant8    = module->addChildModule(("Constant"),    "15",    2, 18, false, false);
  romos::Module *constant9    = module->addChildModule(("Constant"),    "2",     6,  5, false, false);
  romos::Module *constant10   = module->addChildModule(("Constant"),    "4",     6,  7, false, false);
  romos::Module *constant11   = module->addChildModule(("Constant"),    "6",     6,  9, false, false);
  romos::Module *constant12   = module->addChildModule(("Constant"),    "8",     6, 11, false, false);
  romos::Module *constant13   = module->addChildModule(("Constant"),    "10",    6, 13, false, false);
  romos::Module *constant14   = module->addChildModule(("Constant"),    "12",    6, 15, false, false);
  romos::Module *constant15   = module->addChildModule(("Constant"),    "14",    6, 17, false, false);
  romos::Module *constant16   = module->addChildModule(("Constant"),    "16",    6, 19, false, false);
  romos::Module *adderN1      = module->addChildModule(("AdderN"),      "+",    13,  4, false, false);
  romos::Module *audioOutput1 = module->addChildModule(("AudioOutput"), "Out1", 20,  4, false, false);

  module->sortChildModuleArray();

  module->addAudioConnection(constant1,    0, adderN1,     0);
  module->addAudioConnection(constant9,    0, adderN1,     1);
  module->addAudioConnection(constant2,    0, adderN1,     2);
  module->addAudioConnection(constant10,   0, adderN1,     3);
  module->addAudioConnection(constant3,    0, adderN1,     4);
  module->addAudioConnection(constant11,   0, adderN1,     5);
  module->addAudioConnection(constant4,    0, adderN1,     6);
  module->addAudioConnection(constant12,   0, adderN1,     7);
  module->addAudioConnection(constant5,    0, adderN1,     8);
  module->addAudioConnection(constant13,   0, adderN1,     9);
  module->addAudioConnection(constant6,    0, adderN1,     10);
  module->addAudioConnection(constant14,   0, adderN1,     11);
  module->addAudioConnection(constant7,    0, adderN1,     12);
  module->addAudioConnection(constant15,   0, adderN1,     13);
  module->addAudioConnection(constant8,    0, adderN1,     14);
  module->addAudioConnection(constant16,   0, adderN1,     15);
  module->addAudioConnection(adderN1,      0, audioOutput1, 0);

  return module;
}






romos::Module* TestModuleBuilder::createIdentityChain(const rosic::rsString &name, int x, int y, bool polyphonic)
{  
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);

  static const int numIdentities = 20;
  romos::Module *source = module->addChildModule(("Identity"), "Identity", 0,  2, false, false);
  romos::Module *target;
  for(int i = 1; i < numIdentities; i++)
  {
    target = module->addChildModule(("Identity"), "Identity", i,  2, false, false);
    module->addAudioConnection(source, 0, target, 0);
    source = target;
  }

  return module;
}

romos::Module* TestModuleBuilder::createIdentityChainWithFeedback(const rosic::rsString &name, int x, int y, bool polyphonic)
{  
  ContainerModule *module = (ContainerModule*) createIdentityChain(name, x, y, polyphonic);

  romos::Module *source = module->getChildModule(19);
  romos::Module *target = module->getChildModule( 0);

  module->addAudioConnection(source, 0, target, 0);

  return module;
}



romos::Module* TestModuleBuilder::createAdderChain(const rosic::rsString &name, int x, int y, bool polyphonic)
{  
  ContainerModule* module = (ContainerModule*) moduleFactory.createModule("Container", 
    name.asStdString(), x, y, polyphonic);

  static const int numAdders = 20;
  romos::Module *source = module->addChildModule(("Adder"), "Adder", 0,  2, false, false);
  romos::Module *target;
  for(int i = 1; i < numAdders; i++)
  {
    target = module->addChildModule(("Adder"), "Adder", i,  2, false, false);
    module->addAudioConnection(source, 0, target, 0);
    module->addAudioConnection(source, 0, target, 1);
    source = target;
  }

  return module;
}

romos::Module* TestModuleBuilder::createAdderChainWithFeedback(const rosic::rsString &name, int x, int y, bool polyphonic)
{  
  ContainerModule *module = (ContainerModule*) createAdderChain(name, x, y, polyphonic);

  romos::Module *source = module->getChildModule(19);
  romos::Module *target = module->getChildModule( 0);

  module->addAudioConnection(source, 0, target, 0);
  module->addAudioConnection(source, 0, target, 1);

  return module;
}



/*
romos::Module* romos::createMoogFilterModule(ModuleProperties *propertiesToUse)
{
  ContainerModule *moogFilter = new ContainerModule(propertiesToUse);  
  moogFilter->setModuleName("MoogFilter");

  moogFilter->addAudioInputModule(rosic::rsString("In"), 2, 14, false);
  romos::Module *in  = moogFilter->getAudioInputModule(0);

  moogFilter->addAudioInputModule(rosic::rsString("c"), 2, 22, false);
  romos::Module *c  = moogFilter->getAudioInputModule(1);

  moogFilter->addAudioInputModule(rosic::rsString("k"), 2, 26, false);
  romos::Module *k  = moogFilter->getAudioInputModule(2);

  moogFilter->addAudioOutputModule(rosic::rsString("Out"), 54, 4, false);
  romos::Module *out = moogFilter->getAudioOutputModule(0);

  AdderModule *adder = new AdderModule(NULL);
  adder->setPositionXY(8, 13);
  moogFilter->addChildModule(adder);

  Module *leakInt1 = createLeakyIntegratorModule(NULL);
  leakInt1->setPositionXY(12, 14);
  leakInt1->setModuleName("LeakInt1");
  moogFilter->addChildModule(leakInt1);

  Module *leakInt2 = createLeakyIntegratorModule(NULL);
  leakInt2->setPositionXY(22, 14);
  leakInt2->setModuleName("LeakInt2");
  moogFilter->addChildModule(leakInt2);

  Module *leakInt3 = createLeakyIntegratorModule(NULL);
  leakInt3->setPositionXY(32, 14);
  leakInt3->setModuleName("LeakInt3");
  moogFilter->addChildModule(leakInt3);

  Module *leakInt4 = createLeakyIntegratorModule(NULL);
  leakInt4->setPositionXY(42, 14);
  leakInt4->setModuleName("LeakInt4");
  moogFilter->addChildModule(leakInt4);

  moogFilter->addAudioConnection(in,          0, adder,    1);
  moogFilter->addAudioConnection(adder,       0, leakInt1, 0);
  moogFilter->addAudioConnection(leakInt1,    0, leakInt2, 0);
  moogFilter->addAudioConnection(leakInt2,    0, leakInt3, 0);
  moogFilter->addAudioConnection(leakInt3,    0, leakInt4, 0);
  moogFilter->addAudioConnection(leakInt4,    0, out,      0);
  moogFilter->addAudioConnection(c,           0, leakInt1, 1);
  moogFilter->addAudioConnection(c,           0, leakInt2, 1);
  moogFilter->addAudioConnection(c,           0, leakInt3, 1);
  moogFilter->addAudioConnection(c,           0, leakInt4, 1);

  // ..not yet complete - feedback path is missing

  return moogFilter;
}
*/

/*
romos::Module* TestModuleBuilder::createContainerizeTestModule(ModuleProperties *propertiesToUse)
{
  ContainerModule *testModule = new ContainerModule("ContainerizeTest");  

  testModule->addAudioInputModule(rosic::rsString("In1"), 2, 2, false);
  //testModule->getAudioInputModule(0)->setPositionXY(2, 2);  
  romos::Module *in1  = testModule->getAudioInputModule(0);

  //...

  testModule->addAudioOutputModule(rosic::rsString("Out1"), 30, 2, false);
  //testModule->getAudioOutputModule(0)->setPositionXY(30, 2);
  romos::Module *out1 = testModule->getAudioOutputModule(0);

  //...

  ConstantModule *unity = new ConstantModule(NULL);
  unity->setModuleName(1.0);
  unity->setPositionXY(10, 2);
  testModule->addChildModule(unity);

  AdderModule *adder = new AdderModule(NULL);
  adder->setPositionXY(15, 2);
  testModule->addChildModule(adder);

  SubtractororModule *subtractoror = new SubtractororModule(NULL);
  subtractoror->setPositionXY(15, 6);
  testModule->addChildModule(subtractoror);

  MultiplierModule *multiplier = new MultiplierModule(NULL);
  multiplier->setPositionXY(20, 4);
  testModule->addChildModule(multiplier);


  testModule->addAudioConnection(unity,       0, adder,       0);
  testModule->addAudioConnection(unity,       0, subtractoror,  0);
  testModule->addAudioConnection(adder,       0, multiplier,  0);

  return testModule;
}

romos::Module* TestModuleBuilder::createUnContainerizeInnerModule(ModuleProperties *propertiesToUse)
{
  ContainerModule *testModule = new ContainerModule("InnerContainer");  

  testModule->addAudioInputModule(rosic::rsString("IIn1"), 2, 2, false);
  //testModule->getAudioInputModule(0)->setPositionXY(2, 2);  
  romos::Module *in1  = testModule->getAudioInputModule(0);

  testModule->addAudioInputModule(rosic::rsString("IIn2"), 2, 6, false);
  //testModule->getAudioInputModule(1)->setPositionXY(2, 6);  
  romos::Module *in2  = testModule->getAudioInputModule(1);

  testModule->addAudioOutputModule(rosic::rsString("IOut1"), 30, 2, false);
  //testModule->getAudioOutputModule(0)->setPositionXY(30, 2);
  romos::Module *out1 = testModule->getAudioOutputModule(0);

  testModule->addAudioOutputModule(rosic::rsString("IOut2"), 30, 6, false);
  //testModule->getAudioOutputModule(1)->setPositionXY(30, 6);
  romos::Module *out2 = testModule->getAudioOutputModule(1);

  UnitDelayModule *d1 = new UnitDelayModule(NULL);
  d1->setModuleName("D1");
  d1->setPositionXY(15, 2);
  testModule->addChildModule(d1);

  UnitDelayModule *d2 = new UnitDelayModule(NULL);
  d2->setModuleName("D2");
  d2->setPositionXY(15, 6);
  testModule->addChildModule(d2);

  UnitDelayModule *d3 = new UnitDelayModule(NULL);
  d3->setModuleName("D3");
  d3->setPositionXY(15, 10);
  testModule->addChildModule(d3);

  UnitDelayModule *d4 = new UnitDelayModule(NULL);
  d4->setModuleName("D4");
  d4->setPositionXY(15, 14);
  testModule->addChildModule(d4);

  testModule->addAudioConnection(in1, 0, d1,   0);
  testModule->addAudioConnection(in1, 0, d2,   0);
  testModule->addAudioConnection(in2, 0, d3,   0);
  testModule->addAudioConnection(d2,  0, out1, 0);
  testModule->addAudioConnection(d3,  0, out1, 0);
  testModule->addAudioConnection(d3,  0, out2, 0);

  return testModule;
}


romos::Module* TestModuleBuilder::createUnContainerizeTestModule(ModuleProperties *propertiesToUse)
{
  ContainerModule *testModule = new ContainerModule("OuterContainer");  


  testModule->addAudioInputModule(rosic::rsString("OIn1"), 2, 2, false);
  //testModule->getAudioInputModule(0)->setPositionXY(2, 2);  
  romos::Module *in1  = testModule->getAudioInputModule(0);

  testModule->addAudioInputModule(rosic::rsString("OIn2"), 2, 6, false);
  //testModule->getAudioInputModule(1)->setPositionXY(2, 6);  
  romos::Module *in2  = testModule->getAudioInputModule(1);

  testModule->addAudioInputModule(rosic::rsString("OIn3"), 2, 10, false);
  //testModule->getAudioInputModule(2)->setPositionXY(2, 10);  
  romos::Module *in3  = testModule->getAudioInputModule(2);

  testModule->addAudioOutputModule(rosic::rsString("OOut1"), 35, 2, false);
  //testModule->getAudioOutputModule(0)->setPositionXY(35, 2);
  romos::Module *out1 = testModule->getAudioOutputModule(0);

  testModule->addAudioOutputModule(rosic::rsString("OOut2"), 35, 6, false);
  //testModule->getAudioOutputModule(1)->setPositionXY(35, 6);
  romos::Module *out2 = testModule->getAudioOutputModule(1);

  Module *inner = createUnContainerizeInnerModule(NULL);
  inner->setPositionXY(15, 6);
  testModule->addChildModule(inner);

  //AdderModule *adder = new AdderModule(NULL);
  //adder->setPositionXY(15, 2);
  //testModule->addChildModule(adder);


  testModule->addAudioConnection(in1,       0, inner,       0);
  testModule->addAudioConnection(in2,       0, inner,       0);
  testModule->addAudioConnection(in3,       0, inner,       0);
  testModule->addAudioConnection(in2,       0, inner,       1);
  testModule->addAudioConnection(in3,       0, inner,       1);
  testModule->addAudioConnection(inner,     0, out1,        0);


  return testModule;
}

romos::Module* TestModuleBuilder::createMinimizeInputsInnerModule1(ModuleProperties *propertiesToUse)
{
  ContainerModule *testModule = new ContainerModule("InnerContainer");  

  testModule->addAudioInputModule(rosic::rsString("IIn1"), 2, 2, false);
  romos::Module *in1  = testModule->getAudioInputModule(0);

  testModule->addAudioInputModule(rosic::rsString("IIn2"), 2, 6, false);
  romos::Module *in2  = testModule->getAudioInputModule(1);

  testModule->addAudioOutputModule(rosic::rsString("IOut1"), 30, 2, false);
  romos::Module *out1 = testModule->getAudioOutputModule(0);

  testModule->addAudioOutputModule(rosic::rsString("IOut2"), 30, 6, false);
  romos::Module *out2 = testModule->getAudioOutputModule(1);

  UnitDelayModule *d1 = new UnitDelayModule(NULL);
  d1->setModuleName("D1");
  d1->setPositionXY(15, 2);
  testModule->addChildModule(d1);

  UnitDelayModule *d2 = new UnitDelayModule(NULL);
  d2->setModuleName("D2");
  d2->setPositionXY(15, 6);
  testModule->addChildModule(d2);


  testModule->addAudioConnection(in1, 0, d1,   0);
  testModule->addAudioConnection(in2, 0, d2,   0);
  testModule->addAudioConnection(d1,  0, out1, 0);
  testModule->addAudioConnection(d2,  0, out2, 0);

  return testModule;
}

romos::Module* TestModuleBuilder::createMinimizeInputsTestModule1(ModuleProperties *propertiesToUse)
{
  ContainerModule *testModule = new ContainerModule("MinimizeInputs1");  

  testModule->addAudioInputModule(rosic::rsString("OIn1"), 2, 2, false);
  romos::Module *in1  = testModule->getAudioInputModule(0);

  testModule->addAudioInputModule(rosic::rsString("OIn2"), 2, 6, false);
  romos::Module *in2  = testModule->getAudioInputModule(1);

  Module *inner = createMinimizeInputsInnerModule1(NULL);
  inner->setPositionXY(15, 6);
  testModule->addChildModule(inner);

  testModule->addAudioConnection(in1,       0, inner,       0);
  testModule->addAudioConnection(in1,       0, inner,       1);
  testModule->addAudioConnection(in2,       0, inner,       0);
  testModule->addAudioConnection(in2,       0, inner,       1);

  return testModule;
}

romos::Module* romos::createAudioPinSortingInnerModule(ModuleProperties *propertiesToUse)
{
  ContainerModule *testModule = new ContainerModule(propertiesToUse);  
  testModule->setModuleName("PinSortingInner");

  testModule->addAudioInputModule(rosic::rsString("1"), 2,  2, false);
  romos::Module *in1  = testModule->getAudioInputModule(0);

  testModule->addAudioInputModule(rosic::rsString("2"), 2,  6, false);
  romos::Module *in2  = testModule->getAudioInputModule(1);

  testModule->addAudioInputModule(rosic::rsString("3"), 2, 11, false);
  romos::Module *in3  = testModule->getAudioInputModule(2);

  testModule->addAudioOutputModule(rosic::rsString("1p2"), 30,  2, false);
  romos::Module *out1 = testModule->getAudioOutputModule(0);

  testModule->addAudioOutputModule(rosic::rsString("1p3"), 30,  6, false);
  romos::Module *out2 = testModule->getAudioOutputModule(1);

  testModule->addAudioOutputModule(rosic::rsString("2p3"), 30, 10, false);
  romos::Module *out3 = testModule->getAudioOutputModule(2);

  AdderModule *adder1 = new AdderModule(NULL);
  adder1->setPositionXY(15, 2);
  testModule->addChildModule(adder1);

  AdderModule *adder2 = new AdderModule(NULL);
  adder2->setPositionXY(15, 6);
  testModule->addChildModule(adder2);

  AdderModule *adder3 = new AdderModule(NULL);
  adder3->setPositionXY(15, 10);
  testModule->addChildModule(adder3);

  testModule->addAudioConnection(in1,    0, adder1, 0);
  testModule->addAudioConnection(in2,    0, adder1, 1);
  testModule->addAudioConnection(in1,    0, adder2, 0);
  testModule->addAudioConnection(in3,    0, adder2, 1);
  testModule->addAudioConnection(in2,    0, adder3, 0);
  testModule->addAudioConnection(in3,    0, adder3, 1);
  testModule->addAudioConnection(adder1, 0, out1,   0);
  testModule->addAudioConnection(adder2, 0, out2,   0);
  testModule->addAudioConnection(adder3, 0, out3,   0);

  return testModule;
}

romos::Module* romos::createAudioPinSortingTestModule(ModuleProperties *propertiesToUse)
{
  ContainerModule *testModule = new ContainerModule(propertiesToUse);  
  testModule->setModuleName("PinSortTest");

  testModule->addAudioInputModule(rosic::rsString("1"), 2,  2, false);
  romos::Module *in1  = testModule->getAudioInputModule(0);

  testModule->addAudioInputModule(rosic::rsString("2"), 2,  6, false);
  romos::Module *in2  = testModule->getAudioInputModule(1);

  testModule->addAudioInputModule(rosic::rsString("3"), 2, 10, false);
  romos::Module *in3  = testModule->getAudioInputModule(2);

  testModule->addAudioOutputModule(rosic::rsString("1p2"), 30,  2, false);
  romos::Module *out1 = testModule->getAudioOutputModule(0);

  testModule->addAudioOutputModule(rosic::rsString("1p3"), 30,  6, false);
  romos::Module *out2 = testModule->getAudioOutputModule(1);

  testModule->addAudioOutputModule(rosic::rsString("2p3"), 30, 10, false);
  romos::Module *out3 = testModule->getAudioOutputModule(2);

  romos::Module *inner = createAudioPinSortingInnerModule(NULL);
  inner->setPositionXY(12, 5);
  testModule->addChildModule(inner);

  testModule->addAudioConnection(in1,    0, inner, 0);
  testModule->addAudioConnection(in2,    0, inner, 1);
  testModule->addAudioConnection(in3,    0, inner, 2);
  testModule->addAudioConnection(inner,  0, out1,  0);
  testModule->addAudioConnection(inner,  1, out2,  0);
  testModule->addAudioConnection(inner,  2, out3,  0);

  return testModule;
}
*/

}