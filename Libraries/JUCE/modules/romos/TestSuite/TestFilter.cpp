// code seems to be obsolete and superseded by stuff in class TestModuleBuilder
/*
romos::Module* romos::createLeakyIntegrator(int x, int y, const rosic::rsString &name)
{
  ModuleContainer *module = createModuleContainer(x, y, name);
  
  romos::Module *constant1    = module->addChildModule(getTypeId("Constant"),     2,  5, "1",   false);
  romos::Module *audioInput1  = module->addChildModule(getTypeId("AudioInput"),   2,  8, "In",  false);
  romos::Module *audioInput2  = module->addChildModule(getTypeId("AudioInput"),   2, 11, "c",   false);
  romos::Module *subtract1    = module->addChildModule(getTypeId("Subtract"),     8,  6, "-",   false);
  romos::Module *multiply1    = module->addChildModule(getTypeId("Multiply"),    12,  5, "*",   false);
  romos::Module *multiply2    = module->addChildModule(getTypeId("Multiply"),    12,  9, "*",   false);
  romos::Module *add1         = module->addChildModule(getTypeId("Add"),         18,  7, "+",   false);
  romos::Module *identity1    = module->addChildModule(getTypeId("Identity"),    22,  2, "",    false);
  romos::Module *audioOutput1 = module->addChildModule(getTypeId("AudioOutput"), 28,  7, "Out", false);

  module->addAudioConnection(constant1,    0, subtract1,    0);
  module->addAudioConnection(audioInput2,  0, subtract1,    1);
  module->addAudioConnection(subtract1,    0, multiply1,    1);
  module->addAudioConnection(identity1,    0, multiply1,    0);
  module->addAudioConnection(audioInput2,  0, multiply2,    0);
  module->addAudioConnection(audioInput1,  0, multiply2,    1);
  module->addAudioConnection(multiply1,    0, add1,         0);
  module->addAudioConnection(multiply2,    0, add1,         1);
  module->addAudioConnection(add1,         0, identity1,    0);
  module->addAudioConnection(add1,         0, audioOutput1, 0);

  return module;
}


romos::Module* romos::createMovingAverage(int x, int y, const rosic::rsString &name)
{
  ModuleContainer *module = createModuleContainer(x, y, name);
  
  romos::Module *audioInput1  = module->addChildModule(getTypeId("AudioInput"),   2,  2, "In",  false);
  romos::Module *audioInput2  = module->addChildModule(getTypeId("AudioInput"),   2,  5, "b0",  false);
  romos::Module *audioInput3  = module->addChildModule(getTypeId("AudioInput"),   2,  8, "b1",  false);
  romos::Module *unitDelay1   = module->addChildModule(getTypeId("UnitDelay"),    8,  7, "D",   false);
  romos::Module *multiply1    = module->addChildModule(getTypeId("Multiply"),    12,  2, "*",   false);
  romos::Module *multiply2    = module->addChildModule(getTypeId("Multiply"),    12,  7, "*",   false);
  romos::Module *add1         = module->addChildModule(getTypeId("Add"),         16,  4, "+",   false);
  romos::Module *audioOutput1 = module->addChildModule(getTypeId("AudioOutput"), 20,  4, "Out", false);

  module->addAudioConnection(audioInput1,  0, unitDelay1,   0);
  module->addAudioConnection(audioInput2,  0, multiply1,    1);
  module->addAudioConnection(audioInput1,  0, multiply1,    0);
  module->addAudioConnection(audioInput3,  0, multiply2,    1);
  module->addAudioConnection(unitDelay1,   0, multiply2,    0);
  module->addAudioConnection(multiply1,    0, add1,         0);
  module->addAudioConnection(multiply2,    0, add1,         1);
  module->addAudioConnection(add1,         0, audioOutput1, 0);

  return module;
}


romos::Module* romos::createTestFilter1(int x, int y, const rosic::rsString &name)
{
  ModuleContainer *module = createModuleContainer(x, y, name);
  
  romos::Module *audioInput1  = module->addChildModule(getTypeId("AudioInput"),   2,  2, "In",         false);
  romos::Module *audioInput2  = module->addChildModule(getTypeId("AudioInput"),   2,  6, "b0",         false);
  romos::Module *audioInput3  = module->addChildModule(getTypeId("AudioInput"),   2, 10, "b1",         false);
  romos::Module *audioInput4  = module->addChildModule(getTypeId("AudioInput"),   2, 14, "c",          false);
  romos::Module *add1         = module->addChildModule(getTypeId("Add"),         28,  4, "+",          false);
  romos::Module *subtract1    = module->addChildModule(getTypeId("Subtract"),    28,  8, "-",          false);
  romos::Module *multiply1    = module->addChildModule(getTypeId("Multiply"),    28, 12, "*",          false);
  romos::Module *audioOutput1 = module->addChildModule(getTypeId("AudioOutput"), 34,  4, "Sum",        false);
  romos::Module *audioOutput2 = module->addChildModule(getTypeId("AudioOutput"), 34,  8, "Difference", false);
  romos::Module *audioOutput3 = module->addChildModule(getTypeId("AudioOutput"), 34, 12, "Product",    false);

  romos::Module *movingAverage   = module->addChildModule(createMovingAverage(  10,  4, "MovingAverage"));
  romos::Module *leakyIntegrator = module->addChildModule(createLeakyIntegrator(10, 12, "LeakyIntegrator"));

  module->sortChildModuleArray();

  module->addAudioConnection(audioInput1,     0, movingAverage,   0);
  module->addAudioConnection(audioInput2,     0, movingAverage,   1);
  module->addAudioConnection(audioInput3,     0, movingAverage,   2);
  module->addAudioConnection(audioInput1,     0, leakyIntegrator, 0);
  module->addAudioConnection(audioInput4,     0, leakyIntegrator, 1);
  module->addAudioConnection(movingAverage,   0, add1,            0);
  module->addAudioConnection(leakyIntegrator, 0, add1,            1);
  module->addAudioConnection(movingAverage,   0, subtract1,       0);
  module->addAudioConnection(leakyIntegrator, 0, subtract1,       1);
  module->addAudioConnection(movingAverage,   0, multiply1,       0);
  module->addAudioConnection(leakyIntegrator, 0, multiply1,       1);
  module->addAudioConnection(add1,            0, audioOutput1,    0);
  module->addAudioConnection(subtract1,       0, audioOutput2,    0);
  module->addAudioConnection(multiply1,       0, audioOutput3,    0);

  return module;
}
*/