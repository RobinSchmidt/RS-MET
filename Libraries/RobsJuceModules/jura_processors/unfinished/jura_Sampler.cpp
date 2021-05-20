// preliminary - eventually, this should go into the rosic module:
//#include "../../rs_testing/Prototypes/SamplerEngine.cpp"

SamplerModule::SamplerModule(CriticalSection *lockToUse) : AudioModule(lockToUse)
{
  ScopedLock scopedLock(*lock);
  setModuleTypeName("Sampler");
  setModuleName("Sampler");
}

AudioModuleEditor* SamplerModule::createEditor(int type)
{
  return new jura::SamplerEditor(this);
}

void SamplerModule::processBlock(double **inOutBuffer, int numChannels, int numSamples)
{
  //jassert(numChannels == 2);
  //for(int n = 0; n < numSamples; n++)
  //  core->processFrame(&inOutBuffer[0][n], &inOutBuffer[1][n]);
}

void SamplerModule::processStereoFrame(double *dL, double *dR)
{
  //float fL = (float) (*dL);
  //float fR = (float) (*dR);
  //core->processFrame(&fL, &fR);
  //*dL = fL;
  //*dR = fR;
}

void SamplerModule::setSampleRate(double newSampleRate)
{ 
  //core->setSampleRate(newSampleRate); 
}

void SamplerModule::reset()
{ 
  //core->reset();
}

//=================================================================================================

SamplerEditor::SamplerEditor(SamplerModule* SamplerToEdit)
  : AudioModuleEditor(SamplerToEdit->lock) //, SamplerModule(SamplerToEdit)
{
  ScopedLock scopedLock(*lock);

  // ...

  setSize(400, 200);
}

void SamplerEditor::resized()
{
  ScopedLock scopedLock(*lock);
  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom()+4;
  int w = getWidth();
  int h = getHeight() - y;

  // ...
}

/*
Maybe the xml presets for the SamplerModule should contain the filename for an sfz file that sould 
be loaded plus some global settings such as the maximum number of layers, the resampling algo, etc. 
Don't use the xml preset to store the actual sfz opcodes. Preset loading would become a nested 
process the outer level is the xml and the inner level is the sfz. The GUI should provide an 
additional set of load/save widgets for the sfz.
*/