#ifndef romos_UnfinishedModules_h
#define romos_UnfinishedModules_h


// Under construction:
class PhasorPitchDithered : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_3(PhasorPitchDithered);   // Ins: Freq, Min, Max
public:
  virtual void resetVoiceState(int voiceIndex);
protected:
  inline static void updatePhase(PhasorPitchDithered* phasor, double freq, int voiceIndex);
  virtual void allocateMemory();
  virtual void freeMemory();

  double *phases;
  //friend class SineOscillator;
  // ToDo: Maybe embedd a RAPT::rsPhasorPitchDithered<double> object. Or maybe we need an array of
  // them? 

};
class PhasorPitchDitheredTypeInfo : public ModuleTypeInfo
{
public:
  PhasorPitchDitheredTypeInfo() {
    shortName    = "PhsrPiDi";
    fullName     = "PhasorPitchDithered";
    description  = "Pitch dithered phasor with given frequency";
    category     = "Sources";
    createModule =  []()->Module* { return new PhasorPitchDithered; };
    hasHeader = false;
  }
};



#endif 
