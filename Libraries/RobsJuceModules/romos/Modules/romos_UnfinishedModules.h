#ifndef romos_UnfinishedModules_h
#define romos_UnfinishedModules_h


// Under construction:
class PhasorPitchDithered : public AtomicModule
{
  CREATE_COMMON_DECLARATIONS_1(PhasorPitchDithered);   // Ins: Freq
public:
  virtual void resetVoiceState(int voiceIndex);
protected:
  virtual void allocateMemory();
  virtual void freeMemory();
  RAPT::rsPitchDitherOsc<double>* oscs = nullptr;
  bool rangeClosed = false;                            // ToDo: Make a GUI switch for that
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
