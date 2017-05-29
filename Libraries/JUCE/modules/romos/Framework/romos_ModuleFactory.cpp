#include "romos_ModuleFactory.h"
using namespace romos;

//-----------------------------------------------------------------------------------------------------------------------------------------
// module creation:
    
romos::Module* ModuleFactory::createModule(int typeIdentifier, rosic::rsString name, int x, int y, bool polyphonic)
{
  // use default-names, if the name parameter is empty:
  if( name.isEmpty() )
  {
    switch( typeIdentifier )
    {
      case ModuleTypeRegistry::CONSTANT:        name = "0.0"; break;
      case ModuleTypeRegistry::IDENTITY:        name = "";    break;
      case ModuleTypeRegistry::UNARY_MINUS:     name = "-";   break;
      case ModuleTypeRegistry::RECIPROCAL:      name = "1/x"; break;
      case ModuleTypeRegistry::ADDER:           name = "+";   break;
      case ModuleTypeRegistry::SUBTRACTOR:      name = "-";   break;
      case ModuleTypeRegistry::MULTIPLIER:      name = "*";   break;
      case ModuleTypeRegistry::DIVIDER:         name = "/";   break;
      case ModuleTypeRegistry::ADDER_3:         name = "+";   break;
      case ModuleTypeRegistry::ADDER_4:         name = "+";   break;
      case ModuleTypeRegistry::ADDER_5:         name = "+";   break;
      case ModuleTypeRegistry::ADDER_N:         name = "+";   break;
      case ModuleTypeRegistry::UNIT_DELAY:      name = "D";   break;
      case ModuleTypeRegistry::VOICE_COMBINER:  name = "}";   break;
      default: name = ModuleTypeRegistry::getSoleInstance()->getModuleTypeStringFromIdentifier(typeIdentifier);
    }
  }

  romos::Module *newModule = NULL;

  switch( typeIdentifier )
  {
  case ModuleTypeRegistry::TOP_LEVEL_MODULE:          newModule = new TopLevelModule();             break;
  case ModuleTypeRegistry::AUDIO_INPUT:               newModule = new AudioInputModule();           break;
  case ModuleTypeRegistry::AUDIO_OUTPUT:              newModule = new AudioOutputModule();          break;
  case ModuleTypeRegistry::PARAMETER:                 newModule = new ParameterModule();            break;
  case ModuleTypeRegistry::SYSTEM_SAMPLE_RATE:        newModule = new SystemSampleRateModule();     break;
  case ModuleTypeRegistry::NOTE_GATE:                 newModule = new NoteGateModule();             break;
  case ModuleTypeRegistry::NOTE_ON_TRIGGER:           newModule = new NoteOnTriggerModule();        break;
  case ModuleTypeRegistry::NOTE_OFF_TRIGGER:          newModule = new NoteOffTriggerModule();       break;
  case ModuleTypeRegistry::VOICE_KILLER:              newModule = new VoiceKillerModule();          break;
  case ModuleTypeRegistry::VOICE_COMBINER:            newModule = new VoiceCombinerModule();        break;
  case ModuleTypeRegistry::NOTE_FREQUENCY:            newModule = new NoteFrequencyModule();        break;
  case ModuleTypeRegistry::NOTE_VELOCITY:             newModule = new NoteVelocityModule();         break;
  case ModuleTypeRegistry::ADDER:                     newModule = new AdderModule();                break;
  case ModuleTypeRegistry::SUBTRACTOR:                newModule = new SubtractorModule();           break;
  case ModuleTypeRegistry::MULTIPLIER:                newModule = new MultiplierModule();           break;
  case ModuleTypeRegistry::DIVIDER:                   newModule = new DividerModule();              break;


  case ModuleTypeRegistry::ADDER_3:                   newModule = new Adder3Module();               break;
  case ModuleTypeRegistry::ADDER_4:                   newModule = new Adder4Module();               break;
  case ModuleTypeRegistry::ADDER_5:                   newModule = new Adder5Module();               break;
  case ModuleTypeRegistry::ADDER_N:                   newModule = new AdderNModule();               break;

  case ModuleTypeRegistry::CLIPPER:                   newModule = new ClipperModule();              break;
  case ModuleTypeRegistry::SIN_COS:                   newModule = new SinCosModule();               break;
  case ModuleTypeRegistry::CONSTANT:                  newModule = new ConstantModule();             break;
  case ModuleTypeRegistry::CONTAINER:                 newModule = new ModuleContainer();            break;


  case ModuleTypeRegistry::FIRST_ORDER_LOWPASS:       newModule = new FirstOrderLowpass();          break;
  case ModuleTypeRegistry::FIRST_ORDER_FILTER:        newModule = new FirstOrderFilter();           break;
  case ModuleTypeRegistry::BIQUAD:                    newModule = new Biquad();                     break;
  case ModuleTypeRegistry::BIQUAD_DESIGNER:           newModule = new BiquadDesigner();             break;
  case ModuleTypeRegistry::LADDER_FILTER:             newModule = new LadderFilter();               break;

  case ModuleTypeRegistry::IDENTITY:                  newModule = new IdentityModule();             break;
  case ModuleTypeRegistry::UNARY_MINUS:               newModule = new UnaryMinusModule();           break;
  case ModuleTypeRegistry::RECIPROCAL:                newModule = new ReciprocalModule();           break;
  case ModuleTypeRegistry::UNIT_DELAY:                newModule = new UnitDelayModule();            break;

  case ModuleTypeRegistry::PERIODIC_LINEAR_RAMP:      newModule = new PeriodicLinearRamp();         break;
  case ModuleTypeRegistry::WHITE_NOISE:               newModule = new WhiteNoise();                 break;
  case ModuleTypeRegistry::BANDLIMITED_IMPULSE_TRAIN: newModule = new BandlimitedImpulseTrain();    break;
  case ModuleTypeRegistry::BLIT_SAW_OSCILLATOR:       newModule = new BlitSaw();                    break;
  case ModuleTypeRegistry::DUAL_BLIT_SAW_OSCILLATOR:  newModule = new DualBlitSaw();                  break;



  case ModuleTypeRegistry::ENVELOPE_ADSR:             newModule = new EnvelopeADSR();               break;


  // test modules:
  case ModuleTypeRegistry::TEST_GAIN:              newModule = TestModuleBuilder::createGain(           name, x, y, polyphonic);  break;
  case ModuleTypeRegistry::TEST_SUM_DIFF:          newModule = TestModuleBuilder::createSumDiff(        name, x, y, polyphonic);  break;
  case ModuleTypeRegistry::TEST_WRAPPED_SUM_DIFF:  newModule = TestModuleBuilder::createWrappedSumDiff( name, x, y, polyphonic);  break;
  case ModuleTypeRegistry::TEST_SUMMED_DIFFS:      newModule = TestModuleBuilder::createSummedDiffs(    name, x, y, polyphonic);  break;
  case ModuleTypeRegistry::TEST_MOVING_AVERAGE:    newModule = TestModuleBuilder::createMovingAverage(  name, x, y, polyphonic);  break;
  case ModuleTypeRegistry::TEST_LEAKY_INTEGRATOR:  newModule = TestModuleBuilder::createLeakyIntegrator(name, x, y, polyphonic);  break;
  case ModuleTypeRegistry::TEST_FILTER1:           newModule = TestModuleBuilder::createTestFilter1(    name, x, y, polyphonic);  break;
  case ModuleTypeRegistry::TEST_BIQUAD:            newModule = TestModuleBuilder::createBiquadMacro(    name, x, y, polyphonic);  break;
  case ModuleTypeRegistry::TEST_ADDED_CONSTANTS:   newModule = TestModuleBuilder::createAddedConstants( name, x, y, polyphonic);  break;
  case ModuleTypeRegistry::TEST_PIN_SORTING:       newModule = TestModuleBuilder::createPinSortTest(    name, x, y, polyphonic);  break;
  case ModuleTypeRegistry::TEST_BLIP:              newModule = TestModuleBuilder::createBlip(           name, x, y, polyphonic);  break;
  case ModuleTypeRegistry::TEST_POLY_BLIP_STEREO:  newModule = TestModuleBuilder::createPolyBlipStereo( name, x, y, polyphonic);  break;
  case ModuleTypeRegistry::TEST_NOISE_FLUTE:       newModule = TestModuleBuilder::createNoiseFlute(     name, x, y, polyphonic);  break;


  default: 
    {
      DEBUG_BREAK;  // moduleIdentifier not known
      return NULL;
    }
  };

  // the test-containers, once created, should consider themselves as normal containers:
  if( typeIdentifier < ModuleTypeRegistry::NUM_NON_TEST_MODULE_TYPES )
  {
    //newModule->setModuleTypeIdentifier(typeIdentifier);
    newModule->moduleTypeIdentifier = typeIdentifier;
  }
  else
  {
    //newModule->setModuleTypeIdentifier(ModuleTypeRegistry::CONTAINER);
    newModule->moduleTypeIdentifier = ModuleTypeRegistry::CONTAINER;
  }

  newModule->initialize(); // this should set up the number of pins needed, etc.

  newModule->setPositionXY(x, y);
  newModule->setPolyphonic(polyphonic);
  newModule->allocateMemory();

  newModule->setModuleName(name); 
    // must be called after allocateMemory because the Constant fills its I/O arrays with the corresponding value

  newModule->assignProcessingFunctions();
  newModule->resetStateForAllVoices();

  return newModule;
}

void ModuleFactory::deleteModule(romos::Module *moduleToDelete)
{
  moduleToDelete->cleanUp();
  delete moduleToDelete;
}
