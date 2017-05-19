/*******************************************************************************
 The block below describes the properties of this module, and is read by
 the Projucer to automatically generate project code that uses it.
 For details about the syntax and how to create or use a module, see the
 JUCE Module Format.txt file.


 BEGIN_JUCE_MODULE_DECLARATION

  ID:               rosic
  vendor:           RS-MET
  version:          0.0.1
  name:             Rob's Signal Processing Classes
  description:      Library of audio DSP algorithms
  website:          http://www.rs-met.com
  license:          GPL/Commercial

  dependencies:     
  OSXFrameworks:
  iOSFrameworks:

 END_JUCE_MODULE_DECLARATION

*******************************************************************************/

#ifndef ROSIC_H_INCLUDED
#define ROSIC_H_INCLUDED

/** ToDo: re-order the includes in the order in which they depend on each other */

//-------------------------------------------------------------------------------------------------
// new, centralized header includes (when finsihed, delete old ones):

#include <malloc.h>  // for alloca - try to get rid..
//#include <math.h>

// basics:
#include "basics/GlobalDefinitions.h"
#include "basics/GlobalFunctions.h"
#include "basics/rosic_Constants.h"
#include "basics/rosic_ChannelMatrix2x2.h"
#include "basics/rosic_HelperFunctions.h"
#include "basics/rosic_NumberManipulations.h"
#include "basics/rosic_FunctionTemplates.h"
// we need to intersperse some includes from other directories before we can finish the includes 
// from basics (todo: fix the dependency structure/layering):
#include "infrastructure/rosic_MutexLock.h"
#include "math/rosic_CephesDeclarations.h"
#include "math/rosic_RealFunctionEvaluationAlgorithms.h"
#include "math/rosic_IntegerFunctions.h"
#include "math/rosic_ElementaryFunctionsReal.h"
#include "math/rosic_SpecialFunctionsReal.h"
#include "math/rosic_Complex.h"
#include "math/rosic_ComplexFunctions.h"
#include "scripting/rosic_ExpressionEvaluatorFunctions.h"
#include "scripting/rosic_ExpressionEvaluatorComplexFunctions.h"
#include "scripting/rosic_ExpressionEvaluator.h"
// ..now we can finish the includes from "basics":
#include "basics/rosic_Interpolator.h"                // needs ElementaryFunctionsReal
#include "basics/rosic_SampleBuffer.h"                // needs MutexLock
#include "basics/rosic_SamplePlaybackParameters.h"    // needs ElementaryFunctionsReal
#include "basics/rosic_TabulatedFunction.h"           // needs ExpressionEvaluator, Mutexlock
#include "basics/rosic_WarpedAllpassInterpolator.h"   // needs ElementaryFunctionsReal
#include "basics/rosic_WindowDesigner.h"              // needs SpecialFunctionsReal

// datastructures:
#include "datastructures/rosic_Array.h"
#include "datastructures/rosic_String.h"
#include "datastructures/rosic_KeyValueMap.h"
#include "datastructures/rosic_ExtensionsForSTL.h"

// math:
#include "math/rosic_Interpolation.h"
#include "math/rosic_LinearAlgebra.h"
#include "math/rosic_Matrix.h"
#include "math/rosic_Vector.h"
#include "math/rosic_MatrixVectorFunctions.h"
#include "math/rosic_PolynomialAlgorithms.h"
#include "math/rosic_PrimeNumbers.h"
//#include "math/rosic_PrimeArray.h"                  // not required to include
#include "math/rosic_Transformations.h"

// transforms:
#include "transforms/rosic_FourierTransformerRadix2.h"
#include "transforms/rosic_FourierTransformerBluestein.h"
#include "transforms/rosic_SpectralManipulator.h"
#include "transforms/rosic_WaveletTransforms.h"

// numerical:


// neural:


// filters:
#include "filters/rosic_BiquadDesigner.h"
#include "filters/rosic_AllpassChain.h"
#include "filters/rosic_BiquadBase.h"
#include "filters/rosic_FilterAnalyzer.h"
#include "filters/rosic_BiquadCascade.h"
#include "filters/rosic_BiquadMonoDF1.h"
#include "filters/rosic_BiquadStereoDF1.h"
#include "filters/rosic_BiquadStereoDF2.h"
#include "filters/rosic_CombFilter.h"
#include "filters/rosic_DampingFilter.h"
#include "filters/rosic_CombResonator.h"
#include "filters/rosic_ConvolverBruteForce.h"
#include "filters/rosic_ConvolverFFT.h"          // needs transforms folder
#include "filters/rosic_ConvolverPartitioned.h"
#include "filters/rosic_CookbookFilter.h"
#include "filters/rosic_DirectFormFilter.h"
#include "filters/rosic_TwoPoleFilter.h"
#include "filters/rosic_DualTwoPoleFilter.h"
#include "filters/rosic_EllipticQuarterBandFilter.h"
#include "filters/rosic_Equalizer.h"
#include "filters/rosic_EqualizerStereo.h"
#include "filters/rosic_FilterCoefficientConverter.h"
#include "filters/rosic_FiniteImpulseResponseDesigner.h"
#include "filters/rosic_FiniteImpulseResponseFilter.h"
#include "filters/rosic_FourPoleFilter.h"
#include "filters/rosic_OnePoleFilter.h"
#include "filters/rosic_OnePoleFilterStereo.h"
#include "filters/rosic_LadderFilter.h"
#include "filters/rosic_LeakyIntegrator.h"
#include "filters/rosic_PrototypeDesigner.h"
#include "filters/rosic_PoleZeroMapper.h"
#include "filters/rosic_InfiniteImpulseResponseDesigner.h"
#include "filters/rosic_EngineersFilter.h"
#include "filters/rosic_EllipticSubBandFilter.h"
#include "filters/rosic_EllipticSubBandFilterDirectForm.h"
#include "filters/rosic_LinkwitzRileyCrossover.h"
#include "filters/rosic_Crossover4Way.h"
#include "filters/rosic_LowpassHighpass.h"
#include "filters/rosic_LowpassHighpassStereo.h"
#include "filters/rosic_LpfHpfApf.h"
#include "filters/rosic_MultiModeFilter.h"
#include "filters/rosic_NyquistBlocker.h"
#include "filters/rosic_QuadratureNetwork.h"
#include "filters/rosic_SlopeFilter.h"
#include "filters/rosic_TeeBeeFilter.h"
#include "filters/rosic_ToneControl.h"
#include "filters/rosic_TwoPoleBandpass.h"
#include "filters/rosic_VowelFilterStereo.h"
#include "filters/rosic_WarpedBiquadMonoDF1.h"
#include "filters/rosic_WhiteToPinkFilter.h"

// others:
//#include "others/rosic_DemoVersionNoiseEmitter.h"  // not needed anymore
#include "others/rosic_ExponentialSmoother.h"
#include "others/rosic_OverlapAddProcessor.h"
#include "others/rosic_PiecewiseFunction.h"
#include "_third_party/soundtouch/WavFile.h"
#include "infrastructure/rosic_FileInputOutput.h"
#include "others/rosic_Plotter.h"
#include "others/rosic_ProcessorCycleCounter.h"
#include "others/rosic_RandomNumberGenerator01.h"
#include "others/rosic_RandomNumberGenerator02.h"
#include "others/rosic_RoutingMatrix.h"
#include "others/rosic_RandomNumberGenerator01.h"
#include "others/rosic_SlewRateLimiter.h"
#include "others/rosic_SlewRateLimiterLinear.h"
//#include "others/rosic_SlewRateLimiterOld.h"
#include "others/rosic_SpectralProcessor.h"
#include "others/rosic_SpectralEnvelopeProcessor.h"
#include "others/rosic_SpectralFilter.h"
#include "others/rosic_TuningTable.h"
#include "others/rosic_VectorMixer.h"

// (most) generators:
#include "generators/rosic_MipMappedWaveTable.h"
//#include "generators/rosic_MipMappedWaveTableOld.h"
#include "generators/rosic_MipMappedWaveTableStereo.h"
#include "generators/rosic_BlendOscillator.h"
#include "generators/rosic_LorentzSystem.h"
#include "generators/rosic_ModalSynthesizer.h"
#include "generators/rosic_NoiseGenerator.h"
//#include "generators/rosic_NoiseGeneratorOld.h"
//#include "generators/rosic_Oscillator.h"
//#include "generators/rosic_OscillatorBank.h"
#include "generators/rosic_OscillatorStereo.h"
#include "generators/rosic_FourOscSection.h"
#include "generators/rosic_SampleOscillator.h"
#include "generators/rosic_SamplePlayer.h"
#include "generators/rosic_SineOscillator.h"
#include "generators/rosic_SineOscillatorStereo.h"
#include "generators/rosic_SuperOscillator.h"
#include "generators/rosic_TestGenerator.h"

// modulators:
//#include "generators/rosic_MagicCarpetModulator.h"
#include "modulators/rosic_AmpEnvRc.h"
#include "modulators/rosic_AnalogEnvelope.h"
#include "modulators/rosic_AnalogEnvelopeScaled.h"
#include "modulators/rosic_AttackDecayEnvelope.h"
#include "modulators/rosic_BreakpointModulator.h"
#include "modulators/rosic_DecayEnvelope.h"
#include "modulators/rosic_EnvelopeGenerator.h"
#include "modulators/rosic_EnvelopeGenerator2.h"
#include "modulators/rosic_ExponentialRamp.h"
#include "modulators/rosic_PitchEnvRc.h"
#include "modulators/rosic_SampleModulator.h"  
// todo: remove redundant envelope generators

// rendering:
#include "rendering/rosic_AlgorithmicWaveformRenderer.h"
#include "rendering/rosic_MultiSegmentWaveformRenderer.h"
#include "rendering/rosic_StandardWaveformRenderer.h"
#include "rendering/rosic_WaveformBuffer.h"
#include "rendering/rosic_WaveformRenderer.h"
#include "rendering/rosic_NonRealtimeProcesses.h"

// delaylines:
#include "delaylines/rosic_BasicIntegerDelayLine.h"
#include "delaylines/rosic_IntegerDelayLine.h"
#include "delaylines/rosic_AllpassDiffusor.h"
#include "delaylines/rosic_DelayLineStereo.h"
#include "delaylines/rosic_EchoLabDelayLine.h"
#include "delaylines/rosic_FractionalDelayLine.h"
#include "delaylines/rosic_FractionalDelayLineStereo.h"
#include "delaylines/rosic_ModulatedDelayLine.h"
#include "delaylines/rosic_PingPongEcho.h"             // needs StereoPan...or not?


// analysis:
#include "analysis/rosic_CyclicAutoCorrelator.h"
#include "analysis/rosic_EnvelopeFollower.h"
#include "analysis/rosic_LinearPredictor.h"
#include "analysis/rosic_FormantRemover.h"
#include "analysis/rosic_FormantPreserver.h"
#include "analysis/rosic_InstantaneousEnvelopeDetector.h"
#include "analysis/rosic_LevelDetector.h"
//#include "analysis/rosic_OscilloscopeBufferOld.h"
#include "analysis/rosic_PitchDetector.h"
#include "analysis/rosic_ScopeScreenScanner.h"
#include "analysis/rosic_SignalMeasures.h"
#include "analysis/rosic_SpectrumAnalyzer.h"
#include "analysis/rosic_TrackMeter.h"
#include "analysis/rosic_WaveformDisplayBuffer.h"

// dynamics:
#include "dynamics/rosic_BesselFilterForGainSignal.h"
#include "dynamics/rosic_DynamicsProcessorBase.h"
#include "dynamics/rosic_Compressor.h"
#include "dynamics/rosic_Expander.h"
#include "dynamics/rosic_Limiter.h"
#include "dynamics/rosic_NoiseGate.h"
#include "dynamics/rosic_SoftKneeCompressor.h"
#include "dynamics/rosic_SoftKneeExpander.h"

// effects:
#include "effects/rosic_FeedbackDelayNetwork.h"
#include "effects/rosic_FeedbackDelayNetwork8.h"
#include "effects/rosic_FeedbackDelayNetwork16.h"
#include "effects/rosic_AlgoVerb.h"
#include "effects/rosic_AudioToMidi.h"
#include "effects/rosic_BitCrusher.h"
#include "effects/rosic_ModulationEffect.h"
#include "effects/rosic_Vibrato.h"
#include "effects/rosic_Chorus.h"
#include "effects/rosic_CombBank.h"
#include "effects/rosic_CombResonatorStereo.h"
#include "effects/rosic_CombStereoizer.h"
#include "effects/rosic_CompShaper.h"
#include "effects/rosic_Distortion.h"
#include "effects/rosic_EchoLab.h"
#include "effects/rosic_Flanger.h"
#include "effects/rosic_FormantShifter.h"
#include "effects/rosic_FrequencyShifter.h"
#include "effects/rosic_FuncShaper.h"
#include "effects/rosic_Harmonics.h"
#include "effects/rosic_ModulatedAllpass.h"
#include "effects/rosic_Moduluxury.h"
#include "effects/rosic_Noisifier.h"
#include "effects/rosic_Phaser.h"
#include "effects/rosic_DelayPhaser.h"           // needs PingPongEcho and Phaser
#include "effects/rosic_PhaseStereoizer.h"
#include "effects/rosic_PitchShifter.h"
#include "effects/rosic_PitchShifterGrainAdaptive.h"
#include "effects/rosic_Reverb.h"
#include "effects/rosic_RingModulator.h"
#include "effects/rosic_SingleSidebandModulator.h"
#include "effects/rosic_StereoDelay.h"
#include "effects/rosic_StereoPan.h"
#include "effects/rosic_StereoWidth.h"
#include "effects/rosic_Tremolo.h"
#include "effects/rosic_WahWah.h"
#include "effects/rosic_WaveShaper.h"

// infrastructure:
#include "infrastructure/rosic_Module.h"
#include "infrastructure/rosic_EffectModules.h"
#include "infrastructure/rosic_GeneratorModules.h"
#include "infrastructure/rosic_ModulatorModules.h"
#include "infrastructure/rosic_File.h"
#include "infrastructure/rosic_MemoryUser.h"
#include "infrastructure/rosic_MidiNoteEvent.h"
#include "infrastructure/rosic_PolyphonicinstrumentVoice.h"
#include "infrastructure/rosic_Polyphonicinstrument.h"

// sequencing:
#include "sequencing/rosic_AcidPattern.h"
#include "sequencing/rosic_AcidSequencer.h"

// scripting:
#include "scripting/rosic_DspScriptInterpreter.h"
#include "scripting/rosic_DspWorkbench.h"

// some more complex generators that need includes from modulators and rendering:
#include "generators/rosic_WaveTable.h"              // needs WaveformRenderer
#include "modulators/rosic_LowFrequencyOscillator.h" // needs WaveTable
#include "effects/rosic_Quadrifex.h"                 // needs EffectModules
#include "generators/rosic_Quadrigen.h"              // needs BreakpointModulator, infrastructure/*Modules
#include "generators/rosic_VectorSamplePlayer.h"     // needs LowFrequencyOscillator
#include "modulators/rosic_Modulator.h"              // needs LowFrequencyOscillator

// ... good, until here


// instruments:


//-------------------------------------------------------------------------------------------------
// old, de-centralized header includes:

//#include "analysis/rosic_Analysis.h"
//#include "basics/rosic_Basics.h"
//#include "datastructures/rosic_DataStructures.h"
//#include "delaylines/rosic_DelayLines.h"
//#include "dynamics/rosic_Dynamics.h"
//#include "effects/rosic_Effects.h"
//#include "filters/rosic_Filters.h"
//#include "generators/rosic_Generators.h"
//#include "infrastructure/rosic_Infrastructure.h"
#include "instruments/rosic_Instruments.h"
//#include "math/rosic_Math.h"
//#include "modulators/rosic_Modulators.h"
#include "neural/rosic_Neural.h"
#include "numerical/rosic_Numerical.h"
//#include "others/rosic_Others.h"
//#include "rendering/rosic_Rendering.h"
//#include "scripting/rosic_Scripting.h"
//#include "plugins/rosic_PlugIns.h"
//#include "transforms/rosic_Transforms.h"

#endif // #ifndef ROSIC_H_INCLUDED







