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
  license:          Custom

  dependencies:     rapt
  OSXFrameworks:
  iOSFrameworks:

 END_JUCE_MODULE_DECLARATION

*******************************************************************************/

#ifndef ROSIC_H_INCLUDED
#define ROSIC_H_INCLUDED

//-------------------------------------------------------------------------------------------------
// new, centralized header includes (when finished, delete old ones):

#include <rapt/rapt.h>

#include <math.h>
#include <float.h>      // for DBL_MIN, ..
//#include <algorithm>    // for sort
#include <cstring>      // for memcpy on linux/gcc
#include <list>
#include <stdlib.h>     // needed on OSX for atoi, atof, rand, abs
#include <iostream>     // cout is used in GradientBasedMinimizer.cpp (maybe get rid)
#include <climits>      // for INT_MAX on linux
//#include <string.h>
//#include <vector>
//#include <errno.h>
//#include <stddef.h>   // for NULL macro
//#include <functional>
//#include <limits.h>
#include <limits>
#include <atomic>
//#include <stdio.h>
//#include <new>
//#include <iterator>

#ifdef _MSC_VER           // Elans says, this is needed to compile on mac
#include <intrin.h>       // for some reason, it's needed on Elan's machine but not on mine
#endif

#include <emmintrin.h>    // SSE2

//#ifndef __APPLE__
//#include <malloc.h>     // for alloca - try to get rid - alloca is bad, mkay?
//#endif

//#if (defined (_WIN32) || defined (_WIN64))
//#define USE_WINAPI_CRITICAL_SECTION
//#elif (defined (MAC) || defined (LINUX))
//#define USE_PTHREADS_MUTEX
//#endif

// try to get rid of including windows.h - it has all sorts of weird #defines that lead to all
// sorts of compiler errors:
#if defined (USE_WINAPI_CRITICAL_SECTION)
#include <windows.h>
#elif defined (USE_PTHREADS_MUTEX)
#include <pthread.h>
#else
#include <mutex>
#endif

#if defined _MSC_VER
#include <crtdbg.h>
//#pragma warning( disable : 4996 ) // disables 'strcpy' may be unsafe warning in MSVC-2008
//#include <intrin.h> // Or #include <ia32intrin.h> etc.
#endif

// in msvc, set warning level to 3:
#if defined _MSC_VER
#pragma warning(push, 3)
#endif


// basics:
#include "basics/Float32x4.h"
#include "basics/Float64x2.h"
#include "basics/rosic_TemplateInstantiations.h"
#include "basics/GlobalDefinitions.h"
#include "basics/rosic_Constants.h"   // merge with rapt

//#ifdef _MSC_VER                       // i can't figure out, how to make this compile on mac
#include "basics/ComplexFloat64x2.h"
//#endif

#include "basics/rosic_ChannelMatrix2x2.h"
#include "basics/rosic_HelperFunctions.h"       // move to rapt
#include "basics/rosic_NumberManipulations.h"
// we need to intersperse some includes from other directories before we can finish the includes
// from basics (todo: fix the dependency structure/layering):
#include "infrastructure/rosic_MutexLock.h"
#include "math/rosic_CephesDeclarations.h"
#include "math/rosic_RealFunctionEvaluationAlgorithms.h"  // move to rapt
//#include "math/rosic_IntegerFunctions.h"                  // move to rapt
#include "math/rosic_SpecialFunctionsReal.h"              // move to rapt
#include "math/rosic_Complex.h"
#include "math/rosic_ComplexFunctions.h"
#include "_third_party/ExprEval_v3_4/expreval.h"
#include "scripting/rosic_ExpressionEvaluatorFunctions.h"
#include "scripting/rosic_ExpressionEvaluatorComplexFunctions.h"
#include "scripting/rosic_ExpressionEvaluator.h"
// ..now we can finish the includes from "basics":
#include "basics/rosic_Interpolator.h"                // needs ElementaryFunctionsReal, move to rapt
#include "basics/rosic_SampleBuffer.h"                // needs MutexLock
#include "basics/rosic_SamplePlaybackParameters.h"    // needs ElementaryFunctionsReal
#include "basics/rosic_TabulatedFunction.h"           // needs ExpressionEvaluator, Mutexlock
#include "basics/rosic_WarpedAllpassInterpolator.h"   // needs ElementaryFunctionsReal, move to rapt
#include "basics/rosic_WindowDesigner.h"              // needs SpecialFunctionsReal, move to rapt

// datastructures:
#include "datastructures/rosic_Array.h"
#include "datastructures/rosic_String.h"
#include "datastructures/rosic_KeyValueMap.h"
#include "datastructures/rosic_ExtensionsForSTL.h" // // todo: merge with RAPT StandardContainerFunctions

// math:
#include "math/rosic_LinearAlgebra.h"                 // move to / merge with rapt
#include "math/rosic_Matrix.h"                        // move to / merge with rapt
#include "math/rosic_Vector.h"                        // move to / merge with rapt
#include "math/rosic_MatrixVectorFunctions.h"         // move to / merge with rapt
#include "math/rosic_PrimeNumbers.h"                  // move to / merge with rapt
//#include "math/rosic_PrimeArray.h"                  // not required to include
#include "math/rosic_Transformations.h"               // move to / merge with rapt

// transforms:
#include "transforms/rosic_FourierTransformerRadix2.h"      // replace by rapt version
#include "transforms/rosic_FourierTransformerBluestein.h"   // replace by rapt version
#include "transforms/rosic_SpectralManipulator.h"           // move to rapt
#include "transforms/rosic_WaveletTransforms.h"             // move to rapt

// numerical:
#include "numerical/rosic_FunctionObjects.h"
#include "numerical/rosic_GradientBasedMinimizer.h"         // replace by rapt

// neural:
//#include "neural/rosic_MultiLayerPerceptron.h"            // superseded by rapt
//#include "neural/rosic_MultiLayerPerceptronErrorFunction.h"
//#include "neural/rosic_MultiLayerPerceptronTrainer.h"

// filters:
#include "filters/rosic_BiquadDesigner.h"
#include "filters/rosic_AllpassChain.h"
#include "filters/rosic_BiquadBase.h"
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
#include "filters/rosic_TwoPoleFilter.h"
#include "filters/rosic_DualTwoPoleFilter.h"
#include "filters/rosic_EllipticQuarterBandFilter.h"
#include "filters/rosic_Equalizer.h"
#include "filters/rosic_EqualizerStereo.h"
#include "filters/rosic_FiniteImpulseResponseDesigner.h"
#include "filters/rosic_FiniteImpulseResponseFilter.h"
#include "filters/rosic_FourPoleFilter.h"
#include "filters/rosic_LadderFilter.h"
#include "filters/rosic_LadderFilterOld.h"
#include "filters/rosic_LeakyIntegrator.h"
#include "filters/rosic_LowpassHighpass.h"
#include "filters/rosic_LowpassHighpassStereo.h"
#include "filters/rosic_LpfHpfApf.h"
#include "filters/rosic_MultiModeFilter.h"
#include "filters/rosic_NyquistBlocker.h"
#include "filters/rosic_SlopeFilter.h"
#include "filters/rosic_TeeBeeFilter.h"
#include "filters/rosic_ToneControl.h"
#include "filters/rosic_TwoPoleBandpass.h"
#include "filters/rosic_VowelFilterStereo.h"
#include "filters/rosic_WhiteToPinkFilter.h"

// others:
//#include "others/rosic_DemoVersionNoiseEmitter.h"  // not needed anymore
#include "others/rosic_ExponentialSmoother.h"
#include "others/rosic_OverlapAddProcessor.h"
#include "others/rosic_PiecewiseFunction.h"
#include "_third_party/soundtouch/WavFile.h"
#include "infrastructure/rosic_FileInputOutput.h"
//#include "others/rosic_Plotter.h"                  // obsolete - use GNUPlotCPP
#include "others/rosic_ProcessorCycleCounter.h"    // obsolete use PerformanceCounterTSC/QPC/..oh - but it's used by romos...hmmm
#include "others/rosic_RandomNumberGenerator01.h"
#include "others/rosic_RandomNumberGenerator02.h"
#include "others/rosic_RoutingMatrix.h"
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
#include "generators/rosic_MipMappedWaveTableOld.h"   // still needed?
#include "generators/rosic_MipMappedWaveTableStereo.h"
#include "generators/rosic_BlendOscillator.h"
#include "generators/rosic_LorentzSystem.h"
#include "generators/rosic_NoiseGenerator.h"
#include "generators/rosic_NoiseGeneratorOld.h"
#include "generators/rosic_Oscillator.h"
#include "generators/rosic_OscillatorBank.h"
#include "generators/rosic_OscillatorStereo.h"
#include "generators/rosic_FourOscSection.h"
#include "generators/rosic_SampleOscillator.h"
#include "generators/rosic_SamplePlayer.h"
#include "generators/rosic_SineOscillator.h"         // move to rapt
#include "generators/rosic_SineOscillatorStereo.h"   // move to rapt
#include "generators/rosic_SuperOscillator.h"
#include "generators/rosic_TestGenerator.h"

// modulators:
//#include "generators/rosic_MagicCarpetModulator.h"
#include "modulators/rosic_AmpEnvRc.h"
#include "modulators/rosic_AnalogEnvelope.h"
#include "modulators/rosic_AnalogEnvelopeScaled.h"
#include "modulators/rosic_VariousModulators.h"
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
#include "rendering/rosic_TurtleGraphics.h"
#include "rendering/rosic_LindenmayerSystem.h"

// some more complex generators/modulators that need includes from rendering:
#include "generators/rosic_WaveTable.h"              // needs WaveformRenderer
#include "generators/rosic_TurtleSource.h"           // needs TurtleGraphics
#include "generators/rosic_Snowflake.h"              // needs TurtleSource and LindenmayerSystem
#include "modulators/rosic_LowFrequencyOscillator.h" // needs WaveTable
#include "generators/rosic_VectorSamplePlayer.h"     // needs LowFrequencyOscillator
#include "modulators/rosic_Modulator.h"              // needs LowFrequencyOscillator

// delaylines:
#include "delaylines/rosic_BasicIntegerDelayLine.h"
#include "delaylines/rosic_IntegerDelayLine.h"
#include "delaylines/rosic_AllpassDiffusor.h"
#include "delaylines/rosic_DelayLineStereo.h"
#include "delaylines/rosic_EchoLabDelayLine.h"
#include "delaylines/rosic_FractionalDelayLine.h"
#include "delaylines/rosic_FractionalDelayLineStereo.h"
#include "delaylines/rosic_ModulatedDelayLine.h"
#include "delaylines/rosic_PingPongEcho.h"

// analysis:
#include "analysis/rosic_CyclicAutoCorrelator.h"
#include "analysis/rosic_EnvelopeFollower.h"
#include "analysis/rosic_LinearPredictor.h"
#include "analysis/rosic_FormantRemover.h"
#include "analysis/rosic_FormantPreserver.h"
#include "analysis/rosic_InstantaneousEnvelopeDetector.h"
#include "analysis/rosic_LevelDetector.h"
#include "analysis/rosic_OscilloscopeBufferOld.h"
#include "analysis/rosic_PitchDetector.h"
#include "analysis/rosic_ScopeScreenScanner.h"
#include "analysis/rosic_SignalMeasures.h"
#include "analysis/rosic_SpectrumAnalyzer.h"
#include "analysis/rosic_TrackMeter.h"
#include "analysis/rosic_WaveformDisplayBuffer.h"
#include "analysis/rosic_OnsetDetector.h"
#include "analysis/rosic_BeatDetector.h"

// dynamics:
#include "dynamics/rosic_BesselFilterForGainSignal.h"
#include "dynamics/rosic_DynamicsProcessorBase.h"
#include "dynamics/rosic_Compressor.h"
#include "dynamics/rosic_Expander.h"
#include "dynamics/rosic_Limiter.h"
#include "dynamics/rosic_NoiseGate.h"
#include "dynamics/rosic_SoftKneeCompressor.h"
#include "dynamics/rosic_SoftKneeExpander.h"
#include "dynamics/rosic_MultiBandCompressor.h"

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
#include "infrastructure/rosic_ModulationRouting.h"
#include "infrastructure/rosic_Module.h"
#include "infrastructure/rosic_EffectModules.h"
#include "infrastructure/rosic_GeneratorModules.h"
#include "infrastructure/rosic_ModulatorModules.h"
#include "infrastructure/rosic_File.h"
#include "infrastructure/rosic_MemoryUser.h"
#include "infrastructure/rosic_MidiNoteEvent.h"
#include "infrastructure/rosic_PolyphonicInstrumentVoice.h"
#include "infrastructure/rosic_PolyphonicInstrument.h"

// higher level generators/effects that neeed includes from infrastructure:
#include "effects/rosic_Quadrifex.h"           // needs EffectModules
#include "generators/rosic_Quadrigen.h"        // needs infrastructure/*Modules

// sequencing:
#include "sequencing/rosic_AcidPattern.h"
#include "sequencing/rosic_AcidSequencer.h"

// scripting:
#include "scripting/rosic_DspScriptInterpreter.h"
#include "scripting/rosic_DspWorkbench.h"

// instruments:
#include "instruments/rosic_AciDevil.h"
#include "instruments/rosic_KeyShotVoice.h"
#include "instruments/rosic_KeyShot.h"
#include "instruments/rosic_MagicCarpetVoice.h"
#include "instruments/rosic_MagicCarpet.h"
#include "instruments/rosic_Open303.h"
#include "instruments/rosic_QuadrigaVoice.h"
#include "instruments/rosic_Quadriga.h"
#include "instruments/rosic_SimpleSamplerOscSection.h"
#include "instruments/rosic_SimpleSamplerVoice.h"
#include "instruments/rosic_SimpleSampler.h"
#include "instruments/rosic_StraightlinerVoice.h"
#include "instruments/rosic_Straightliner.h"
#include "instruments/rosic_WorkhorseOscSection.h"
#include "instruments/rosic_WorkhorseVoice.h"
#include "instruments/rosic_Workhorse.h"

// unfinished:
#include "unfinished/rosic_Polyphony.h"
#include "unfinished/rosic_EllipseOscillator.h"
#include "unfinished/rosic_MiscUnfinished.h"
#include "unfinished/rosic_QuadSource.h"
#include "unfinished/rosic_DualFilter.h"
#include "unfinished/rosic_NewSynth.h"

#include "unfinished/rosic_ModalFilters.h"
#include "unfinished/rosic_ModalFilterBank.h"
#include "unfinished/rosic_ModalSynth.h"

// third party:
#include "_third_party/SampleTailExtender/SampleTailExtender.h"

// restore warning level in msvc:
#if defined _MSC_VER
#pragma warning(pop)
#endif

#endif // #ifndef ROSIC_H_INCLUDED
