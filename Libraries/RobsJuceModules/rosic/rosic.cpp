#ifdef ROSIC_H_INCLUDED
/* When you add this cpp file to your project, you mustn't include it in a file where you've
already included any other headers - just put it inside a file on its own, possibly with your config
flags preceding it, but don't include anything else. That also includes avoiding any automatic prefix
header files that the compiler may be using. */
#error "Incorrect use of JUCE cpp file"
#endif

#include "rosic.h"

//=================================================================================================

/** The cpp files are included in the order in which they depend on each other. The only folder
from which we have not included the .cpp files is the legacy folder.

ToDo: reorder them in the library accordingly, such that files in one library folder depend only on
other files in folders that are considered "above" in the hierarchy.
In the future, rosic should be made dependent on the RAPT library and whereever it makes sense, the
rosic-class should be turned into an template instatiation of  RAPT class template.
jura_processors should depend on the rosic module and grab its DSP code from there
rename modules:

jura_framework: rs_framework
jura_processors: rs_audio_processors
...namespace name should be rs (but do all of this only after dragging in all old plugin code and
integrating it into the Chainer.

rosic: rs_dsp (this should never depend on any juce class/module)
 -..or maybe not - maybe it's good to be in its own namespace because we have a Module class here
  and also an AudioModule class in jura...maybe keep it named rosic for the time being
 -or maybe we use the rs prefix everywhere
 -maybe call this modules rs_dsp_realtime and make another rs_dsp_offline (which depends on the
  realtime module)
*/

// in msvc, set warning level to 3:
#if defined _MSC_VER
#pragma warning(push, 3)
#endif

//=================================================================================================
// rosic includes

// basics (but we needed to intersperse some stuff from other folders)
#include "basics/rosic_TemplateInstantiations.cpp"
#include "basics/rosic_ChannelMatrix2x2.cpp"
#include "basics/rosic_HelperFunctions.cpp"
#include "basics/rosic_Interpolator.cpp"
#include "infrastructure/rosic_MutexLock.cpp"                // used by SampleBuffer
#include "basics/rosic_SampleBuffer.cpp"
#include "basics/rosic_SamplePlaybackParameters.cpp"
#include "math/rosic_RealFunctionEvaluationAlgorithms.cpp"   // used by SpecialFunctionsReal
#include "math/rosic_SpecialFunctionsReal.cpp"               // used by ComplexFunctions?
#include "math/rosic_Complex.cpp"                            // used by ComplexFunctionsAlgorithms
#include "math/rosic_ComplexFunctionsAlgorithms.cpp"
#include "math/rosic_ComplexFunctions.cpp"                   // used by ExpressionEvaluator
#include "scripting/rosic_ExpressionEvaluator.cpp"           // used by TabulatedFunction
#include "basics/rosic_TabulatedFunction.cpp"                // needs ExpressionEvaluator
#include "basics/rosic_WarpedAllpassInterpolator.cpp"
#include "basics/rosic_WindowDesigner.cpp"

// datastructures
#include "datastructures/rosic_String.cpp"

// math (some of the cpp files in this folder are already included in the basics section):
//#include "math/rosic_IntegerFunctions.cpp"
#include "math/rosic_LinearAlgebra.cpp"
#include "math/rosic_Vector.cpp"
#include "math/rosic_Matrix.cpp"
#include "math/rosic_MatrixVectorFunctions.cpp"
#include "math/rosic_PrimeNumbers.cpp"
#include "math/rosic_Transformations.cpp"

// transforms
#include "transforms/rosic_FourierTransformerRadix2.cpp"
#include "transforms/rosic_FourierTransformerBluestein.cpp"
#include "transforms/rosic_WaveletTransforms.cpp"
#include "transforms/rosic_SpectralManipulator.cpp"  // maybe move to other folder

// filters
#include "filters/rosic_AllpassChain.cpp"
#include "filters/rosic_BiquadBase.cpp"
#include "filters/rosic_BiquadDesigner.cpp"
#include "filters/rosic_BiquadMonoDF1.cpp"
#include "filters/rosic_BiquadStereoDF1.cpp"
#include "filters/rosic_BiquadStereoDF2.cpp"
#include "filters/rosic_CombFilter.cpp"
#include "filters/rosic_DampingFilter.cpp"
#include "filters/rosic_CombResonator.cpp"
#include "filters/rosic_ConvolverBruteForce.cpp"
#include "filters/rosic_ConvolverFFT.cpp"
#include "filters/rosic_ConvolverPartitioned.cpp"
#include "filters/rosic_CookbookFilter.cpp"
#include "filters/rosic_TwoPoleFilter.cpp"
#include "filters/rosic_DualTwoPoleFilter.cpp"
#include "filters/rosic_EllipticQuarterBandFilter.cpp"
#include "filters/rosic_Equalizer.cpp"
#include "filters/rosic_EqualizerStereo.cpp"
#include "filters/rosic_FiniteImpulseResponseDesigner.cpp"
#include "filters/rosic_FiniteImpulseResponseFilter.cpp"
#include "filters/rosic_FourPoleFilter.cpp"
#include "filters/rosic_LadderFilter.cpp"
#include "filters/rosic_LadderFilterOld.cpp"
#include "filters/rosic_LeakyIntegrator.cpp"
#include "filters/rosic_LowpassHighpass.cpp"
#include "filters/rosic_LowpassHighpassStereo.cpp"
#include "filters/rosic_LpfHpfApf.cpp"
#include "filters/rosic_MultiModeFilter.cpp"
#include "filters/rosic_NyquistBlocker.cpp"
#include "filters/rosic_SlopeFilter.cpp"
#include "filters/rosic_TeeBeeFilter.cpp"
#include "filters/rosic_ToneControl.cpp"
#include "filters/rosic_TwoPoleBandpass.cpp"
#include "filters/rosic_VowelFilterStereo.cpp"
#include "filters/rosic_WhiteToPinkFilter.cpp"

// rendering
//#include "rendering/rosic_AdditveWaveformRenderer.cpp"    // obsolete?
#include "rendering/rosic_AlgorithmicWaveformRenderer.cpp"
#include "rendering/rosic_MultiSegmentWaveformRenderer.cpp"
#include "rendering/rosic_NonRealtimeProcesses.cpp"         // should be renamed
#include "rendering/rosic_StandardWaveformRenderer.cpp"
#include "rendering/rosic_WaveformBuffer.cpp"
#include "rendering/rosic_WaveformRenderer.cpp"
#include "rendering/rosic_TurtleGraphics.cpp"
#include "rendering/rosic_LindenmayerSystem.cpp"

// generators
#include "generators/rosic_MipMappedWaveTable.cpp"
#include "generators/rosic_MipMappedWaveTableOld.cpp"       // can this be removed?
#include "generators/rosic_MipMappedWaveTableStereo.cpp"
#include "generators/rosic_BlendOscillator.cpp"
#include "generators/rosic_LorentzSystem.cpp"
#include "generators/rosic_NoiseGenerator.cpp"
#include "generators/rosic_NoiseGeneratorOld.cpp"           // can this be removed?
#include "generators/rosic_Oscillator.cpp"
#include "generators/rosic_OscillatorBank.cpp"
#include "generators/rosic_OscillatorStereo.cpp"
#include "generators/rosic_FourOscSection.cpp"
#include "generators/rosic_SampleOscillator.cpp"
#include "generators/rosic_SamplePlayer.cpp"
#include "generators/rosic_SineOscillator.cpp"
#include "generators/rosic_SineOscillatorStereo.cpp"
#include "generators/rosic_SuperOscillator.cpp"
#include "generators/rosic_TestGenerator.cpp"
#include "generators/rosic_WaveTable.cpp"
#include "generators/rosic_TurtleSource.cpp"
#include "generators/rosic_Snowflake.cpp"

// modulators
//#include "modulators/MagicCarpetModulator.cpp" // needs MagicCarpetDefinitions.h - where is this? legacy?
#include "modulators/rosic_AmpEnvRc.cpp"
#include "modulators/rosic_AnalogEnvelope.cpp"
#include "modulators/rosic_AnalogEnvelopeScaled.cpp"
#include "modulators/rosic_VariousModulators.cpp"
#include "modulators/rosic_BreakpointModulator.cpp"
#include "modulators/rosic_DecayEnvelope.cpp"
#include "modulators/rosic_EnvelopeGenerator.cpp"
#include "modulators/rosic_EnvelopeGenerator2.cpp"  // do we need this?
#include "modulators/rosic_ExponentialRamp.cpp"
#include "others/rosic_SlewRateLimiter.cpp"         // move to basics
#include "others/rosic_SlewRateLimiterLinear.cpp"
//#include "others/rosic_SlewRateLimiterOld.cpp"    // not used anymore
#include "modulators/rosic_LowFrequencyOscillator.cpp" // needs SlewRateLimiter, WaveTable
#include "modulators/rosic_PitchEnvRc.cpp"
#include "modulators/rosic_SampleModulator.cpp"
#include "modulators/rosic_Modulator.cpp"              // needs SampleModulator

// others
//#include "others/rosic_DemoVersionNoiseEmitter.cpp"    // may not be needed
#include "others/rosic_RandomNumberGenerator01.cpp"
#include "others/rosic_RandomNumberGenerator02.cpp"
//#include "others/rosic_KeyGenerator.cpp"             // remove...or keep only the "Validator" part
#include "others/rosic_OverlapAddProcessor.cpp"
#include "others/rosic_PiecewiseFunction.cpp"
//#include "others/rosic_Plotter.cpp"                  // obsolete - use GNUPlotCPP now
#include "others/rosic_ProcessorCycleCounter.cpp"    // obsolete
#include "others/rosic_RoutingMatrix.cpp"
#include "others/rosic_SpectralProcessor.cpp"
#include "others/rosic_SpectralEnvelopeProcessor.cpp"
#include "others/rosic_SpectralFilter.cpp"
#include "others/rosic_TuningTable.cpp"
#include "others/rosic_VectorMixer.cpp"

// analysis
#include "analysis/rosic_CyclicAutoCorrelator.cpp"
#include "analysis/rosic_EnvelopeFollower.cpp"    // needs SlewRateLimiter
#include "analysis/rosic_LinearPredictor.cpp"
#include "analysis/rosic_FormantRemover.cpp"
#include "analysis/rosic_FormantPreserver.cpp"
#include "analysis/rosic_InstantaneousEnvelopeDetector.cpp"
#include "analysis/rosic_LevelDetector.cpp"
#include "analysis/rosic_OscilloscopeBufferOld.cpp"
#include "analysis/rosic_PitchDetector.cpp"           // may be buggy - compare to RSLib version
#include "analysis/rosic_SignalMeasures.cpp"
#include "analysis/rosic_SpectrumAnalyzer.cpp"
#include "analysis/rosic_TrackMeter.cpp"
#include "analysis/rosic_WaveformDisplayBuffer.cpp"
#include "analysis/rosic_ScopeScreenScanner.cpp"
#include "analysis/rosic_OnsetDetector.cpp"
#include "analysis/rosic_BeatDetector.cpp"

// delaylines
#include "delaylines/rosic_BasicIntegerDelayLine.cpp"
#include "delaylines/rosic_IntegerDelayLine.cpp"
#include "delaylines/rosic_AllpassDiffusor.cpp"
#include "delaylines/rosic_DelayLineStereo.cpp"
#include "delaylines/rosic_EchoLabDelayLine.cpp"
#include "delaylines/rosic_FractionalDelayLine.cpp"
#include "delaylines/rosic_FractionalDelayLineStereo.cpp"
#include "delaylines/rosic_ModulatedDelayLine.cpp"
#include "delaylines/rosic_PingPongEcho.cpp"

// dynamics
#include "dynamics/rosic_BesselFilterForGainSignal.cpp"
#include "dynamics/rosic_DynamicsProcessorBase.cpp"
#include "dynamics/rosic_Compressor.cpp"
#include "dynamics/rosic_Expander.cpp"
#include "dynamics/rosic_Limiter.cpp"
#include "dynamics/rosic_NoiseGate.cpp"
#include "dynamics/rosic_SoftKneeCompressor.cpp"
#include "dynamics/rosic_SoftKneeExpander.cpp"
#include "dynamics/rosic_MultiBandCompressor.cpp"
// where's the Leveller?
// make a dynamics processor with freely adjustable curve

// effects
#include "effects/rosic_FeedbackDelayNetwork.cpp"
#include "effects/rosic_FeedbackDelayNetwork8.cpp"
#include "effects/rosic_FeedbackDelayNetwork16.cpp"
#include "effects/rosic_AlgoVerb.cpp"
#include "effects/rosic_AudioToMidi.cpp"
#include "effects/rosic_BitCrusher.cpp"
#include "effects/rosic_ModulationEffect.cpp"
#include "effects/rosic_Vibrato.cpp"
#include "effects/rosic_Chorus.cpp"
#include "effects/rosic_CombBank.cpp"
#include "effects/rosic_CombResonatorStereo.cpp"
#include "effects/rosic_CombStereoizer.cpp"
#include "effects/rosic_CompShaper.cpp"
#include "effects/rosic_Phaser.cpp"
#include "effects/rosic_DelayPhaser.cpp"
#include "effects/rosic_Distortion.cpp"
#include "effects/rosic_EchoLab.cpp"
#include "effects/rosic_Flanger.cpp"
#include "effects/rosic_FormantShifter.cpp"
#include "effects/rosic_FrequencyShifter.cpp"
#include "effects/rosic_FuncShaper.cpp"
#include "effects/rosic_Harmonics.cpp"
#include "effects/rosic_ModulatedAllpass.cpp"
#include "effects/rosic_Moduluxury.cpp"
#include "effects/rosic_Noisifier.cpp"
#include "effects/rosic_PhaseStereoizer.cpp"
#include "effects/rosic_PitchShifter.cpp"
#include "effects/rosic_PitchShifterGrainAdaptive.cpp"
#include "effects/rosic_Reverb.cpp"
#include "effects/rosic_RingModulator.cpp"
#include "effects/rosic_SingleSidebandModulator.cpp"
#include "effects/rosic_StereoDelay.cpp"
#include "effects/rosic_StereoPan.cpp"
#include "effects/rosic_StereoWidth.cpp"
#include "effects/rosic_Tremolo.cpp"
#include "effects/rosic_WahWah.cpp"
#include "effects/rosic_WaveShaper.cpp"

// infrastructure
#include "infrastructure/rosic_ModulationRouting.cpp"
#include "infrastructure/rosic_Module.cpp"
#include "infrastructure/rosic_EffectModules.h"
#include "infrastructure/rosic_GeneratorModules.cpp"
#include "infrastructure/rosic_ModulatorModules.cpp"
#include "infrastructure/rosic_File.cpp"
#include "infrastructure/rosic_FileInputOutput.cpp"
#include "infrastructure/rosic_MemoryUser.cpp"        // where is this used?
#include "infrastructure/rosic_MidiNoteEvent.cpp"
#include "infrastructure/rosic_PolyphonicInstrument.cpp"
#include "infrastructure/rosic_PolyphonicInstrumentVoice.cpp"
// copy MonoSynth baseclass from ChaosGenerator here ...but maybe PolyphonicInstrumentVoice has
// already the required functionality and more?

// numerical (maybe move to math)
#include "numerical/rosic_FunctionObjects.cpp"
#include "numerical/rosic_GradientBasedMinimizer.cpp"

// neural
//#include "neural/rosic_MultiLayerPerceptron.cpp"
//#include "neural/rosic_MultiLayerPerceptronErrorFunction.cpp"
//#include "neural/rosic_MultiLayerPerceptronTrainer.cpp"

// sequencing
#include "sequencing/rosic_AcidPattern.cpp"
#include "sequencing/rosic_AcidSequencer.cpp"

// offline
// ...include code from RSLib for offline processes (Resampler, PitchFlattener, etc.) here...
// ...also the BeatDetector stuff

// scripting
//#include "scripting/rosic_AngelScriptInterpreter.cpp"      // not yet used anywhere
#include "scripting/rosic_DspScriptInterpreter.cpp"
#include "scripting/rosic_DspWorkbench.cpp"

// instruments
#include "instruments/rosic_AciDevil.cpp"
#include "instruments/rosic_KeyShotVoice.cpp"
#include "instruments/rosic_KeyShot.cpp"
#include "instruments/rosic_MagicCarpetVoice.cpp"
#include "instruments/rosic_MagicCarpet.cpp"
#include "instruments/rosic_Open303.cpp"
#include "instruments/rosic_QuadrigaVoice.cpp"
#include "instruments/rosic_Quadriga.cpp"
#include "instruments/rosic_SimpleSamplerOscSection.cpp"
#include "instruments/rosic_SimpleSamplerVoice.cpp"
#include "instruments/rosic_SimpleSampler.cpp"
#include "instruments/rosic_StraightlinerVoice.cpp"
#include "instruments/rosic_Straightliner.cpp"
#include "instruments/rosic_WorkhorseOscSection.cpp"
#include "instruments/rosic_WorkhorseVoice.cpp"
#include "instruments/rosic_Workhorse.cpp"

// these do not really fit into the directory order - they are kind-of higher-level classes, maybe
// they should be in the instruments section...or something:
#include "effects/rosic_Quadrifex.cpp"  // needs the effects wrapped into "Module" subclasses
#include "generators/rosic_Quadrigen.cpp"             // needs Module, RoutingMatrix
#include "generators/rosic_VectorSamplePlayer.cpp"    // needs LowFrequencyOscillator (in modulators)
                                                      // and VectorMixer

#include "unfinished/rosic_MusiciansFilter.cpp"
#include "unfinished/rosic_Polyphony.cpp"
#include "unfinished/rosic_EllipseOscillator.cpp"
#include "unfinished/rosic_QuadSource.cpp"
#include "unfinished/rosic_DualFilter.cpp"
#include "unfinished/rosic_NewSynth.cpp"
#include "unfinished/rosic_MiscUnfinished.cpp"
#include "unfinished/rosic_ModalFilters.cpp"
#include "unfinished/rosic_ModalFilterBank.cpp"
#include "unfinished/rosic_ModalSynth.cpp"

#include "unfinished/rosic_AudioStream.cpp"

#include "unfinished/sampler/rosic_SamplerTools.cpp"
#include "unfinished/sampler/rosic_SfzCodeBook.cpp"
#include "unfinished/sampler/rosic_SamplerData.cpp"
#include "unfinished/sampler/rosic_SamplerEffectCores.cpp"
#include "unfinished/sampler/rosic_SamplerProcessors.cpp"
#include "unfinished/sampler/rosic_SamplerPlayers.cpp"
#include "unfinished/sampler/rosic_SamplerEngine.cpp"


// third party:
#include "_third_party/SampleTailExtender/FFT.cpp"
//#include "_third_party/SampleTailExtender/libs/kiss_fft130/kiss_fft.c" // clashes with already included kiss_fft_ v1_2_6 - can we get rid of that?
#include "_third_party/SampleTailExtender/HarmonicAnalyser.cpp"
#include "_third_party/SampleTailExtender/SampleTailExtender.cpp"


// restore warning level in msvc:
#if defined _MSC_VER
#pragma warning(pop)
#endif
