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
// RAPT instantiations:
#include <rapt/rapt.cpp>

// We request some explicit instantiations here - later, when we add modules to the jura framework
// which use these classes, they may be deleted. At the moment, they are needed for Elan's
// Chaosfly but are nowhere instantiatied within jura. It's not a very elegant solution, but it's
// supposed to be temporary anyway:

//template double RAPT::rsLinToLin(double x, double inMin, double inMax, double outMin, double outMax);

// Ouura FFT instantiations (trying to fix linker errors on mac):
template void RAPT::bitrv2conj(int, int*, double*);
template void RAPT::bitrv2(int, int*, double*);
template void RAPT::makect(int, int*, double*);
template void RAPT::makewt(int, int*, double*);
template void RAPT::cftbsub(int, double*, double*);
template void RAPT::cftfsub(int, double*, double*);
template void RAPT::rftbsub(int, double*, int, double*);
template void RAPT::rftfsub(int, double*, int, double*);
// ...hmm - but it doesn't seem to help - linker errors persist
// ...let's try these:
template void RAPT::cdft(int, int, double *, int *, double *);
template void RAPT::rdft(int, int, double *, int *, double *);
template void RAPT::ddct(int, int, double *, int *, double *);
template void RAPT::ddst(int, int, double *, int *, double *);
template void RAPT::dfct(int, double *, double *, int *, double *);
template void RAPT::dfst(int, double *, double *, int *, double *);
// ...nope - doesn't make difference
// hmm - the error says: Undefined symbols for architecture i386:
// ...ok - linker error gone after changing the fft4g.cpp and using the RAPT version in rosic, too
// ...try what happens now, if we delete these instantiations again (and maybe revert fft4g to its
// old state - maybe it was using the rapt version in rosic that made the difference...)



template class RAPT::rsMatrixView<double>;
template class RAPT::rsMatrixNew<double>;

template class RAPT::rsNodeBasedFunction<double>;
template class RAPT::rsParametricBellFunction<double>;
template class RAPT::rsPositiveBellFunctions<double>;
template class RAPT::rsPositiveSigmoids<double>;
template class RAPT::rsNormalizedSigmoids<double>;
template class RAPT::rsParametricSigmoid<double>;
template class RAPT::rsSinCosTable<double>;
template class RAPT::rsScaledAndShiftedSigmoid<double>;
template class RAPT::rsEllipse<double>;
template class RAPT::rsRotationXY<double>;
template class RAPT::rsRotationXYZ<double>;
template class RAPT::rsCoordinateMapper<double>;
template class RAPT::rsCoordinateMapper2D<double>;

template class RAPT::rsFourierTransformerRadix2<double>;


template class RAPT::rsStateVariableFilter<double, double>;

template class RAPT::rsImage<float>;
template class RAPT::rsAlphaMask<float>;
template class RAPT::rsImagePainter<float, float, float>;

// for PhaseScope:
template class RAPT::rsAlphaMask<double>;
template class RAPT::rsScopeScreenScanner<float>;  // do we need this?
template class RAPT::rsScopeScreenScanner<double>;
template class RAPT::rsPhaseScopeBuffer<double, float, double>;
template class RAPT::rsPhaseScopeBuffer2<double, float, double>;

// needed for the release build of Chaosfly on Linux - without them, apparently the compiler
// generates the classes only partially - some member functions are missing probably because they
// not called from anywhere inside jura or rosic:
template double RAPT::rsAbs(double x);
template class RAPT::rsBreakpointModulator<double>;

//template class RAPT::rsBiquadCascade<double, double>;
template class RAPT::rsBiquadCascade<rsFloat64x2, double>;

template class RAPT::rsPrototypeDesigner<double>;
template class RAPT::rsPoleZeroMapper<double>;
template class RAPT::rsFilterCoefficientConverter<double>;
template class RAPT::rsInfiniteImpulseResponseDesigner<double>;
template class RAPT::rsEngineersFilter<double, double>;
template class RAPT::rsEngineersFilter<rsFloat64x2, double>;
template struct RAPT::rsFilterSpecificationZPK<double>;
template struct RAPT::rsFilterSpecificationBA<double>;

template class RAPT::rsSmoothingFilter<double, double>;
template class RAPT::rsLadderFilter<double, double>;
template class RAPT::rsLadderFilter<rsFloat64x2, double>;

#ifdef _MSC_VER
template class RAPT::rsLadderFilter<rsFloat64x2, rsFloat64x2>;
// does not compile on mac because of std::complex<rsFloat64x2>
#endif

//template class RAPT::rsLadderFilter<std::complex<double>, double>; // needed for TestPluginJUCE

template class RAPT::rsPhasorFilter<double, double>;
template class RAPT::rsPhasorStateMapper<double>;
// todo: get rid of directly using rapt classes in jura and/or products - create instantiations for
// double in rosic and use these instantiations only

template class RAPT::rsBouncillator<double>;
template class RAPT::rsRayBouncer<double>;
template class RAPT::rsRayBouncerDriver<double>;
template class RAPT::rsLissajousOscillator3D<double>;
template class RAPT::rsEllipseOscillator<double>;
template class RAPT::rsTriSawOscillator<double>;

template class RAPT::rsMultiBandSplitter<double, double>;
template class RAPT::rsMultiBandSplitter<rsFloat64x2, double>;

template class RAPT::rsHalfWaveSaturator<double, double>;
template class RAPT::rsSaturator<double, double>;
template class RAPT::rsSlewRateLimiterLinear<double, double>;
//template class RAPT::rsBreakpointModulator<double>;

template class RAPT::rsOnePoleFilter<double, double>;

template class RAPT::rsModalFilter<double, double>;
template class RAPT::rsNonlinearModalFilter<double, double>;
template class RAPT::rsModalFilterBank<double, double>;
template class RAPT::rsModalFilterWithAttack2<double, double>;

//template class RAPT::rsStateVariableFilter<double, double>;
template class RAPT::rsPhonoFilter<double, double>;
template class RAPT::rsMovingAverage<double, double>;


//template class RAPT::rsPrototypeDesigner<double>;
//template class RAPT::rsInfiniteImpulseResponseDesigner<double>;
//template class RAPT::rsEngineersFilter<double, double>;

template class RAPT::rsLadderFilter2<double, double>;
template class RAPT::rsLadderFilterZDF<double, double>;
template class RAPT::rsLadderResoShaped<double, double>;
template class RAPT::rsLadderResoShaped2<double, double>;
template class RAPT::rsLadderFilterFeedbackSaturated<double, double>;
template class RAPT::rsResoReplacer<double, double>;
template class RAPT::rsResoReplacerPhaseBumped<double, double>;
template class RAPT::rsFakeResonanceFilter<double, double>;
template class RAPT::rsLadderMystran<double, double>;

template class RAPT::rsDelayLine<double, double>;
template class RAPT::rsFractionalDelayLine<double, double>;


template class RAPT::rsInstantaneousFundamentalEstimator<double>; // rename

template class RAPT::rsZeroCrossingPitchDetector<double>;
template class RAPT::rsAutoCorrelationPitchDetector<double>;

template class RAPT::rsPhaseVocoder<double>;

template class RAPT::rsDoublePendulum<double, double>;

template class RAPT::rsResampler<double, double>;
template class RAPT::rsTimeWarper<double, double>;
//template class RAPT::rsInstantaneousFundamentalEstimator<double>;
template class RAPT::rsCycleMarkFinder<double>;
template class RAPT::rsVariableSpeedPlayer<double, double>;
template class RAPT::rsPhaseLockedCrossfader<double, double>;


// hmm...it seems, we need all these explicit instantiations anyway - maybe clean up the build
// system...rename the files to raptJuceModule.h/cpp




//=================================================================================================
// rosic includes

// basics (but we needed to intersperse some stuff from other folders)
#include "basics/GlobalFunctions.cpp"
#include "basics/rosic_ChannelMatrix2x2.cpp"
#include "basics/rosic_Constants.cpp"                        // empty
#include "basics/rosic_FunctionTemplates.cpp"                // empty
#include "basics/rosic_HelperFunctions.cpp"
#include "basics/rosic_Interpolator.cpp"
#include "basics/rosic_NumberManipulations.cpp"              // empty
#include "infrastructure/rosic_MutexLock.cpp"                // used by SampleBuffer
#include "basics/rosic_SampleBuffer.cpp"
#include "basics/rosic_SamplePlaybackParameters.cpp"
#include "math/rosic_ElementaryFunctionsReal.cpp"            // used by SpecialFunctionsReal?
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
#include "datastructures/rosic_Array.cpp"                    // empty
#include "datastructures/rosic_ExtensionsForSTL.cpp"         // empty
#include "datastructures/rosic_KeyValueMap.cpp"              // empty
#include "datastructures/rosic_String.cpp"

// math (some of the cpp files in this folder are already included in the basics section):
#include "math/rosic_IntegerFunctions.cpp"
#include "math/rosic_LinearAlgebra.cpp"
#include "math/rosic_Vector.cpp"
#include "math/rosic_Matrix.cpp"
#include "math/rosic_MatrixVectorFunctions.cpp"
#include "math/rosic_Interpolation.cpp"
#include "math/rosic_PolynomialAlgorithms.cpp"
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
#include "filters/rosic_FilterAnalyzer.cpp"
#include "filters/rosic_BiquadCascade.cpp"
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
#include "filters/rosic_PrototypeDesigner.cpp"
#include "filters/rosic_PoleZeroMapper.cpp"
#include "filters/rosic_FilterCoefficientConverter.cpp"
#include "filters/rosic_InfiniteImpulseResponseDesigner.cpp"
#include "filters/rosic_EngineersFilter.cpp"
#include "filters/rosic_EngineersFilterOld.cpp"
#include "filters/rosic_LinkwitzRileyCrossOver.cpp"
#include "filters/rosic_CrossOver4Way.cpp"
#include "filters/rosic_DirectFormFilter.cpp"
#include "filters/rosic_TwoPoleFilter.cpp"
#include "filters/rosic_DualTwoPoleFilter.cpp"
#include "filters/rosic_EllipticQuarterBandFilter.cpp"
#include "filters/rosic_EllipticSubBandFilter.cpp"
#include "filters/rosic_EllipticSubBandFilterDirectForm.cpp"
#include "filters/rosic_Equalizer.cpp"
#include "filters/rosic_EqualizerStereo.cpp"
#include "filters/rosic_FiniteImpulseResponseDesigner.cpp"
#include "filters/rosic_FiniteImpulseResponseFilter.cpp"
#include "filters/rosic_FourPoleFilter.cpp"
#include "filters/rosic_OnePoleFilter.cpp"
#include "filters/rosic_OnePoleFilterStereo.cpp"
#include "filters/rosic_LadderFilter.cpp"
#include "filters/rosic_LadderFilterOld.cpp"
#include "filters/rosic_LeakyIntegrator.cpp"
#include "filters/rosic_LowpassHighpass.cpp"
#include "filters/rosic_LowpassHighpassStereo.cpp"
#include "filters/rosic_LpfHpfApf.cpp"
#include "filters/rosic_MultiModeFilter.cpp"
#include "filters/rosic_NyquistBlocker.cpp"
#include "filters/rosic_QuadratureNetwork.cpp"
#include "filters/rosic_SlopeFilter.cpp"
#include "filters/rosic_TeeBeeFilter.cpp"
#include "filters/rosic_ToneControl.cpp"
#include "filters/rosic_TwoPoleBandpass.cpp"
#include "filters/rosic_VowelFilterStereo.cpp"
#include "filters/rosic_WarpedBiquadMonoDF1.cpp"            // empty
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
#include "generators/rosic_ModalSynthesizer.cpp"            // replace with code from RSLib
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
#include "others/rosic_ProcessorCycleCounter.cpp"
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

#include "unfinished/rosic_Polyphony.cpp"
#include "unfinished/rosic_EllipseOscillator.cpp"
#include "unfinished/rosic_QuadSource.cpp"
#include "unfinished/rosic_DualFilter.cpp"
#include "unfinished/rosic_NewSynth.cpp"

// restore warning level in msvc:
#if defined _MSC_VER
#pragma warning(pop)
#endif
