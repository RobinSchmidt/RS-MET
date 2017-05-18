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

// ... good, until here

// modulators:

// rendering:

// delaylines needs: filters, generators


// some more complex generators that need includes from modulators and rendering:
#include "generators/rosic_Quadrigen.h"             // needs BreakpointModulator, infrastructure/*Modules
#include "generators/rosic_VectorSamplePlayer.h"    // needs LowFrequencyOscillator
#include "generators/rosic_WaveTable.h"             // needs WaveformRenderer

//-------------------------------------------------------------------------------------------------
// old, de-centralized header includes:

#include "analysis/rosic_Analysis.h"
//#include "basics/rosic_Basics.h"
//#include "datastructures/rosic_DataStructures.h"
#include "delaylines/rosic_DelayLines.h"
#include "dynamics/rosic_Dynamics.h"
#include "effects/rosic_Effects.h"
//#include "filters/rosic_Filters.h"
//#include "generators/rosic_Generators.h"
#include "infrastructure/rosic_Infrastructure.h"
#include "instruments/rosic_Instruments.h"
//#include "math/rosic_Math.h"
#include "modulators/rosic_Modulators.h"
#include "neural/rosic_Neural.h"
#include "numerical/rosic_Numerical.h"
//#include "others/rosic_Others.h"
#include "rendering/rosic_Rendering.h"
#include "scripting/rosic_Scripting.h"
//#include "plugins/rosic_PlugIns.h"
//#include "transforms/rosic_Transforms.h"

#endif // #ifndef ROSIC_H_INCLUDED







