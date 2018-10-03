#ifndef rosic_TemplateInstantiations_h
#define rosic_TemplateInstantiations_h

// This file contains typedefs for explicit template instantiations for templates from the RAPT 
// library. Sometimes we also create a subclass of a template call from RAPT in order to provide
// additional functionality such as converting between two doubles and rsFloat64x2 etc. to make it
// more convenient to use the classes

namespace rosic
{

typedef RAPT::rsSinCosTable<double> rsSinCosTableD;

typedef RAPT::rsSmoothingFilter<double, double> rsSmoothingFilterDD;
typedef RAPT::rsLadderFilter<double, double> rsLadderDD;
typedef RAPT::rsLadderFilter<rsFloat64x2, double> rsLadderD2D;

typedef RAPT::rsRayBouncer<double> rsRayBouncerD;
typedef RAPT::rsRayBouncerDriver<double> rsRayBouncerDriverD;



typedef RAPT::rsPositiveSigmoids<double> rsPositiveSigmoidsD;
typedef RAPT::rsNormalizedSigmoids<double> rsNormalizedSigmoidsD;
typedef RAPT::rsParametricSigmoid<double> rsParametricSigmoidD;

typedef RAPT::rsFourierTransformerRadix2<double> rsFourierTransformerRadix2D;


typedef RAPT::rsSaturator<double, double> rsSaturatorDD;
typedef RAPT::rsSlewRateLimiterLinear<double, double> rsSlewRateLimiterLinearDD;
typedef RAPT::rsBreakpointModulator<double> rsBreakpointModulatorD;

typedef RAPT::rsBasicDelayLine<double> rsBasicDelayLineD;
typedef RAPT::rsDelayLine<double, double> rsDelayLineDD;
typedef RAPT::rsFractionalDelayLine<double, double> rsFractionalDelayLineDD;

typedef RAPT::rsOnePoleFilter<double, double> rsOnePoleFilterDD;


typedef RAPT::rsNonlinearModalFilter<double, double> rsNonlinearModalFilterDD;
typedef RAPT::rsModalFilter<double, double> rsModalFilterDD;
typedef RAPT::rsModalFilterWithAttack<double, double> rsModalFilterWithAttackDD;
typedef RAPT::rsModalFilterWithAttack2<double, double> rsModalFilterWithAttack2DD;
typedef RAPT::rsModalFilterBank<double, double> rsModalFilterBankDD;


typedef RAPT::rsStateVariableFilter<double, double> rsStateVariableFilterDD;
typedef RAPT::rsPhonoFilter<double, double> rsPhonoFilterDD;
typedef RAPT::rsMovingAverage<double, double> rsMovingAverageDD;

typedef RAPT::rsPrototypeDesigner<double> rsPrototypeDesignerD;
typedef RAPT::rsPoleZeroMapper<double> rsPoleZeroMapperD;
typedef RAPT::rsFilterCoefficientConverter<double> rsFilterCoefficientConverterD;
typedef RAPT::rsInfiniteImpulseResponseDesigner<double> rsInfiniteImpulseResponseDesignerD;

//typedef RAPT::rsEngineersFilter<double, double> rsEngineersFilterDD;
//typedef RAPT::rsEngineersFilter<rsFloat64x2, double> rsEngineersFilterD2D;

//typedef RAPT::rsEllipticSubBandFilterDirectForm<double, double> rsEllipticSubBandFilterDirectFormMono;
typedef RAPT::rsEllipticSubBandFilterDirectForm<double, double> rsSubBandFilterMono;


typedef RAPT::rsLadderFilter2<double, double> rsLadderFilter2DD;  // the version from RSLib
typedef RAPT::rsLadderFilterZDF<double, double> rsLadderFilterZDFDD;
typedef RAPT::rsLadderResoShaped<double, double> rsLadderResoShapedDD;
typedef RAPT::rsLadderResoShaped2<double, double> rsLadderResoShaped2DD;
typedef RAPT::rsLadderFilterFeedbackSaturated<double, double> rsLadderFeedbackSaturatedDD;
typedef RAPT::rsResoReplacer<double, double> rsResoReplacerDD;
typedef RAPT::rsResoReplacerPhaseBumped<double, double> rsResoReplacerPhaseBumpedDD;
typedef RAPT::rsFakeResonanceFilter<double, double> rsFakeResonanceFilterDD;
typedef RAPT::rsLadderMystran<double, double> rsLadderMystranDD;


typedef RAPT::rsInstantaneousFundamentalEstimator<double> rsInstantaneousFundamentalEstimatorD;
typedef RAPT::rsZeroCrossingPitchDetector<double> rsZeroCrossingPitchDetectorD;
typedef RAPT::rsAutoCorrelationPitchDetector<double> rsAutoCorrelationPitchDetectorD;

typedef RAPT::rsPhaseVocoder<double> rsPhaseVocoderD;

typedef RAPT::rsDoublePendulum<double, double> rsDoublePendulumDD;

typedef RAPT::rsResampler<double, double> rsResamplerDD;
typedef RAPT::rsTimeWarper<double, double> rsTimeWarperDD;
typedef RAPT::rsInstantaneousFundamentalEstimator<double> rsInstantaneousFundamentalEstimatorD;
typedef RAPT::rsCycleMarkFinder<double> rsCycleMarkFinderD;
typedef RAPT::rsVariableSpeedPlayer<double, double> rsVariableSpeedPlayerDD;
typedef RAPT::rsPhaseLockedCrossfader<double, double>  rsPhaseLockedCrossfaderDD;

}


#endif