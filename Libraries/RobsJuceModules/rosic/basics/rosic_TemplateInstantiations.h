#ifndef rosic_TemplateInstantiations_h
#define rosic_TemplateInstantiations_h

/** This file contains typedefs for explicit template instantiations for templates from the RAPT
library. Sometimes we also create a subclass of a class template from RAPT in order to provide
additional functionality to make it more convenient to use the classes. These convenience functions
are for tasks like converting between two doubles and rsFloat64x2 in case of SIMD-type
instantiations and/or providing suitable callback target functions in the case where the underlying
RAPT class doesn't conform to the interface required by our callback system in jura etc. This code
is all really ugly administrative clutter - don't look at it, if you want to avoid eye cancer! */

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

typedef RAPT::rsBiquadCascade<double, double> rsBiquadCascadeDD;
typedef RAPT::rsDirectFormFilter<double, double> rsDirectFormFilterDD;

typedef RAPT::rsFilterAnalyzer<double> rsFilterAnalyzerD;

typedef RAPT::rsPrototypeDesigner<double> rsPrototypeDesignerD;
typedef RAPT::rsPoleZeroMapper<double> rsPoleZeroMapperD;
typedef RAPT::rsFilterCoefficientConverter<double> rsFilterCoefficientConverterD;
typedef RAPT::rsInfiniteImpulseResponseDesigner<double> rsInfiniteImpulseResponseDesignerD;
typedef RAPT::rsQuadratureNetwork<double, double> rsQuadratureNetwork;
//typedef RAPT::rsLinkwitzRileyCrossOver<double, double> rsLinkwitzRileyCrossOver;
//typedef RAPT::rsCrossOver4Way<double, double> rsCrossOver4Way;
 // crossover does not yet work in ToolChain due to stereo/mono handling

//typedef RAPT::rsEngineersFilter<double, double> rsEngineersFilterDD;
//typedef RAPT::rsEngineersFilter<rsFloat64x2, double> rsEngineersFilterD2D;

//typedef RAPT::rsEllipticSubBandFilterDirectForm<double, double> rsEllipticSubBandFilterDirectFormMono;

typedef RAPT::rsEllipticSubBandFilter<double, double> rsSubBandFilterMonoBQ; // biquad cascade
typedef RAPT::rsEllipticSubBandFilterDirectForm<double, double> rsSubBandFilterMono; // direct form - rename to DF


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

typedef RAPT::rsSpectrogramProcessor<double> rsSpectrogramD;

typedef RAPT::rsDoublePendulum<double, double> rsDoublePendulumDD;

typedef RAPT::rsResampler<double, double> rsResamplerDD;
typedef RAPT::rsTimeWarper<double, double> rsTimeWarperDD;
typedef RAPT::rsInstantaneousFundamentalEstimator<double> rsInstantaneousFundamentalEstimatorD;
typedef RAPT::rsCycleMarkFinder<double> rsCycleMarkFinderD;
typedef RAPT::rsVariableSpeedPlayer<double, double> rsVariableSpeedPlayerDD;
typedef RAPT::rsPhaseLockedCrossfader<double, double>  rsPhaseLockedCrossfaderDD;
typedef RAPT::rsEnvelopeExtractor<double> rsEnvelopeExtractorD;


//typedef RAPT::rsBlepOscArray<double, RAPT::rsBlepReadyOscBase<double>, RAPT::rsPolyBlep1<double, double>>
//  rsOscArrayPolyBlep1;

// maybe move this into its own file in the generators folder - it's grown a bit big for keeping it
// here...
class rsOscArrayPolyBlep1 : public RAPT::rsBlepOscArray<double, RAPT::rsBlepReadyOscBase<double>,
                                                        RAPT::rsPolyBlep1<double, double>>
{
public:

  typedef RAPT::rsBlepOscArray<double, RAPT::rsBlepReadyOscBase<double>,
    RAPT::rsPolyBlep1<double, double>> Base;

//  typedef public RAPT::rsBlepOscArray<double, RAPT::rsBlepReadyOscBase<double>,
//    RAPT::rsPolyBlep1<double, double>> Base;

  using Base::Base; // to inherit baseclass constructors with arguments


  rsOscArrayPolyBlep1()
  {
    ratioGenerator = new RAPT::rsRatioGenerator<double>;
    setMaxDensity(32);
    ratioGenerator->setPrimeTable(&primeTable);
  }

  ~rsOscArrayPolyBlep1()
  {
    delete ratioGenerator;
  }


  void setSampleRate(double newRate)
  {
    sampleRate = newRate;
    Base::setReferenceIncrement(frequency/sampleRate);
  }

  void setFrequency(double newFreq)
  {
    frequency = newFreq;
    Base::setReferenceIncrement(frequency/sampleRate);
  }

  void setDetunePercent(double newDetune)
  {
    Base::setDetune(0.01*newDetune);
  }

  void setMaxDensity(int newMaximum)
  {
    primeTable.resize(newMaximum+1); // +1 because of the primeSqrtDiff algo
    RAPT::rsFillPrimeTable(&primeTable[0], (RAPT::rsUint32) primeTable.size());
    Base::setMaxDensity(newMaximum);
  }

  void setFrequencyDistribution(int newDistribution)
  {
    Base::setFrequencyDistribution((RAPT::rsRatioGenerator<double>::RatioKind) newDistribution);
  }
  // needed as callback target for the valueChangeCallback in jura::Parameter - converts the
  // incoming integer (passed from the callback system) to the strongly typed enumeration type
  // required by RAPT::rsBlepOscArray



protected:

  //RAPT::rsRatioGenerator<double>* ratioGenerator = nullptr;


  std::vector<RAPT::rsUint32> primeTable;

  double sampleRate = 1, frequency = 0;

};

class rsOnePoleFilterStereo : public RAPT::rsOnePoleFilter<rsFloat64x2, double>
{
public:
  inline void getSampleFrameStereo(double *inL, double *inR, double *outL, double *outR)
  {
    rsFloat64x2 tmp = getSample(rsFloat64x2(*inL, *inR));
    *outL = tmp[0];
    *outR = tmp[1];
  }
};
// not yet in use

class rsEngineersFilterMono : public RAPT::rsEngineersFilter<double, double>
{
public:
  inline double getSample(double in) { return getSampleDirect2(in); }
  // maybe implement in RAPT::rsEngineersFilter and don't subclass here but use a typedef'd
  // explicit instantiation
};

class rsEngineersFilterStereo : public RAPT::rsEngineersFilter<rsFloat64x2, double>
{
public:
  inline void getSampleFrameStereo(double* left, double* right)
  {
    //rsFloat64x2 tmp = getSample(rsFloat64x2(*left, *right));
    rsFloat64x2 tmp = getSampleDirect2(rsFloat64x2(*left, *right));
    *left  = tmp[0];
    *right = tmp[1];
  }
};

class rsCrossOver4WayStereo : public RAPT::rsCrossOver4Way<rsFloat64x2, double>
{
public:
  inline void processSampleFrameStereo(double* inOut)
  {
    processSampleFrame(rsCastPointer(inOut));
  }
};






}


#endif
