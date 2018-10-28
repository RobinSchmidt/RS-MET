#include <rapt/rapt.cpp>

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

// Interpolation:
template void RAPT::fitCubicWithDerivative(double x1, double x2, double y1, double y2, double yd1,
  double yd2, double *a3, double *a2, double *a1, double *a0);
template void RAPT::fitCubicWithDerivativeFixedX(double y0, double y1, double yd0, double yd1,
  double *a3, double *a2, double *a1, double *a0);
template void RAPT::fitQuinticWithDerivativesFixedX(double y0, double y1, double yd0, double yd1,
  double ydd0, double ydd1, double *a5, double *a4, double *a3, double *a2, double *a1,
  double *a0);
template void RAPT::getHermiteCoeffsM(double *y0, double *y1, double *a, int M);
template void RAPT::getHermiteCoeffs1(double *y0, double *y1, double *a);
template void RAPT::getHermiteCoeffs2(double *y0, double *y1, double *a);
template void RAPT::getHermiteCoeffs3(double *y0, double *y1, double *a);
template double RAPT::getDelayedSampleAsymmetricHermiteM(double d, double *y, int M, double shape);
template double RAPT::getDelayedSampleAsymmetricHermite1(double d, double *y, double shape);
template double RAPT::getDelayedSampleLinear(double d, double *y);





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


//template class RAPT::rsOnePoleFilter<double, double>;
//template class RAPT::rsOnePoleFilter<rsFloat64x2, double>;
template class RAPT::rsSmoothingFilter<double, double>;

template class RAPT::rsBiquadCascade<double, double>;
template class RAPT::rsBiquadCascade<rsFloat64x2, double>;
template class RAPT::rsDirectFormFilter<double, double>;

template class RAPT::rsFilterAnalyzer<double>;

template class RAPT::rsPrototypeDesigner<double>;
template class RAPT::rsPoleZeroMapper<double>;
template class RAPT::rsFilterCoefficientConverter<double>;
template class RAPT::rsInfiniteImpulseResponseDesigner<double>;
template class RAPT::rsQuadratureNetwork<double, double>;
//template class RAPT::rsLinkwitzRileyCrossOver<double, double>;
//template class RAPT::rsCrossOver4Way<double, double>;
template class RAPT::rsCrossOver4Way<rsFloat64x2, double>;
template class RAPT::rsEngineersFilter<double, double>;
template class RAPT::rsEngineersFilter<rsFloat64x2, double>;

template class RAPT::rsEllipticSubBandFilter<double, double>;
template class RAPT::rsEllipticSubBandFilterDirectForm<double, double>;

template struct RAPT::rsFilterSpecificationZPK<double>;
template struct RAPT::rsFilterSpecificationBA<double>;


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
template class RAPT::rsOnePoleFilter<rsFloat64x2, double>;

template class RAPT::rsModalFilter<double, double>;
template class RAPT::rsNonlinearModalFilter<double, double>;
template class RAPT::rsModalFilterBank<double, double>;
template class RAPT::rsModalFilterWithAttack2<double, double>;

//template class RAPT::rsStateVariableFilter<double, double>;
template class RAPT::rsPhonoFilter<double, double>;
template class RAPT::rsMovingAverage<double, double>;

template class RAPT::rsLadderFilter2<double, double>;
template class RAPT::rsLadderFilterZDF<double, double>;
template class RAPT::rsLadderResoShaped<double, double>;
template class RAPT::rsLadderResoShaped2<double, double>;
template class RAPT::rsLadderFilterFeedbackSaturated<double, double>;
template class RAPT::rsResoReplacer<double, double>;
template class RAPT::rsResoReplacerPhaseBumped<double, double>;
template class RAPT::rsFakeResonanceFilter<double, double>;
template class RAPT::rsLadderMystran<double, double>;

template std::complex<double> RAPT::dcGainNormalizer(
  const std::complex<double>* zeros, size_t numZeros, 
  const std::complex<double>* poles, size_t numPoles);


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
template class RAPT::rsEnvelopeExtractor<double>;
