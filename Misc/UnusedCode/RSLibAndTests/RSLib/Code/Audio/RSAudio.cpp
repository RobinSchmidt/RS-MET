// unity-build file for RSLib's Audio module (depends on RSMath.cpp)


#include "RSAudio.h"

#include "Misc/AudioFunctions.cpp"
#include "Misc/BandwidthConverter.cpp"
#include "Misc/MiscAudio.cpp"
#include "Misc/Saturator.cpp"
#include "Misc/SlewRateLimiter.cpp"
#include "Misc/SlewRateLimiterLinear.cpp"
#include "Misc/WindowFunctions.cpp"
#include "Misc/Interpolator.cpp"

#include "Filters/FilterAnalyzer.cpp"
#include "Filters/PrototypeDesigner.cpp"
#include "Filters/PoleZeroMapper.cpp"
#include "Filters/FilterCoefficientConverter.cpp"
#include "Filters/InfiniteImpulseResponseDesigner.cpp"
#include "Filters/Biquad.cpp"
#include "Filters/BiquadCascade.cpp"
#include "Filters/EngineersFilter.cpp"
#include "Filters/FilterDesignFormulas.cpp"
#include "Filters/FourPoleFilter.cpp"
#include "Filters/OnePoleFilter.cpp"
#include "Filters/PhonoFilter.cpp"
#include "Filters/StateVariableFilter.cpp"
#include "Filters/LadderFilter.cpp"
#include "Filters/MovingAverage.cpp"
#include "Filters/ModalFilterBank.cpp"
#include "Filters/FakeResonanceFilter.cpp"

#include "Physics/DoublePendulum.cpp"

#include "Delay/DelayLine.cpp"

#include "Analysis/AutoCorrelationPitchDetector.cpp"
#include "Analysis/CyclicAutoCorrelator.cpp"
#include "Analysis/EnvelopeFollower.cpp"
#include "Analysis/FormantRemover.cpp"
#include "Analysis/LinearPredictor.cpp"
#include "Analysis/ZeroCrossingPitchDetector.cpp"

#include "Modulators/BreakpointModulator.cpp"

#include "Sampling/Resampler.cpp"
#include "Sampling/SampleManipulation.cpp"
#include "Sampling/PhaseVocoder.cpp"
