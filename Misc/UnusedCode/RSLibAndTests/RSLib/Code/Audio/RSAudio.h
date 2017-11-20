#ifndef RS_RSAUDIO_H
#define RS_RSAUDIO_H


#include "../Math/RSMath.h"

#include "Misc/AudioFunctions.h"
#include "Misc/BandwidthConverter.h"
#include "Misc/MiscAudio.h"
#include "Misc/Saturator.h"
#include "Misc/SlewRateLimiter.h"
#include "Misc/SlewRateLimiterLinear.h"
#include "Misc/WindowFunctions.h"
#include "Misc/Interpolator.h"


#include "Delay/DelayLine.h"

#include "Filters/FilterAnalyzer.h"
#include "Filters/FilterDesignFormulas.h"
#include "Filters/PrototypeDesigner.h"
#include "Filters/PoleZeroMapper.h"
#include "Filters/FilterCoefficientConverter.h"
#include "Filters/InfiniteImpulseResponseDesigner.h"
#include "Filters/Biquad.h"
#include "Filters/BiquadCascade.h"
#include "Filters/EngineersFilter.h"
#include "Filters/MovingAverage.h"
#include "Filters/ModalFilterBank.h"
#include "Filters/OnePoleFilter.h"
#include "Filters/FourPoleFilter.h"
#include "Filters/LadderFilter.h"
#include "Filters/StateVariableFilter.h"
#include "Filters/PhonoFilter.h"
#include "Filters/FakeResonanceFilter.h"

#include "Physics/DoublePendulum.h"

#include "Analysis/ResponseGetters.inl"
#include "Analysis/LinearPredictor.h"
#include "Analysis/FormantRemover.h"
#include "Analysis/EnvelopeFollower.h"
#include "Analysis/CyclicAutoCorrelator.h"
#include "Analysis/AutoCorrelationPitchDetector.h"
#include "Analysis/ZeroCrossingPitchDetector.h"

#include "Modulators/BreakpointModulator.h"

#include "Sampling/Resampler.h"
#include "Sampling/SampleManipulation.h"
#include "Sampling/PhaseVocoder.h"


#endif
