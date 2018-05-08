#include "Unfinished.h"

namespace RAPT
{

#include "Data/MultiArray.cpp"

#include "Math/Matrix.cpp"
#include "Math/FunctionObjects.cpp"
#include "Math/GradientBasedMinimizer.cpp"
#include "Math/MultiLayerPerceptron.cpp"
#include "Math/FourierTransformer.cpp"
#include "Math/Interpolation.cpp"
#include "Math/NumberTheory.cpp"
#include "Math/NumericCalculus.cpp"
#include "Math/Transforms.cpp"
#include "Math/Statistics.cpp"
#include "Math/GeometricFunctions.cpp"
#include "Math/ModularInteger.cpp"

#include "MiscAudio/Interpolator.cpp"
#include "MiscAudio/DelayLine.cpp"
#include "MiscAudio/MiscAudio.cpp"  // may have to be included later (needs higher level stuff)
#include "MiscAudio/Saturator.cpp"
#include "MiscAudio/SlewRateLimiter.cpp"
#include "MiscAudio/SlewRateLimiterLinear.cpp"
#include "MiscAudio/WindowFunctions.cpp"
#include "MiscAudio/AudioFunctions.cpp"
#include "MiscAudio/BandwidthConverter.cpp"
#include "MiscAudio/DoublePendulum.cpp"
#include "MiscAudio/ResponseGetters.cpp"

#include "Filters/BandSplitter.cpp"
#include "Filters/FilterDesignFormulas.cpp"
#include "Filters/MovingAverage.cpp"
#include "Filters/PhonoFilter.cpp"
#include "Filters/ModalFilterBank.cpp"
#include "Filters/LadderFilter.cpp"
#include "Filters/FakeResonanceFilter.cpp"
#include "Filters/Biquad.cpp"

// under construction:
#include "Analysis/EnvelopeFollower.cpp" // file is actually empty anyway
#include "Analysis/LinearPredictor.cpp"
#include "Analysis/FormantRemover.cpp"
#include "Analysis/CyclicAutoCorrelator.cpp"
#include "Analysis/AutoCorrelationPitchDetector.cpp"
#include "Analysis/ZeroCrossingPitchDetector.cpp"
//#include "Analysis/ResponseGetters.cpp"

// under construction:
#include "Sampling/PhaseVocoder.cpp"
#include "Sampling/Resampler.cpp"
#include "Sampling/SampleManipulation.cpp"



}