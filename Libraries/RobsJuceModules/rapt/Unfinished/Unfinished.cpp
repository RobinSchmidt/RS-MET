#include "Unfinished.h"

namespace RAPT
{

#include "Data/MultiArray.cpp"

#include "Math/Matrix.cpp"
#include "Math/FunctionObjects.cpp"
#include "Math/GradientBasedMinimizer.cpp"
#include "Math/MultiLayerPerceptron.cpp"
#include "Math/FourierTransformer.cpp"
#include "Math/NumberTheory.cpp"
#include "Math/Transforms.cpp"
#include "Math/Statistics.cpp"
#include "Math/GeometricFunctions.cpp"
#include "Math/ModularInteger.cpp"
#include "Math/RationalFunctionTools.cpp"
#include "Math/RationalFunction.cpp"

#include "MiscAudio/Interpolator.cpp"
#include "MiscAudio/DelayLine.cpp"
#include "MiscAudio/MiscAudio.cpp"  // may have to be included later (needs higher level stuff)
#include "MiscAudio/Saturator.cpp"
#include "MiscAudio/AudioFunctions.cpp"
#include "MiscAudio/BandwidthConverter.cpp"
#include "MiscAudio/DoublePendulum.cpp"
#include "MiscAudio/ResponseGetters.cpp"
#include "MiscAudio/BlitBlepBlamp.cpp"
#include "MiscAudio/BlepBlampOscs.cpp"
#include "MiscAudio/OscArrays.cpp"

#include "Filters/BandSplitter.cpp"
#include "Filters/FilterDesignFormulas.cpp"
#include "Filters/MovingAverage.cpp"
#include "Filters/PhonoFilter.cpp"
#include "Filters/ModalFilterBank.cpp"
#include "Filters/LadderFilter.cpp"
#include "Filters/FakeResonanceFilter.cpp"
#include "Filters/Biquad.cpp"
#include "Filters/NonUniformFilter.cpp"

// the new polyphony stuff:
#include "MiscAudio/Polyphony.cpp"
#include "MiscAudio/AttackDecayEnvelope.cpp"

// under construction:
#include "Analysis/LinearPredictor.cpp"
#include "Analysis/FormantRemover.cpp"
#include "Analysis/CyclicAutoCorrelator.cpp"
#include "Analysis/AutoCorrelationPitchDetector.cpp"
#include "Analysis/ZeroCrossingPitchDetector.cpp"
//#include "Analysis/ResponseGetters.cpp"

// under construction:
#include "Sampling/MiscUnfinished.cpp"
#include "Sampling/SampleManipulation.cpp"



}