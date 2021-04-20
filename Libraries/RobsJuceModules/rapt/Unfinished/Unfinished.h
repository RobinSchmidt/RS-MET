#ifndef RAPT_UNFINISHED_H_INCLUDED
#define RAPT_UNFINISHED_H_INCLUDED

/* Here are classes and functions that are not yet finished and not yet ready for use in production
code. It's a sort of construction yard for new code. */

namespace RAPT
{

#include "Data/Flags.h"
//#include "Data/MultiArray.h"     // obsolete - moved to _Deprecated folder

#include "Math/Vector.h"
#include "Math/Matrix.h"              // this is deprecated
//#include "Math/MatrixUnOptimized.h" // simpler, can perhaps be optimized via move-constructor
//#include "Math/FunctionObjects.h"
//#include "Math/GradientBasedMinimizer.h"
//#include "Math/MultiLayerPerceptron.h"
#include "Math/DifferentialEquationSystem.h"
#include "Math/NumberTheory.h"
#include "Math/Transforms.h"
#include "Math/Statistics.h"  // merge with the Statistics.h aleady present in the not "Unfinished" Math folder
#include "Math/GeometricFunctions.h"
#include "Math/Point2D.h"
#include "Math/Polygon2D.h"
#include "Math/Rectangle2D.h"
#include "Math/Triangle2D.h"
#include "Math/AffineTransform2D.h"
#include "Math/ModularInteger.h"
//#include "Math/RationalFunctionTools.h"
#include "Math/RationalFunction.h"

// still missing math files from RSLib:  BigInt/BigFloat (should go into rosic)

#include "MiscAudio/Interpolator.h"
#include "MiscAudio/DelayLine.h"
#include "MiscAudio/MiscAudio.h"  // may have to be included later (needs higher level stuff)
#include "MiscAudio/Saturator.h"
#include "MiscAudio/AudioFunctions.h" // merge with other AudioFunctions.h file
#include "MiscAudio/BandwidthConverter.h"
#include "MiscAudio/DoublePendulum.h" // will go to Physics
#include "MiscAudio/ResponseGetters.h"
#include "MiscAudio/BlitBlepBlamp.h"
#include "MiscAudio/BlepBlampOscs.h"
#include "MiscAudio/OscArrays.h"

#include "Filters/BandSplitter.h"
#include "Filters/FilterDesignFormulas.h"
#include "Filters/MovingAverage.h"
#include "Filters/PhonoFilter.h"
#include "Filters/ModalFilterBank.h"
#include "Filters/LadderFilter.h"
#include "Filters/FakeResonanceFilter.h"
#include "Filters/Biquad.h"
#include "Filters/NonUniformFilter.h"

// the new polyphony stuff:
#include "MiscAudio/Polyphony.h"
#include "MiscAudio/AttackDecayEnvelope.h"

#include "Analysis/LinearPredictor.h"
#include "Analysis/FormantRemover.h"
#include "Analysis/CyclicAutoCorrelator.h"
#include "Analysis/AutoCorrelationPitchDetector.h"
#include "Analysis/ZeroCrossingPitchDetector.h"
//#include "Analysis/ResponseGetters.h" // maybe that code should be in the test suite

// under construction:
#include "Sampling/MiscUnfinished.h"
#include "Sampling/SampleManipulation.h"

}

#endif