#ifdef RS_TESTING_H_INCLUDED
/* When you add this cpp file to your project, you mustn't include it in a file where you've
already included any other headers - just put it inside a file on its own, possibly with your config
flags preceding it, but don't include anything else. That also includes avoiding any automatic prefix
header files that the compiler may be using. */
#error "Incorrect use of JUCE cpp file"
#endif

#include "rs_testing.h"

#include "TestTools/Utilities/FileWriting.cpp"
#include "TestTools/Utilities/PerformanceTestTools.cpp"
#include "TestTools/Utilities/TestInputCreation.cpp"
#include "TestTools/Utilities/TestUtilities.cpp"
#include "Prototypes/OscDrivers.cpp"
#include "Prototypes/PartialDifferentialEquations.cpp"
#include "Prototypes/Prototypes.cpp"
#include "Prototypes/FiniteAutomaton.cpp"
#include "Prototypes/ParticleBouncer.cpp"
#include "Prototypes/ParticleSystem.cpp"
#include "Prototypes/Tensor.cpp"

#include "TestTools/DSPPlotters.cpp"
#include "TestTools/Plotting.cpp"

#include "Experiments/MiscExperiments.cpp"

#include "Misc/TestClasses.cpp"



#include "RaptInstantiations.cpp"

// move these instantioans into RaptInstantiations.cpp:

template struct RAPT::rsFilterSpecificationZPK<float>;
template struct RAPT::rsFilterSpecificationBA<float>;

template class rsTableLinBlep<double, double>;
template class rsTableMinBlep<double, double>;
template class rsBlepReadyOsc<double>;

template class rsRatioGenerator<double>;
template class rsBlepOscArray<double, rsBlepReadyOscBase<double>, rsPolyBlep1<double, double>>;
template class rsBlepOscArray<double, rsBlepReadyOscBase<double>, rsPolyBlep2<double, double>>;
template class rsBlepOscArray<double, rsBlepReadyOscBase<double>, rsTableMinBlep<double, double>>;

template class rsSyncPhasor<double, rsPolyBlep1<double, double>>;
template class rsSyncPhasor<double, rsPolyBlep2<double, double>>;
template class rsSyncPhasor<double, rsTableMinBlep<double, double>>;

template class rsSyncOsc<double, rsPolyBlep2<double, double>>;
template class rsSyncOsc<double, rsTableMinBlep<double, double>>;

template class rsDualBlepOsc<double, rsPolyBlep2<double, double>>;
template class rsDualBlepOsc<double, rsTableMinBlep<double, double>>;


template class rsHeatEquation1D<double>;