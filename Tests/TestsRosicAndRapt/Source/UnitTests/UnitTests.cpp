#include "UnitTests.h"

// unity build file for the unit tests (todo: split into two files for rosic and rapt - or...for 
// rapt just use the existing file...but maybe rename it to UnitTestsRapt.cpp)

#include "../rapt_tests/UnitTests/DataUnitTests.cpp"
#include "../rapt_tests/UnitTests/FilterUnitTests.cpp"
#include "../rapt_tests/UnitTests/ImageUnitTests.cpp"
#include "../rapt_tests/UnitTests/MathUnitTests.cpp"
#include "../rapt_tests/UnitTests/DrawingUnitTests.cpp"
#include "../rapt_tests/UnitTests/MiscUnitTests.cpp"
#include "../rapt_tests/UnitTests/SortAndSearchTests.cpp"
#include "../rapt_tests/UnitTests/BufferFunctionTests.cpp"

#include "../rapt_tests/UnitTests/UnitTestsRapt.cpp"  // this should become the unity build file for rapt unit tests





#include "../rosic_tests/rosic_AnalysisTests.cpp"
#include "../rosic_tests/rosic_BasicsTests.cpp"
#include "../rosic_tests/rosic_CorrectnessTests.cpp"
#include "../rosic_tests/rosic_EffectsTests.cpp"
#include "../rosic_tests/rosic_FileTests.cpp"
#include "../rosic_tests/rosic_FilterTests.cpp"
#include "../rosic_tests/rosic_GeneratorsTests.cpp"
#include "../rosic_tests/rosic_MathTests.cpp"
#include "../rosic_tests/rosic_ModulatorsTests.cpp"
#include "../rosic_tests/rosic_NonRealTimeTests.cpp"
#include "../rosic_tests/rosic_NumericalTests.cpp"
#include "../rosic_tests/rosic_OthersTests.cpp"
#include "../rosic_tests/rosic_StringTests.cpp"

// get rid:
using namespace RAPT;
using namespace rosic;

#include "../rosic_tests/PortedFromRSLib/UnitTests/FilterTests.cpp"
#include "../rosic_tests/PortedFromRSLib/UnitTests/MiscAudioTests.cpp"
#include "../rosic_tests/PortedFromRSLib/UnitTests/ModalTests.cpp"
#include "../rosic_tests/PortedFromRSLib/UnitTests/NumberManipulationsTests.cpp"
#include "../rosic_tests/PortedFromRSLib/UnitTests/PitchDetectorTests.cpp"
#include "../rosic_tests/PortedFromRSLib/UnitTests/TypeSizeTests.cpp"


//#include "../rosic_tests/UnitTestsRosic.cpp"  // should become unity build file for rosic unit tests