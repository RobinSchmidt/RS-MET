#ifndef rosic_CorrectnessTests_h
#define rosic_CorrectnessTests_h

#include "datastructures/rosic_StringTests.h"
#include "infrastructure/rosic_FileTests.h"
#include "effects/rosic_EffectsTests.h"
#include "filters/rosic_FilterTests.h"
#include "generators/rosic_GeneratorsTests.h"
#include "modulators/rosic_ModulatorsTests.h"
#include "analysis/rosic_AnalysisTests.h"
#include "basics/rosic_BasicsTests.h"
#include "math/rosic_MathTests.h"
#include "numerical/rosic_NumericalTests.h"
#include "nonrealtime/rosic_NonRealtimeTests.h"
#include "others/rosic_OthersTests.h"

namespace rotes
{

  void testAllRosicClasses();
  void testRosicAnalysis(); 
  void testRosicBasics(); 
  void testRosicFile(); 
  void testRosicEffects(); 
  void testRosicFilter(); 
  void testRosicGenerators(); 
  void testRosicModulators(); 
  void testRosicMath(); 
  void testRosicNumerical(); 
  void testRosicNonRealTime(); 
  void testRosicOthers(); 

}

#endif 


