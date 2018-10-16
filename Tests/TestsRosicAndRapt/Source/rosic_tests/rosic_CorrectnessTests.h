#ifndef rosic_CorrectnessTests_h
#define rosic_CorrectnessTests_h

#include "rosic_StringTests.h"
#include "rosic_FileTests.h"
#include "rosic_EffectsTests.h"
#include "rosic_FilterTests.h"
#include "rosic_GeneratorsTests.h"
#include "rosic_ModulatorsTests.h"
#include "rosic_AnalysisTests.h"
#include "rosic_BasicsTests.h"
#include "rosic_MathTests.h"
#include "rosic_NumericalTests.h"
#include "rosic_NonRealTimeTests.h"
#include "rosic_OthersTests.h"

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
