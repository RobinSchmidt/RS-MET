#ifndef rosic_FilterTests_h
#define rosic_FilterTests_h

//#include "../datastructures/rosic_StringTests.h"

namespace rotes
{
  void testLadderFilter();
  void testModalFilter();
  void testModalFilterWithAttack();
  void testBiquadPhasePlot();
  void testFiniteImpulseResponseDesigner();
  void testConvolverPartitioned();
  void testFiniteImpulseResponseFilter();
  void testFilterAnalyzer();
  void testBiquadCascade();
  void testCrossover4Way();
  void testCrossover4Way2();
  void testSlopeFilter();
  void testPrototypeDesigner();
  void testLowpassToLowshelf();
  void testBesselPrototypeDesign();
  void testPapoulisPrototypeDesign();
  void testEngineersFilter();
  void testPoleZeroMapping();
  void highOrderFilterPolesAndZeros();
}

#endif 