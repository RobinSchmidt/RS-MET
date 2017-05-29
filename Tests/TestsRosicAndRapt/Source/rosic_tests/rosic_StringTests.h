#ifndef rosic_StringTests_h
#define rosic_StringTests_h

//#include "../../../rosic/datastructures/rosic_String.h"
//#include "../../../rosic/rosic.h"
#include "rosic/rosic.h"

namespace rotes
{

  void testRosicString(); // all tests for rosic::String


  //void testStringComparisons();
  void testCharacterComparisons();
  void testStringBufferCopying();
  void testStringIntConversions(int numIterations = 10000); 
  void testStringDoubleConversions(); 

  void testStringDoubleConversionsRandom(int numIterations = 10000); 
  void testStringDoubleConversionsSpecialValues(); 
  void testStringDoubleConversionsDenormals(); 
  void testStringDoubleConversionsLarge(); 
  void testStringDoubleConversionsGeometricProgression(double start, double factor);


  //void testStringConcatenation();
  //void testStringComparison();...


  rosic::rsString createStringWithAllCharacters();
  rosic::rsString createStringWithAllPrintableCharacters();
}

#endif 