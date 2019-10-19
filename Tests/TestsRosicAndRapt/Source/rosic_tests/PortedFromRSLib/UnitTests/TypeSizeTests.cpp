//#include "TypeSizeTests.h"

#undef min  // some silly include header on windows (minwindef.h) defines min/max as macros
#undef max

bool testTypeSizes()
{
  bool testResult = true;

  // note: we first assign the sizes to variables (instead of driectly using sizeof inside the comparison) to be able to inspect the actual
  // sizes in the debugger when something goes wrong

  // singed integers:
  int int8Size  = sizeof(rsInt8);
  int int16size = sizeof(rsInt16);
  int int32Size = sizeof(rsInt32);
  int int64size = sizeof(rsInt64);
  testResult &= ( int8Size  ==  1 );
  testResult &= ( int16size ==  2 );
  testResult &= ( int32Size ==  4 );
  testResult &= ( int64size ==  8 );

  // unsigned integers:
  int uint8Size  = sizeof(rsUint8);
  int uint16size = sizeof(rsUint16);
  int uint32Size = sizeof(rsUint32);
  int uint64size = sizeof(rsUint64);
  testResult &= ( uint8Size  ==  1 );
  testResult &= ( uint16size ==  2 );
  testResult &= ( uint32Size ==  4 );
  testResult &= ( uint64size ==  8 );

  // floats:
  int float32Size = sizeof(rsFloat32);
  int float64size = sizeof(rsFloat64);
  //int float80size = sizeof(float80);  // 12 on Win32, 16 on Linux x86
  testResult &= ( float32Size ==  4 );
  testResult &= ( float64size ==  8 );
  //testResult &= ( float80size == 12 );


  //
  rsInt32 rsInt32Min = RS_MIN(rsInt32);  // -2147483648
  rsInt32 rsInt32Max = RS_MAX(rsInt32);  //  2147483647
  testResult &= rsInt32Min == -(int)2147483648; // msvc gives a warning here - why?
  testResult &= rsInt32Max ==  2147483647;

  //double doubleMin = RS_MIN(double);

  return testResult;
}





