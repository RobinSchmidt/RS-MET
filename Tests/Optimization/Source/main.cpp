/** This project is for experimenting a bit with various optimization techniques and tricks. */


#include "OptimizationStuff.h"




int main(int argc, char* argv[])
{
  // The compiler spits out a main.asm file into the folder x64/Debug or x64/Release for 
  // inspection.
  //
  // ToDo:
  // -Try to figure out, what difference it makes (if any) to declare function parameters const vs
  //  not doing so

  static const int N = 10;
  double a[N], b[N], c[N], d[N];
  for(int n = 0; n < N; n++)
    a[n] = b[n] = n;

  add1(N, a, b, c);
  add2(N, a, b, d);

  // Check the memory footprints of the different implementations of a class that has 4 arrays of
  // double as members:
  int sizeV = sizeof(std::vector<double>);  // 32
  int size1 = sizeof(FourArrays1);          // 128 = 4*32
  int size2 = sizeof(FourArrays2);          // 40  = 4*8  + 8
  int size3 = sizeof(FourArrays3);          // 16  = 1*8  + 8, will not grow with number of arrays
  // Clearly, std::vector loses by a big margin. On the plus side, the implmentation will be easier
  // to read and debug. Generally, a vector/array needs to store: pointer,size,capacity or 
  // start,end,capacity. But that would be 3 bytes, but it apparently has 4 so it must store a 
  // fourth value. In the standard library of Visual Studio,  std::vector seems to have a data 
  // member of tye _Complressed_pair:
  //
  //   _Compressed_pair<_Alty, _Scary_val> _Mypair;
  //
  // and I have no idea what that is and why one would implement a vector liek that. I guess, the
  // code is generated. It certainly doesn't look like it was written with readability in mind.

  // Test my own idea for a string class with short string optimization:
  MyString str;



  return 0;
  //return(EXIT_SUCCESS);
}

// maybe test on: https://godbolt.org/
// this is really useful to see the assembly code generated from short code snippets - here's a 
// video about it: 
// https://www.youtube.com/watch?v=kIoZDUd5DKw
// https://www.youtube.com/watch?v=4_HL3PH4wDg
// https://www.youtube.com/watch?v=1u_ku_OJPDg

/*
Resources:

https://www.youtube.com/watch?v=o4-CwDo2zpg Fastware - Andrei Alexandrescu
https://www.youtube.com/watch?v=Qq_WaiwzOtI CppCon 2014: Andrei Alexandrescu "Optimization Tips - Mo' Hustle Mo' Problems"


This has a function to convert from double to int faster than the usual way:
https://stackoverflow.com/questions/28668348/how-expensive-is-it-to-convert-between-int-and-double


Fast and Small C++ - When Efficiency Matters - Andreas Fertig - CppCon 2024
https://www.youtube.com/watch?v=rNl591__9zY


-My idea for short string optimization: In short mode, use the last byte of the 24 for the length
 of the short string. When subtracting also the space for the zero-terminator, we may store short
 strings up to a length of 22 characters, if this wokrs out:

class rsString
{


  bool is_long() const { str.longStr.data != this; }  // would this work? ...NOPE!

  size_t size() const
  {
     if(is_long()
       return str.longStr.size;
     else
       return (size_t) str.shortStr.data[23];
  }

  size_t capacity()
  {
     if(is_long()
       return str.longStr.capacity;
     else
       return 22;
  }

 

private:

  struct LongStr
  {
    char*  data;
    size_t size;
    size_t capacity;
  };

  struct ShortStr
  {
    char data[24];
  }

  union Str
  {
    LongStr longStr;
    ShortStr shortStr;
  }


  Str str;

}

*/