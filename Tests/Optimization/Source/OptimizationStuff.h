#pragma once

#include <vector>


//-------------------------------------------------------------------------------------------------
//
// Various variations of a function with 2 parameters and a return value of the same type. We want 
// to figure out, if it makes any difference for the compiled code if we pass the arguments by 
// value, pointer or reference, if constness makes a difference, if using output parameters makes a
// difference etc. We want to inspect the generated assembly code. It ends up in x64/Release when
// compiled in Release mode (which is the relevant mode here)
//
// ToDo: Check, if the compiler inlines the bodies. If so, try to avoid it. We are interested in
// figuring out potential differences in the function call overhead for various ways of passing
// parameters.

double add1(double x, double y)
{
  return x+y;
}

double add2(double& x, double& y)
{
  return x+y;
}

double add3(const double& x, const double& y)
{
  return x+y;
}

double add4(const double* x, const double* y)
{
  return *x + *y;
}

void add5(const double& x, const double& y, double& r)
{
  r = x+y;
}

void add1(int N, double* in1, double* in2, double* out)
{
  for(int n = 0; n < N; n++)
    out[n] = in1[n] + in2[n];
}

void add2(const int& N, const double* in1, const double* in2, double* out)
{
  for(int n = 0; n < N; n++)
    out[n] = in1[n] + in2[n];
}


//-------------------------------------------------------------------------------------------------
//
// 3 variations of a class has 4 arrays of double as data members. They have wildy different memory
// footprints. The naive implementation FourArrays1 has just 4 members of type std::vector. This 
// has a memory footprint of 128 byte. The variant FourArrays2 stores the length just once and no
// capacity and thereby saves quite a lot of memory. It uses 40 bytes. Finally FourArrays3 has just
// one pointer and implements the 4 arrays as parts of that allocated memory block. It has the 
// additonal potential advantage that the 4 arrays are guaranteed to be in successive memory
// blocks. And the best: the size would not grow further, if we would need more than 4 arrays.
// Debugging-wise, the std::vector based version is most convenient but when memory should be saved
// then one might opt for variant 3. This could be a good strategy for rsBiquadCascade. It has 
// these a, b, x, y arrays. We could make the objects a lot smaller by that strategy. And biquad 
// cascades are very important and can be used in all sorts of places, so optimizing their memory
// footprint might be a good idea.

class FourArrays1
{

public:

  FourArrays1(size_t N) : a(N), b(N), c(N), d(N) {}

  size_t getLength() const { return a.size(); }

  // Read/Write Accessors:
  double* getA() { return &a[0]; }
  double* getB() { return &b[0]; }
  double* getC() { return &c[0]; }
  double* getD() { return &d[0]; }


protected:

  std::vector<double> a, b, c, d;

};

class FourArrays2
{

public:

  FourArrays2(size_t initialLength)
  {
    length = initialLength;
    a = new double[length];
    b = new double[length];
    c = new double[length];
    d = new double[length];
  }

  ~FourArrays2()
  {
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
  }

  size_t getLength() const { return length; }

  // Read/Write Accessors:
  double* getA() { return a; }
  double* getB() { return b; }
  double* getC() { return c; }
  double* getD() { return d; }


protected:

  double *a, *b, *c, *d;
  size_t length;

};

class FourArrays3
{
public:

  FourArrays3(size_t initialLength)
  {
    length = initialLength;
    data = new double[4*length];
  }

  ~FourArrays3()
  {
    delete[] data;
  }

  size_t getLength() const { return length; }

  // Read/Write Accessors:
  double* getA() { return data;            }
  double* getB() { return &data[  length]; }
  double* getC() { return &data[2*length]; }
  double* getD() { return &data[3*length]; }


protected:

  double* data;
  size_t  length;

};


//-------------------------------------------------------------------------------------------------
//
// My own idea for implementing a string with short string optimization. I want to try, if it is
// workable.
//
// Oh - nope! it isn't! Nevermind! The idea Was a brainfart!
//
// But I have another idea:
// -use the last byte in the short strinng as follows:
//  -5 of the bits store the length of the short string, the last bit is a flag that indicates 
//   short mode. That implies that capacity must always be even

class MyString
{

public:


  MyString()
  {
    //memset(&str, 0, sizeof(Str));  // Verify!
    // Maybe try to avoid calling memset. Maybe it's more efficient to use
    // str.longStr.data = nullptr;
    // str.longStr.size = 0;
    // str.longStr.capacity = 0;

    // But no! we should init it to the empty string - which is of course a short string!

    str.shortStr.data[0]  = '\0';
    str.shortStr.data[23] = 0b00000001;
  }



  bool is_short() const 
  {  
    return str.shortStr.data[23] & 0b00000001;
    // Alternative: check if str.longStr.capacity is odd
  } 



private:

  static const size_t maxShortLength = 22;
  // The data array for the short string has space for 24 chars. The last char is reserved for 
  // storing the flag for short strings (last bit, i.e. bit index 7 in the byte) and the length of
  // the short string (bits 2..6), i.e. shifting the last bit right by 1 should give the length
  // of the short string
  //
  //  01234567
  //    12345

  struct LongStr
  {
    char*  data;
    size_t size;
    size_t capacity;
  };

  struct ShortStr
  {
    char data[24];
  };

  union Str
  {
    LongStr  longStr;
    ShortStr shortStr;
  };

  Str str;

};

// See:
// https://stackoverflow.com/questions/21694302/what-are-the-mechanics-of-short-string-optimization-in-libc
// Maybe it's pointless to try to reinvent this wheel. It has already been done. std::string should
// be fine.
//
// https://cppdepend.com/blog/understanding-small-string-optimization-sso-in-stdstring/