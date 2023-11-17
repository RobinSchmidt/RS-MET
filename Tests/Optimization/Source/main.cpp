#include <vector>

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



// 3 variations of a class has 4 arrays of double as data members. They have wildy different memory
// footprints. The naive implementation SomeArrays1 has hust 4 members of type std::vector. This 
// has a memory footprint of 128 byte. The variant SomeAttays2 stores the length just once and no
// capacity and thereby saves quite a lot of memory. It uses 40 bytes. Finally SomeArrays3 has just
// one pointer and implements the 4 arrays as parts of that allocated memory block. It has the 
// additonal potential advantage that the 4 arrays are guaranteed to be in successive memory
// blocks. And the best: the size would not grow further, if we would need more than 4 arrays.
// Debugging-wise, the std::vector based version is most convenient but when memory should be saved
// then one might opt for variant 3. This could be a good strategy for rsBiquadCascade. It has 
// these a, b, x, y arrays. We could make the objects a lot smaller by that strategy. And biquad 
// cascades are very important and can be used in all sorts of places, so optimizing their memory
// footprint might be a good idea.

class SomeArrays1
{

public:

  SomeArrays1(size_t N) : a(N), b(N), c(N), d(N) {}

  size_t getLength() const { return a.size(); }

  // Read/Write Accessors:
  double* getA() { return &a[0]; }
  double* getB() { return &b[0]; }
  double* getC() { return &c[0]; }
  double* getD() { return &d[0]; }


protected:

  std::vector<double> a, b, c, d;

};

class SomeArrays2
{

public:

  SomeArrays2(size_t initialLength)
  {
    length = initialLength;
    a = new double[length];
    b = new double[length];
    c = new double[length];
    d = new double[length];
  }

  ~SomeArrays2()
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

class SomeArrays3
{
public:

  SomeArrays3(size_t initialLength)
  {
    length = initialLength;
    data = new double[4*length];
  }

  ~SomeArrays3()
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



int main(int argc, char* argv[])
{
  // The compiler spits out a main.asm file into the folder x64/Debug or x64/Release

  // try to figure out, what difference it makes (if any) to declare function parameters const vs
  // not doing so

  static const int N = 10;
  double a[N], b[N], c[N], d[N];
  for(int n = 0; n < N; n++)
    a[n] = b[n] = n;

  add1(N, a, b, c);
  add2(N, a, b, d);

  // Check the memory footprints of the different implementations of a class that has 4 arrays of
  // double as members:
  int sizeV = sizeof(std::vector<double>);  // 32
  int size1 = sizeof(SomeArrays1);          // 128 = 4*32
  int size2 = sizeof(SomeArrays2);          // 40  = 4*8  + 8
  int size3 = sizeof(SomeArrays3);          // 16  = 1*8  + 8, will not grow with number of arrays
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

*/