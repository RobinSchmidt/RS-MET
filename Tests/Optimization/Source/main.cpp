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


int main(int argc, char* argv[])
{
  // try to figure out, what difference it makes (if any) to declare function parameters const vs
  // not doing so

  static const int N = 10;
  double a[N], b[N], c[N], d[N];
  for(int n = 0; n < N; n++)
    a[n] = b[n] = n;

  add1(N, a, b, c);
  add2(N, a, b, d);

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