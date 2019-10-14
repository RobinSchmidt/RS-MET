

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