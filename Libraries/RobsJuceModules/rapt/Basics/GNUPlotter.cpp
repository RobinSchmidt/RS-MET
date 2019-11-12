#include "GNUPlotter.h"
#include <fstream>
#include <iostream>
#include <algorithm>
using namespace std;

GNUPlotter::GNUPlotter()
{
  setDataPrecision(8);

  graphStyles.resize(1);
  graphStyles[0] = std::string("lines");    // standard lines, width 1 - maybe change to 1.5 or 2

  // installation path of the GNUPlot executable:
  //gnuplotPath = "C:/Program Files/gnuplot/bin/gnuplot.exe";
  //gnuplotPath = "C:/Program Files/gnuplot/bin/wgnuplot.exe";
  gnuplotPath = "C:/Octave/Octave-5.1.0.0/mingw64/bin/gnuplot.exe";

  // paths for data file and command batchfile:
  dataPath    = "E:/Temp/gnuplotData.dat";
  commandPath = "E:/Temp/gnuplotCommands.txt";   // this path may not contain whitepaces

  initialize();                                  // initializes data- and commandfile
}

GNUPlotter::~GNUPlotter()
{
  delete[] formatString;
}

void GNUPlotter::initialize()
{
  initFile(dataPath);
  initFile(commandPath);
  addDefaultCommands();
}

// convenience functions:

template<class T>
void GNUPlotter::plotComplexArrayReIm(const T* x, const std::complex<T>* z, int N)
{
  std::vector<T> re(N), im(N);
  for(int i = 0; i < N; i++) {
    re[i] = z[i].real();
    im[i] = z[i].imag();
  }
  GNUPlotter plt;
  plt.addDataArrays(N, x, &re[0]);
  plt.addDataArrays(N, x, &im[0]);
  plt.plot();
}

template<class T>
void GNUPlotter::plotComplexArrayReIm(const std::complex<T>* z, int N)
{
  std::vector<T> x(N);
  rangeLinear(&x[0], N, T(0), T(N-1));
  plotComplexArrayReIm(&x[0], &z[0], N);
}
template void GNUPlotter::plotComplexArrayReIm(const std::complex<int>* z, int N);
template void GNUPlotter::plotComplexArrayReIm(const std::complex<float>* z, int N);
template void GNUPlotter::plotComplexArrayReIm(const std::complex<double>* z, int N);

template <class T>
void GNUPlotter::plotCurve2D(const std::function<T(T)>& fx, const std::function<T(T)>& fy,
  int Nt, T tMin, T tMax)
{
  GNUPlotter p;
  p.addDataCurve2D(fx, fy, Nt, tMin, tMax);
  p.plot();
}
// explicit instantiations for int, float and double:
template void GNUPlotter::plotCurve2D(const std::function<int(int)>& fx, const std::function<int(int)>& fy, int Nt, int tMin, int tMax);
template void GNUPlotter::plotCurve2D(const std::function<float(float)>& fx, const std::function<float(float)>& fy, int Nt, float tMin, float tMax);
template void GNUPlotter::plotCurve2D(const std::function<double(double)>& fx, const std::function<double(double)>& fy, int Nt, double tMin, double tMax);

template <class T>
void GNUPlotter::plotSurface(
  const std::function<T(T, T)>& fx,
  const std::function<T(T, T)>& fy,
  const std::function<T(T, T)>& fz,
  int Nu, T uMin, T uMax, int Nv, T vMin, T vMax)
{
  GNUPlotter p;
  p.addDataSurface(fx, fy, fz, Nu, uMin, uMax, Nv, vMin, vMax);
  p.addCommand("set hidden3d");                  // don't draw hidden lines
  //p.addCommand("set view 20,50");                // set up perspective
  //p.addCommand("set lmargin 0");                 // margin between plot and left border
  //p.addCommand("set tmargin 0");                 // margin between plot and top border
  //p.addCommand("set ztics 0.5");                 // density of z-axis tics
  p.plot3D();                                    // invoke GNUPlot
}
template void GNUPlotter::plotSurface(const std::function<int(int, int)>& fx, const std::function<int(int, int)>& fy, const std::function<int(int, int)>& fz, int Nu, int uMin, int uMax, int Nv, int vMin, int vMax);
template void GNUPlotter::plotSurface(const std::function<float(float, float)>& fx, const std::function<float(float, float)>& fy, const std::function<float(float, float)>& fz, int Nu, float uMin, float uMax, int Nv, float vMin, float vMax);
template void GNUPlotter::plotSurface(const std::function<double(double, double)>& fx, const std::function<double(double, double)>& fy, const std::function<double(double, double)>& fz, int Nu, double uMin, double uMax, int Nv, double vMin, double vMax);

template<class T>
void GNUPlotter::plotVectorField2D(const function<T(T, T)>& fx, const function<T(T, T)>& fy,
  int Nx, T xMin, T xMax, int Ny, T yMin, T yMax)
{
  GNUPlotter p;
  p.addVectorField2D(fx, fy, Nx, xMin, xMax, Ny, yMin, yMax);
  p.addCommand("set palette rgbformulae 30,31,32 negative");
  p.plot();
  // maybe try a black background (and invert the colormap)
}
template void GNUPlotter::plotVectorField2D(const function<int(int, int)>& fx, const function<int(int, int)>& fy, int Nx, int xMin, int xMax, int Ny, int yMin, int yMax);
template void GNUPlotter::plotVectorField2D(const function<float(float, float)>& fx, const function<float(float, float)>& fy, int Nx, float xMin, float xMax, int Ny, float yMin, float yMax);
template void GNUPlotter::plotVectorField2D(const function<double(double, double)>& fx, const function<double(double, double)>& fy, int Nx, double xMin, double xMax, int Ny, double yMin, double yMax);


template<class T>
void GNUPlotter::plotComplexVectorField(const function<complex<T>(complex<T>)>& f,
  int Nr, T rMin, T rMax, int Ni, T iMin, T iMax, bool conj)
{
  T sign = T(1); if(conj) sign = T(-1);
  std::function<T(T, T)> fx, fy;
  fx = [&] (T re, T im) { return        real(f(complex<T>(re, im))); };
  fy = [&] (T re, T im) { return sign * imag(f(complex<T>(re, im))); };
  plotVectorField2D(fx, fy, Nr, rMin, rMax, Ni, iMin, iMax);
}
template void GNUPlotter::plotComplexVectorField(const function<complex<int>(complex<int>)>& f, int Nr, int rMin, int rMax, int Ni, int iMin, int iMax, bool conj);
template void GNUPlotter::plotComplexVectorField(const function<complex<float>(complex<float>)>& f, int Nr, float rMin, float rMax, int Ni, float iMin, float iMax, bool conj);
template void GNUPlotter::plotComplexVectorField(const function<complex<double>(complex<double>)>& f, int Nr, double rMin, double rMax, int Ni, double iMin, double iMax, bool conj);



// plotting:


/*
template <class T>
static void GNUPlotter::plot(int N, T *x, T *y1, T *y2, T *y3, T *y4, T *y5, T *y6, T *y7, T *y8,
  T *y9)
{
  GNUPlotter plt;
  plt.addDataArrays(N, x, y1, y2, y3, y4, y5, y6, y7, y8, y9);
  plt.plot();
}
template void GNUPlotter::plot(int N, double *y1, double *y2, double *y3, double *y4,
  double *y5, double *y6, double *y7, double *y8, double *y9);
*/

void GNUPlotter::plot()
{
  addPlotCommand(false);
  invokeGNUPlot();
}

void GNUPlotter::plot3D()
{
  addPlotCommand(true);
  invokeGNUPlot();
}

template <class T>
void GNUPlotter::plotFunctionTables(int N, const T *x, const T *y1, const T *y2, const T *y3, 
  const T *y4, const T *y5, const T *y6, const T *y7, const T *y8, const T *y9)
{
  addDataArrays(N, x, y1, y2, y3, y4, y5, y6, y7, y8, y9);
  plot();
}

template <class T>
void GNUPlotter::plotArrays(int N, const T *y1, const T *y2, const T *y3, const T *y4, 
  const T *y5, const T *y6, const T *y7, const T *y8, const T *y9)
{
  T *x = new T[N];
  rangeLinear(x, N, T(0), T(N-1));
  plotFunctionTables(N, x, y1, y2, y3, y4, y5, y6, y7, y8, y9);
  delete[] x;
}
// explicit instantiations for double, float and int:
template void GNUPlotter::plotArrays(int N, const double *y1, const double *y2, const double *y3, 
  const double *y4, const double *y5, const double *y6, const double *y7, const double *y8, 
  const double *y9);
template void GNUPlotter::plotArrays(int N, const float *y1, const float *y2, const float *y3, 
  const float *y4, const float *y5, const float *y6, const float *y7, const float *y8, 
  const float *y9);
template void GNUPlotter::plotArrays(int N, const int *y1, const int *y2, const int *y3, 
  const int *y4, const int *y5, const int *y6, const int *y7, const int *y8, const int *y9);

template <class T>
void GNUPlotter::plotFunctions(int N, T *x, T (*f0)(T), T (*f1)(T), T (*f2)(T), T (*f3)(T),
  T (*f4)(T), T (*f5)(T), T (*f6)(T), T (*f7)(T), T (*f8)(T), T (*f9)(T))
{
  addDataFunctions(N, x, f0, f1, f2, f3, f4, f5, f6, f7, f8, f9);
  plot();
}

template <class T>
void GNUPlotter::plotFunctions(int N, T xMin, T xMax, T (*f0)(T), T (*f1)(T), T (*f2)(T),
  T (*f3)(T), T (*f4)(T), T (*f5)(T), T (*f6)(T), T (*f7)(T), T (*f8)(T), T (*f9)(T))
{
  T *x = new T[N];
  rangeLinear(x, N, xMin, xMax);
  plotFunctions(N, x, f0, f1, f2, f3, f4, f5, f6, f7, f8, f9);
  delete[] x;
}
template void GNUPlotter::plotFunctions(int N, double xMin, double xMax, double (*f0)(double),
  double (*f1)(double), double (*f2)(double), double (*f3)(double), double (*f4)(double),
  double (*f5)(double), double (*f6)(double), double (*f7)(double), double (*f8)(double),
  double (*f9)(double));
template void GNUPlotter::plotFunctions(int N, float xMin, float xMax, float (*f0)(float),
  float (*f1)(float), float (*f2)(float), float (*f3)(float), float (*f4)(float),
  float (*f5)(float), float (*f6)(float), float (*f7)(float), float (*f8)(float),
  float (*f9)(float));
template void GNUPlotter::plotFunctions(int N, int xMin, int xMax, int (*f0)(int),
  int (*f1)(int), int (*f2)(int), int (*f3)(int), int (*f4)(int),
  int (*f5)(int), int (*f6)(int), int (*f7)(int), int (*f8)(int),
  int (*f9)(int));

template <class T>
void GNUPlotter::plotSurface(int Nx, int Ny, T *x, T *y, T **z)
{
  addDataMatrix(Nx, Ny, x, y, z);
  plot3D();
}
template void GNUPlotter::plotSurface(int Nx, int Ny, int *x, int *y, int **z);
template void GNUPlotter::plotSurface(int Nx, int Ny, float *x, float *y, float **z);
template void GNUPlotter::plotSurface(int Nx, int Ny, double *x, double *y, double **z);

template <class T>
void GNUPlotter::plotBivariateFunction(int Nx, int Ny, T *x, T *y, T (*f)(T, T))
{
  addDataBivariateFunction(Nx, Ny, x, y, f);
  plot3D();

  //int i, j;
  //T **z = new T*[Nx];
  //for(i = 0; i < Nx; i++)
  //  z[i] = new T[Ny];
  //for(i = 0; i < Nx; i++)
  //{
  //  for(j = 0; j < Ny; j++)
  //    z[i][j] = f(x[i], y[j]);
  //}
  //plotSurface(Nx, Ny, x, y, z);
  //for(i = 0; i < Nx; i++)
  //  delete[] z[i];
  //delete[] z;
}

template <class T>
void GNUPlotter::plotBivariateFunction(int Nx, T xMin, T xMax, int Ny, T yMin, T yMax,
  T (*f)(T, T))
{
  addDataBivariateFunction(Nx, xMin, xMax, Ny, yMin, yMax, f);
  plot3D();

  //T *x = new T[Nx];
  //T *y = new T[Ny];
  //GNUPlotter::rangeLinear(x, Nx, xMin, xMax);
  //GNUPlotter::rangeLinear(y, Ny, yMin, yMax);
  //plotBivariateFunction(Nx, Ny, x, y, f);
  //delete[] x;
  //delete[] y;
}
template void GNUPlotter::plotBivariateFunction(int Nx, double xMin, double xMax, int Ny,
  double yMin, double yMax, double (*f)(double, double));
template void GNUPlotter::plotBivariateFunction(int Nx, float xMin, float xMax, int Ny, float yMin,
  float yMax, float (*f)(float, float));
template void GNUPlotter::plotBivariateFunction(int Nx, int xMin, int xMax, int Ny, int yMin,
  int yMax, int (*f)(int, int));

// style setup:

void GNUPlotter::addCommand(string command)
{
  // factor out into a function "withNewLine":
  const char *n = "\n";
  if(command.back() != *n)
    command += "\n";

  std::ofstream out(commandPath, std::ofstream::app);    // open commandfile for append
  out << command;  // later out << withNewLine(command)
  out.close();
}

void GNUPlotter::addDefaultCommands()
{
  addCommand("# Default Settings:");
  addCommand("set zero 0"); // by default, GNUPlot uses 1.e-8 as threshold - numbers smaller than
                            // that in absolute value may be interpreted as zero and not plotted
                            // properly
  setGrid();
  setGraphColors("000000", "0000FF", "007700", "BB0000", "9900FF");
  //addCommand("set linetype cycle 2"); // doesn't work
  addCommand("\n# Custom Settings:"); // subsequent commands appear in "Custom Settings" section
}

void GNUPlotter::setAxisLabels(std::string x, std::string y, std::string z)
{
  if( !x.empty() ) addCommand("set xlabel \"" + x + "\"\n");
  if( !y.empty() ) addCommand("set ylabel \"" + y + "\"\n");
  if( !z.empty() ) addCommand("set zlabel \"" + z + "\"\n");
}

void GNUPlotter::setLegends(CSR l0, CSR l1, CSR l2, CSR l3, CSR l4, CSR l5, CSR l6, CSR l7,
  CSR l8, CSR l9)
{
  addCommand("set key opaque box"); // graphs shouldn't obscure legends
  setStringVector(graphTitles, l0, l1, l2, l3, l4, l5, l6, l7, l8, l9);
}

void GNUPlotter::setLegends(CVR<string> legends)
{
  graphTitles = legends;
}

void GNUPlotter::setGraphColors(CSR c0, CSR c1, CSR c2, CSR c3, CSR c4, CSR c5, CSR c6, CSR c7,
  CSR c8, CSR c9)
{
  vector<string> v;
  setStringVector(v, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
  setGraphColors(v);
}

void GNUPlotter::setGraphColors(CVR<string> c)
{
  for(unsigned int i = 0; i < c.size(); i++)
    setGraphColor(i+1, c[i]);
}

void GNUPlotter::setGraphColor(unsigned int index, CSR color)
{
  addCommand("set lt " + s(index) + " lc rgb \"#" + color + "\"" ); // lt: linetype, lc: linecolor
}

void GNUPlotter::setDashType(unsigned int index, CSR type)
{
  addCommand("set lt " + s(index) + " dt " + type); // lt: linetype, dt: dashtype
}

void GNUPlotter::setGraphStyles(CSR s0, CSR s1, CSR s2, CSR s3, CSR s4, CSR s5, CSR s6, CSR s7,
  CSR s8, CSR s9)
{
  setStringVector(graphStyles, s0, s1, s2, s3, s4, s5, s6, s7, s8, s9);
}

void GNUPlotter::setGraphStyles(CVR<string> styles)
{
  graphStyles = styles;
}

void GNUPlotter::setGrid(bool x, bool y, bool x2, bool y2, bool z)
{
  std::string s = "set grid ";
  if(!x)  s += "no"; s += "xtics ";
  if(!y)  s += "no"; s += "ytics ";
  if(!z)  s += "no"; s += "ztics ";
  if(!x2) s += "no"; s += "x2tics ";
  if(!y2) s += "no"; s += "y2tics \n";
  addCommand(s);
}

void GNUPlotter::setLogScale(string axes, double base, bool shouldBeLogarithmic)
{
  string s;
  if( !shouldBeLogarithmic )
    s += "un";
  s += "set logscale " + axes + "\n";
  addCommand(s);
}

void GNUPlotter::setRange(double xMin, double xMax, double yMin, double yMax, double zMin,
  double zMax)
{
  if( xMin < xMax ) addCommand("set xrange[" + s(xMin) + ":" + s(xMax) + "]\n");
  if( yMin < yMax ) addCommand("set yrange[" + s(yMin) + ":" + s(yMax) + "]\n");
  if( zMin < zMax ) addCommand("set zrange[" + s(zMin) + ":" + s(zMax) + "]\n");
}

void GNUPlotter::setPixelSize(unsigned int width, unsigned int height)
{
  addCommand("set terminal wxt size " + s(width) +  "," + s(height) + "\n");
}

void GNUPlotter::setTitle(std::string title)
{
  if( !title.empty()  ) addCommand("set title \""  + title  + "\"\n");
}

void GNUPlotter::addAnnotation(double x, double y, CSR text, CSR options)
{
  addCommand("set label at " + s(x) + "," + s(y) + " \"" + text + "\" " + options + "\n");
}

// data setup:

void GNUPlotter::setDataPrecision(unsigned int n)
{
  delete[] formatString;
  std::string str = "% 0" + s(n) + "." + s(n) + "le";
  formatString = toZeroTerminatedString(str);
}

template <class T>
void GNUPlotter::addDataBlockLineColumn(const vector<vector<vector<T>>>& d)
{
  ofstream out(dataPath, ofstream::app);
  for(size_t i = 0; i < d.size(); i++) {           // loop over blocks
    for(size_t j = 0; j < d[i].size(); j++) {      // loop over rows in current block
      for(size_t k = 0; k < d[i][j].size(); k++)   // loop over columns in current row
        out << sd(d[i][j][k]) + " ";
      out << "\n";
    }
    out << "\n";
  }
  out << "\n";
  out.close();
  dataInfo.push_back(DataInfo(d.size(), d[0].size())); // keep track of written data
}
// explicit instantiations for double, float and int:
template void GNUPlotter::addDataBlockLineColumn(const vector<vector<vector<double>>>& d);
template void GNUPlotter::addDataBlockLineColumn(const vector<vector<vector<float>>>& d);
template void GNUPlotter::addDataBlockLineColumn(const vector<vector<vector<int>>>& d);

template <class T>
void GNUPlotter::addDataBlockColumnLine(const std::vector<std::vector<std::vector<T>>>& d)
{
  ofstream out(dataPath, ofstream::app);
  for(size_t i = 0; i < d.size(); i++) {           // loop over blocks
    for(size_t j = 0; j < d[i][0].size(); j++) {   // loop over columns
      for(size_t k = 0; k < d[i].size(); k++)      // loop over lines
        out << sd(d[i][k][j]) + " ";
      out << "\n";
    }
    out << "\n";
  }
  out << "\n";
  out.close();
  dataInfo.push_back(DataInfo(d.size(), d[0][0].size()));
}
template void GNUPlotter::addDataBlockColumnLine(const vector<vector<vector<double>>>& d);
template void GNUPlotter::addDataBlockColumnLine(const vector<vector<vector<float>>>& d);
template void GNUPlotter::addDataBlockColumnLine(const vector<vector<vector<int>>>& d);

template <class T>
void GNUPlotter::addData(int numBlocks, int *blockLengths, int numColumns, T **data)
{
  ofstream out(dataPath, ofstream::app);
  int offset = 0;
  for(int i = 0; i < numBlocks; i++)           // loop over blocks
  {
    for(int j = 0; j < blockLengths[i]; j++)   // loop over rows in current block
    {
      for(int k = 0; k < numColumns; k++)      // loop over columns in current row
        out << sd(data[k][j+offset]) + " ";
      out << "\n";
    }
    offset += blockLengths[i];
    out << "\n";
  }
  out << "\n";
  out.close();
  dataInfo.push_back(DataInfo(numBlocks, numColumns)); // keep track of written data
}

template <class T>
void GNUPlotter::addData(int numRows, int numColumns, T **data)
{
  addData(1, &numRows, numColumns, data);
}

template <class T>
void GNUPlotter::addDataComplex(const vector<complex<T>>& d)
{
  ofstream out(dataPath, ofstream::app);
  for(size_t i = 0; i < d.size(); i++)
  {
    out << sd(d[i].real()) + " " + sd(d[i].imag()) + " ";
    out << "\n";
  }
  out << "\n\n";
  out.close();
  dataInfo.push_back(DataInfo(1, 2));
}
template void GNUPlotter::addDataComplex(const vector<complex<int>>& d);
template void GNUPlotter::addDataComplex(const vector<complex<float>>& d);
template void GNUPlotter::addDataComplex(const vector<complex<double>>& d);

template <class T>
void GNUPlotter::addDataArrays(int N, T *x, int M, T **y)
{
  ofstream out(dataPath, ofstream::app);
  for(int i = 0; i < N; i++)
  {
    out << sd(x[i]) + " ";
    for(int j = 0; j < M; j++)
      out << sd(y[j][i]) + " ";
    out << "\n";
  }
  out << "\n\n";
  out.close();
  dataInfo.push_back(DataInfo(1, M+1));
}
template void GNUPlotter::addDataArrays(int N, int *x, int M, int **y);
template void GNUPlotter::addDataArrays(int N, float *x, int M, float **y);
template void GNUPlotter::addDataArrays(int N, double *x, int M, double **y);

template <class T>
void GNUPlotter::addDataArrays(int N, const T *c0, const T *c1, const T *c2, const T *c3, 
  const T *c4, const T *c5, const T *c6, const T *c7, const T *c8, const T *c9)
{
  const vector<const T*> v = collectLeadingNonNullArguments(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
  const T* a[10];
  for(size_t i = 0; i < v.size(); i++)
    a[i] = v[i];
  addData(N, (int) v.size(), a);
}

template <class T>
void GNUPlotter::addDataFunctions(int N, T *x, T (*f0)(T), T (*f1)(T), T (*f2)(T), T (*f3)(T),
  T (*f4)(T), T (*f5)(T), T (*f6)(T), T (*f7)(T), T (*f8)(T), T (*f9)(T))
{
  int i, j;
  vector<T(*)(T)> f = collectLeadingNonNullArguments(f0, f1, f2, f3, f4, f5, f6, f7, f8, f9);

  // allocate data matrix (1st index is column, 2nd ist row):
  int M = (int)f.size()+1;
  T **d = new T*[M];
  for(i = 0; i < M; i++)
    d[i] = new T[N];

  // copy x-axis into 1st column:
  for(j = 0; j < N; j++)
    d[0][j] = x[j];

  // use the M-1 other columns for the different functions:
  for(i = 1; i < M; i++)
  {
    for(j = 0; j < N; j++)
      d[i][j] = f[i-1](x[j]); // d[i][j] = apply function pointer f[i-1] to x[j]
  }

  // add dataset and clean up memory:
  addData(N, M, d);
  for(i = 0; i < M; i++)
    delete[] d[i];
  delete[] d;
}

template <class T>
void GNUPlotter::addDataFunctions(int N, T xMin, T xMax, T (*f0)(T), T (*f1)(T), T (*f2)(T),
  T (*f3)(T), T (*f4)(T), T (*f5)(T), T (*f6)(T), T (*f7)(T), T (*f8)(T), T (*f9)(T))
{
  T *x = new T[N];
  rangeLinear(x, N, xMin, xMax);
  addDataFunctions(N, x, f0, f1, f2, f3, f4, f5, f6, f7, f8, f9);
  delete[] x;
}
template void GNUPlotter::addDataFunctions(int N, double xMin, double xMax,
  double (*f0)(double), double (*f1)(double), double (*f2)(double), double (*f3)(double),
  double (*f4)(double), double (*f5)(double), double (*f6)(double), double (*f7)(double),
  double (*f8)(double), double (*f9)(double));
template void GNUPlotter::addDataFunctions(int N, float xMin, float xMax,
  float (*f0)(float), float (*f1)(float), float (*f2)(float), float (*f3)(float),
  float (*f4)(float), float (*f5)(float), float (*f6)(float), float (*f7)(float),
  float (*f8)(float), float (*f9)(float));
template void GNUPlotter::addDataFunctions(int N, int xMin, int xMax,
  int (*f0)(int), int (*f1)(int), int (*f2)(int), int (*f3)(int),
  int (*f4)(int), int (*f5)(int), int (*f6)(int), int (*f7)(int),
  int (*f8)(int), int (*f9)(int));

template <class T>
void GNUPlotter::addDataGrid(int Nx, int Ny, T *x, T *y, T **z)
{
  ofstream out(dataPath, ofstream::app);
  for(int i = 0; i < Nx; i++)
  {
    for(int j = 0; j < Ny; j++)
      out << sd(x[i]) + " " + sd(y[j]) + " " + sd(z[i][j]) + "\n";
    out << "\n";
  }
  out << "\n";
  out.close();
  dataInfo.push_back(DataInfo(Nx, 3));
}
template void GNUPlotter::addDataGrid(int Nx, int Ny, int *x, int *y, int **z);
template void GNUPlotter::addDataGrid(int Nx, int Ny, float *x, float *y, float **z);
template void GNUPlotter::addDataGrid(int Nx, int Ny, double *x, double *y, double **z);

template <class T>
void GNUPlotter::addDataMatrix(int Nx, int Ny, T *x, T *y, T **z)
{
  ofstream out(dataPath, ofstream::app);
  int i, j;
  out << sd(T(Nx));
  for(i = 0; i < Nx; i++)
    out << " " + sd(x[i]);
  out << "\n";
  for(j = 0; j < Ny; j++) {
    out << sd(y[j]);
    for(i = 0; i < Nx; i++)
      out << " " + sd(z[i][j]);
    out << "\n";
  }
  out << "\n\n";
  out.close();
  dataInfo.push_back(DataInfo(1, Nx+1, "matrix"));
}
// todo: allow also flat matrices - i.e. z is not pointer-to-pointer but just a regular pointer
// ...maybe allow the user to select row-major and column-major formats - we should allow a format
// that is compatible with LaPack (flat, column-major)

template <class T>
void GNUPlotter::addDataCurve2D(const std::function<T(T)>& fx, const std::function<T(T)>& fy,
  int Nt, T tMin, T tMax, bool writeT)
{
  vector<T> t(Nt), x(Nt), y(Nt);
  rangeLinear(&t[0], Nt, tMin, tMax);
  for(int i = 0; i < Nt; i++) {
    x[i] = fx(t[i]);
    y[i] = fy(t[i]);
  }
  if(writeT) addDataArrays(Nt, &t[0], &x[0], &y[0]); // write t-values into 1st column
  else       addDataArrays(Nt,        &x[0], &y[0]); // ...or don't

  // the t values may be useful for using color to indicate the t-value along the 
  // curve? or maybe use tick-marks along the curve - for example, if t in in 0..1 ste markers
  // at 0.1, 0.2, 0.9, 1.0 - maybe these markers could be little arrows to also convey the 
  // orientation of the curve? in any case, the markers should not appear on each t - that would
  // be far too dense - maybe every 20th sample or so would be appropriate
  // how about drawing tangent vectors along the curve? make a function addDataCurveWithTangents2D
  // compute the actual tangents numerically
}



template <class T>
void GNUPlotter::addDataSurface(
  const function<T(T, T)>& fx, const function<T(T, T)>& fy, const function<T(T, T)>& fz,
  int Nu, T uMin, T uMax, int Nv, T vMin, T vMax)
{
  // The outer index runs over the indices for parameter u, the middle index runs over v and the 
  // innermost vector index runs from 0...2 giving a 3-vector containing x, y, z coordinates for 
  // each point:
  vector<vector<vector<T>>> d;                   // doubly nested vector of data
  T uStep = (uMax-uMin) / T(Nu-1);               // step size for u
  T vStep = (vMax-vMin) / T(Nv-1);               // step size for v
  d.resize(Nu);                                  // we have Nu blocks of data
  for(int i = 0; i < Nu; i++) {                  // loop over the data blocks
    d[i].resize(Nv);                             // each block has Nv lines/datapoints
    T u = uMin + T(i) * uStep;                   // value of parameter u
    for(int j = 0; j < Nv; j++) {                // loop over lines in current block
      T v = vMin + T(j) * vStep;                 // value of parameter v
      d[i][j].resize(3);                         // each datapoint has 3 columns/dimensions
      d[i][j][0] = fx(u,v);                      // x = fx(u,v)
      d[i][j][1] = fy(u,v);                      // y = fy(u,v)
      d[i][j][2] = fz(u,v);                      // z = fz(u,v)
    }
  }
  addData(d);
}

template <class T>
void GNUPlotter::addDataBivariateFunction(int Nx, int Ny, T *x, T *y, 
  const std::function<T(T, T)>& f)

{
  int i, j;
  T **z = new T*[Nx];
  for(i = 0; i < Nx; i++)
    z[i] = new T[Ny];
  for(i = 0; i < Nx; i++)
  {
    for(j = 0; j < Ny; j++)
      z[i][j] = f(x[i], y[j]);
  }
  addDataMatrix(Nx, Ny, x, y, z);
  for(i = 0; i < Nx; i++)
    delete[] z[i];
  delete[] z;
}

template <class T>
void GNUPlotter::addDataBivariateFunction(int Nx, int Ny, T *x, T *y, T (*f)(T, T))
{
  addDataBivariateFunction(Nx, Ny, x, y, std::function<T(T, T)>(f));
}

template <class T>
void GNUPlotter::addDataBivariateFunction(int Nx, T xMin, T xMax, int Ny, T yMin, T yMax,
  const std::function<T(T, T)>& f)
{
  T *x = new T[Nx];
  T *y = new T[Ny];
  GNUPlotter::rangeLinear(x, Nx, xMin, xMax);
  GNUPlotter::rangeLinear(y, Ny, yMin, yMax);
  addDataBivariateFunction(Nx, Ny, x, y, f);
  delete[] x;
  delete[] y;
}

template <class T>
void GNUPlotter::addDataBivariateFunction(int Nx, T xMin, T xMax, int Ny, T yMin, T yMax,
  T (*f)(T, T))
{
  addDataBivariateFunction(Nx, xMin, xMax, Ny, yMin, yMax, std::function<T(T, T)>(f));
}

template <class T>
void GNUPlotter::addDataVectorField2D(const function<T(T, T)>& fx, const function<T(T, T)>& fy,
  int Nx, T xMin, T xMax, int Ny, T yMin, T yMax)
{
  int Nv = Nx*Ny;                                // number of vectors to draw
  vector<T> x(Nv), y(Nv), dx(Nv), dy(Nv), c(Nv); // arrays to hold our data
  T xStep = (xMax-xMin) / T(Nx-1);               // step size for x
  T yStep = (yMax-yMin) / T(Ny-1);               // step size for y
  T arrowLength = min(xStep, yStep);             // length of arrows to draw
  T s;                                           // length scaler
  for(int i = 0; i < Nx; i++) {                  // loop over x-samples
    for(int j = 0; j < Ny; j++) {                // loop over y-samples
      int k = i*Ny + j;                          // current index in data arrays
      x[k]  = xMin + T(i) * xStep;               // x coordinate (for tail of vector)
      y[k]  = yMin + T(j) * yStep;               // y coordinate (for tail of vector)
      dx[k] = fx(x[k], y[k]);                    // vector's x component
      dy[k] = fy(x[k], y[k]);                    // vector's y component
      c[k]  = (T)hypot(dx[k], dy[k]);            // store length in c for use as color
      if(c[k] != T(0)) s = arrowLength / c[k];   // compute length scaler...
      else             s = T(0);                 // ...catch div-by-0
      dx[k] *= s;                                // adjust the length...
      dy[k] *= s;                                // ...of the vector
    }
  }
  addDataArrays(Nv, &x[0], &y[0], &dx[0], &dy[0], &c[0]);
  // maybe allow the vector field to be evaluated at an arbitrary set of sample points - not 
  // necessarily all on a rectangular grid - we may also want a polar grid...maybe the function 
  // should take a matrix of vectors
  // addDataVectorField2D(fx, fy, Nx, Ny, T** x, T** y, T** dx, T** dy, T** c = nullptr)
  // but maybe allow also flat matrices - i.e. not pointer-to-pointer

  // maybe make not all of them the same length but scale only those which would be too long and
  // leave the shorter ones as is
}

template<class T>
void GNUPlotter::addDataFieldLine2D(const std::function<T(T, T)>& fx, const std::function<T(T, T)>& fy,
  T x0, T y0, T stepSize, int numPoints, int oversampling)
{
  int N = numPoints * oversampling;
  T   h = stepSize  / oversampling;
  std::vector<T> x(N), y(N);
  x[0] = x0;
  y[0] = y0;
  for(int i = 1; i < N; i++) {
    x[i] = x[i-1] + h * fx(x[i-1], y[i-1]);
    y[i] = y[i-1] + h * fy(x[i-1], y[i-1]);
  }
  decimate(&x[0], N, &x[0], oversampling);
  decimate(&y[0], N, &y[0], oversampling);
  addDataArrays(numPoints, &x[0], &y[0]);
}
// maybe factor out the code for the solver - separate it from addData - 
// void solveInitialValueProblem(fx, fy, x0, y0, h, N, bool backward = false) ...the backward 
// option is supposed to solve for negative time - this may be useful for drawing field-lines
// bidirectionally, whose initial point is chosen in the middle of the line/trajectory - the 
// equation becomes x[i] = x[i-1] + sign * h * fx(x[i-1], y[i-1]); and similar for y, sign =+-1

void GNUPlotter::addGraph(CSR descriptor)
{
  graphDescriptors.push_back(descriptor);
}

// addData + addGraph:

template <class T>
void GNUPlotter::addVectorField2D(const function<T(T, T)>& fx, const function<T(T, T)>& fy,
  int Nx, T xMin, T xMax, int Ny, T yMin, T yMax)
{
  addDataVectorField2D(fx, fy, Nx, xMin, xMax, Ny, yMin, yMax);
  addGraph(string("index ") + to_string(dataInfo.size()-1) + 
    string(" using 1:2:3:4:5 with vectors head filled size 0.08,15 ls 2 lc palette notitle"));
}

template<class T>
void GNUPlotter::addFieldLine2D(const std::function<T(T, T)>& fx, const std::function<T(T, T)>& fy,
  T x0, T y0, T stepSize, int numPoints, int oversampling)
{
  addDataFieldLine2D(fx, fy, x0, y0, stepSize, numPoints, oversampling);
  addGraph("index " + to_string(dataInfo.size()-1) + " using 1:2 with lines lt 1 notitle");
}
template void GNUPlotter::addFieldLine2D(const std::function<float(float, float)>& fx, const std::function<float(float, float)>& fy, float x0, float y0, float stepSize, int numPoints, int oversampling);
template void GNUPlotter::addFieldLine2D(const std::function<double(double, double)>& fx, const std::function<double(double, double)>& fy, double x0, double y0, double stepSize, int numPoints, int oversampling);

// inquiry:

string GNUPlotter::getDataPath()
{
  return dataPath;
}

// misc:

template <class T>
void GNUPlotter::rangeLinear(T *x, int N, T min, T max)
{
  double factor = (max-min) / (double)(N-1);
  for(int i = 0; i < N; i++)
    x[i] = (T)(factor * T(i) + min);
}

template <class T>
void GNUPlotter::rangeLogarithmic(T *x, int N, T min, T max)
{
  double a = log(max/min) / T(N-1);
  for(int i = 0; i < N; i++)
    x[i] = (T) (min*exp(a*i));
}
template void GNUPlotter::rangeLogarithmic(float *x, int N, float min, float max);
template void GNUPlotter::rangeLogarithmic(double *x, int N, double min, double max);

template <class T>
void GNUPlotter::decimate(T* x, int Nx, T* y, int factor)
{
  int Ny = Nx / factor;
  for(int i = 0; i < Ny; i++)
    y[i] = x[i*factor];
}

void GNUPlotter::clearCommandFile()
{
  initFile(commandPath);
}

void GNUPlotter::invokeGNUPlot()
{
  // create the callstring and invoke GNUPlot:
  //string callString = "\"" + gnuplotPath + "\" " + commandPath + " -";
  //string callString = "\"" + gnuplotPath + "\" " + commandPath;
  string callString = "\"" + gnuplotPath + "\" " + commandPath + " -persist";

  // wrapping gnuplotPath into quotes is required to handle installation paths with whitespaces,
  // but it doesn't seem to be possible to handle a commandPath with whitespaces (i tried
  // wrapping it into quotes too and wrapping the whole callString into quotes - none of that
  // works)
  // the minus at the end prevents gnuplot from immediately closing

  systemCall(callString);
  // is it possible to call GNUPlot in a separate process? - this is actually already the case
  // ...but we have to close it to call it again - i.e. we can't do multiple plots at once
  int dummy = 0;
}

// internal functions:

void GNUPlotter::initFile(const std::string &path)
{
  ofstream out(path);
  out.close();
}

char* GNUPlotter::toZeroTerminatedString(const std::string &s)
{
  size_t N = s.size();
  char *zs = new char[N+1];
  for(size_t i = 0; i < N; i++)
    zs[i] = s[i];
  zs[N] = '\0';
  return zs;
}

void GNUPlotter::systemCall(const std::string &callString)
{
  char *cString = toZeroTerminatedString(callString);
  system(cString);
  delete[] cString;
}

std::string GNUPlotter::getStyleString(unsigned int i)
{
  std::string s;
  //if(dataInfo[i].type == "matrix")
  //  s += " nonuniform matrix";
  s += " w ";
  s += graphStyles[i % graphStyles.size()] + " ";
  s += getGraphLegend(i);
  return s;
}

std::string GNUPlotter::getGraphLegend(unsigned int i)
{
  if( i >= graphTitles.size() )
    return "notitle";
  else
  {
    if( graphTitles[i].empty() )
      return "notitle";
    return "t \"" + graphTitles[i] + "\"";
  }
}

std::string GNUPlotter::s(unsigned int x)
{
  //return std::to_string((_Longlong)x); // does not compile with gcc on linux
  return std::to_string((long long)x);
}

std::string GNUPlotter::s(double x)
{
  char cString[32];

#ifdef _MSC_VER
  sprintf_s(cString, 32, "%-.16g", x);
#else
  snprintf(cString, 32, "%-.16g", x);
#endif

  return std::string(cString);
}

std::string GNUPlotter::sd(double x)
{
  char cString[32];

#ifdef _MSC_VER
  sprintf_s(cString, 32, formatString, x);
#else
  snprintf(cString, 32, formatString, x);
#endif

  return std::string(cString);
}

std::string GNUPlotter::sd(int x)
{
  return to_string(x);
}

// internal functions

void GNUPlotter::addPlotCommand(bool splot)
{
  int i;
  string pc;
  addCommand("\n# Plotting:");
  if( splot == true )
    pc = "splot \\\n";  // 3D plots
  else
    pc = "plot \\\n";   // 2D plots
  if(graphDescriptors.empty())
    generateGraphDescriptors(splot);
  for(i = 0; i < (int)graphDescriptors.size()-1; i++)
    pc += "'" + dataPath + "' " + graphDescriptors[i] + ",\\\n";
  pc += "'" + dataPath + "' " + graphDescriptors[i];
  addCommand(pc);
}

void GNUPlotter::generateGraphDescriptors(bool splot)
{
  unsigned int i, j, k;
  if(splot == false) {
    k = 0;
    for(i = 0; i < dataInfo.size(); i++)
    {
      if(dataInfo[i].type == "matrix")
      {
        addGraph("i " + s(i) + " nonuniform matrix" + getStyleString(k));
        k++;
      }
      else
      {
        if(dataInfo[i].numColumns == 1)  // this is new -> test...
        {
          //addGraph("i " + s(i) + " u 1:1" + getStyleString(k));
          addGraph("i " + s(i) + getStyleString(k));
          k++;
        }
        // entered only, if condition above is false:
        for(j = 0; j < (unsigned int)dataInfo[i].numColumns-1; j++)
        {
          addGraph("i " + s(i) + " u 1:" + s(j+2) + getStyleString(k));
          k++;
        }
      }
    }
  }
  else
  {
    for(i = 0; i < dataInfo.size(); i++)
    {
      if(dataInfo[i].type == "matrix")
        addGraph("i " + s(i) + " nonuniform matrix" + getStyleString(i));
      else
        addGraph("i " + s(i) + getStyleString(i));
    }
  }
}

// // old:
//void GNUPlotter::generateGraphDescriptors(bool splot)
//{
//  unsigned int i, j, k;
//  if(splot == false) {
//    k = 0;
//    for(i = 0; i < dataInfo.size(); i++)
//    {
//      if(dataInfo[i].type == "matrix")
//      {
//        addGraph("i " + s(i) + " nonuniform matrix" + getStyleString(k));
//        k++;
//      }
//      else
//      {
//        for(j = 0; j < (unsigned int)dataInfo[i].numColumns-1; j++)
//        {
//          addGraph("i " + s(i) + " u 1:" + s(j+2) + getStyleString(k));
//          k++;
//        }
//      }
//    }
//  }
//  else
//  {
//    for(i = 0; i < dataInfo.size(); i++)
//    {
//      if(dataInfo[i].type == "matrix")
//        addGraph("i " + s(i) + " nonuniform matrix" + getStyleString(i));
//      else
//        addGraph("i " + s(i) + getStyleString(i));
//    }
//  }
//}




// argument handling

template<class T>
T GNUPlotter::nullValue(T)
{
  return T(0);
}

string GNUPlotter::nullValue(string)
{
  return "";
}

template<class T>
const vector<T> GNUPlotter::collectLeadingNonNullArguments(const T a0, const T a1, const T a2, 
  const T a3, const T a4, const T a5, const T a6, const T a7, const T a8, const T a9)
{
  T null = nullValue(a0);
  vector<T> v;
  if(a0 != null) v.push_back(a0); else return v;
  if(a1 != null) v.push_back(a1); else return v;
  if(a2 != null) v.push_back(a2); else return v;
  if(a3 != null) v.push_back(a3); else return v;
  if(a4 != null) v.push_back(a4); else return v;
  if(a5 != null) v.push_back(a5); else return v;
  if(a6 != null) v.push_back(a6); else return v;
  if(a7 != null) v.push_back(a7); else return v;
  if(a8 != null) v.push_back(a8); else return v;
  if(a9 != null) v.push_back(a9); else return v;
  return v;
}

template<class T>
void GNUPlotter::append(vector<T>& v, const vector<T>& appendix)
{
  v.reserve(v.size() + appendix.size());
  for(size_t i = 0; i < appendix.size(); i++)
    v.push_back(appendix[i]);
}

template<class T>
vector<vector<T>> GNUPlotter::wrapIntoVectors(int N, const T *a0, const T *a1, const T *a2, 
  const T *a3, const T *a4, const T *a5, const T *a6, const T *a7, const T *a8, const T *a9)
{
  vector<T*> pointers = collectLeadingNonNullArguments(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
  vector<vector<T>> v;
  v.resize(pointers.size());
  for(int i = 0; i < v.size(); i++)
  {
    v[i].resize(N);
    for(int j = 0; j < N; j++)
      v[i][j] = pointers[i][j];
  }
  return v;
}

void GNUPlotter::addToStringVector(vector<string>& v, CSR s0, CSR s1, CSR s2, CSR s3, CSR s4,
  CSR s5, CSR s6, CSR s7, CSR s8, CSR s9)
{
  append(v, collectLeadingNonNullArguments(s0, s1, s2, s3, s4, s5, s6, s7, s8, s9));
}

void GNUPlotter::setStringVector(vector<string>& v, CSR s0, CSR s1, CSR s2, CSR s3, CSR s4, CSR s5,
  CSR s6, CSR s7, CSR s8, CSR s9)
{
  v.clear();
  addToStringVector(v, s0, s1, s2, s3, s4, s5, s6, s7, s8, s9);
}

/*
ToDo:
-maybe move the explicit template instantiations to another file...that would reduce clutter in 
 this implementation file - but would make the library harder to use - the user would have to deal
 with more files...so it's probably not such a good idea...simple use is more important than
 nice looking code

*/