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
  gnuplotPath = "C:/Program Files/gnuplot/bin/gnuplot.exe";
  //gnuplotPath = "C:/Program Files/gnuplot/bin/wgnuplot.exe";

  // paths for data file and command batchfile:
  dataPath    = "E:/Temp/gnuplotData.dat";
  commandPath = "E:/Temp/gnuplotCommands.txt";   // this path may not contain whitepaces

  // initialize data- and commandfile:
  initFile(dataPath);
  initFile(commandPath);
  addDefaultCommands();
}

GNUPlotter::~GNUPlotter()
{
  delete[] formatString;
}

// plotting:

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
void GNUPlotter::plotFunctionTables(int N, T *x, T *y1, T *y2, T *y3, T *y4, T *y5, T *y6, T *y7,
  T *y8, T *y9)
{
  addDataArrays(N, x, y1, y2, y3, y4, y5, y6, y7, y8, y9);
  plot();
}

template <class T>
void GNUPlotter::plotArrays(int N, T *y1, T *y2, T *y3, T *y4, T *y5, T *y6, T *y7, T *y8, 
  T *y9)
{
  T *x = new T[N];
  rangeLinear(x, N, T(0), T(N-1));
  plotFunctionTables(N, x, y1, y2, y3, y4, y5, y6, y7, y8, y9);
  delete[] x;
}
// explicit instantiations for double, float and int:
template void GNUPlotter::plotArrays(int N, double *y1, double *y2, double *y3, double *y4, 
  double *y5, double *y6, double *y7, double *y8, double *y9);
template void GNUPlotter::plotArrays(int N, float *y1, float *y2, float *y3, float *y4, float *y5, 
  float *y6, float *y7, float *y8, float *y9);
template void GNUPlotter::plotArrays(int N, int *y1, int *y2, int *y3, int *y4, int *y5, 
  int *y6, int *y7, int *y8, int *y9);

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
void GNUPlotter::addDataArrays(int N, T *c0, T *c1, T *c2, T *c3, T *c4, T *c5, T *c6, T *c7,
  T *c8, T *c9)
{
  vector<T*> v = collectLeadingNonNullArguments(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
  T* a[10];
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
  for(j = 0; j < Ny; j++)
  {
    out << sd(y[j]);
    for(i = 0; i < Nx; i++)
      out << " " + sd(z[i][j]);
    out << "\n";
  }
  out << "\n\n";
  out.close();
  dataInfo.push_back(DataInfo(1, Nx+1, "matrix"));
}

template <class T>
void GNUPlotter::addDataBivariateFunction(int Nx, int Ny, T *x, T *y, T (*f)(T, T))
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
void GNUPlotter::addDataBivariateFunction(int Nx, T xMin, T xMax, int Ny, T yMin, T yMax, 
  T (*f)(T, T))
{
  T *x = new T[Nx];
  T *y = new T[Ny];
  GNUPlotter::rangeLinear(x, Nx, xMin, xMax);
  GNUPlotter::rangeLinear(y, Ny, yMin, yMax);
  addDataBivariateFunction(Nx, Ny, x, y, f);
  delete[] x;
  delete[] y;
}

void GNUPlotter::addGraph(CSR descriptor)
{
  graphDescriptors.push_back(descriptor);
}

// inquiry:

string GNUPlotter::getDataPath()
{
  return dataPath;
}

// misc:

template <class T>
static void GNUPlotter::rangeLinear(T *x, int N, T min, T max)
{
  double factor = (max-min) / (double)(N-1);
  for(int i = 0; i < N; i++)
    x[i] = (T)(factor * T(i) + min);
}

template <class T>
static void GNUPlotter::rangeLogarithmic(T *x, int N, T min, T max)
{
  double a = log(max/min) / T(N-1);
  for(int i = 0; i < N; i++)
    x[i] = (T) (min*exp(a*i));
}
template void GNUPlotter::rangeLogarithmic(float *x, int N, float min, float max);
template void GNUPlotter::rangeLogarithmic(double *x, int N, double min, double max);

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
  return std::to_string((_Longlong)x);
}

std::string GNUPlotter::s(double x)
{
  char cString[32];
  sprintf_s(cString, 32, "%-.16g", x);
  return std::string(cString);
}

std::string GNUPlotter::sd(double x)
{
  char cString[32];
  sprintf_s(cString, 32, formatString, x);
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
vector<T> GNUPlotter::collectLeadingNonNullArguments(T a0, T a1, T a2, T a3, T a4, T a5, T a6,
  T a7, T a8, T a9)
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
vector<vector<T>> GNUPlotter::wrapIntoVectors(int N, T *a0, T *a1, T *a2, T *a3, T *a4, T *a5, 
  T *a6, T *a7, T *a8, T *a9)
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
