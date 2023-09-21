#include "GNUPlotter.h"
#include <fstream>
#include <iostream>
#include <algorithm>
//using namespace std; // bad idea when including it in rapt.cpp

GNUPlotter::GNUPlotter()
{
  setDataPrecision(8);

  graphStyles.resize(1);
  graphStyles[0] = std::string("lines");    // standard lines, width 1 - maybe change to 1.5 or 2

  // installation path of the GNUPlot executable:
  gnuplotPath = "C:/Program Files/gnuplot/bin/gnuplot.exe";
  //gnuplotPath = "C:/Program Files/gnuplot/bin/wgnuplot.exe";

  // paths for data file and command batchfile:
  dataPath    = "C:/Temp/gnuplotData.dat";
  commandPath = "C:/Temp/gnuplotCommands.txt";   // this path may not contain whitepaces

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
void GNUPlotter::plotVectorField2D(
  const std::function<T(T, T)>& fx, const std::function<T(T, T)>& fy,
  int Nx, T xMin, T xMax, int Ny, T yMin, T yMax)
{
  GNUPlotter p;
  p.addVectorField2D(fx, fy, Nx, xMin, xMax, Ny, yMin, yMax);

  p.addCommand("set palette rgbformulae 30,31,32 negative");  
  // use setColorPalette or get rid entirely!

  p.plot();
  // maybe try a black background (and invert the colormap)
}
template void GNUPlotter::plotVectorField2D(const std::function<int(int, int)>& fx, const std::function<int(int, int)>& fy, int Nx, int xMin, int xMax, int Ny, int yMin, int yMax);
template void GNUPlotter::plotVectorField2D(const std::function<float(float, float)>& fx, const std::function<float(float, float)>& fy, int Nx, float xMin, float xMax, int Ny, float yMin, float yMax);
template void GNUPlotter::plotVectorField2D(const std::function<double(double, double)>& fx, const std::function<double(double, double)>& fy, int Nx, double xMin, double xMax, int Ny, double yMin, double yMax);


template<class T>
void GNUPlotter::plotComplexVectorField(const std::function<std::complex<T>(std::complex<T>)>& f,
  int Nr, T rMin, T rMax, int Ni, T iMin, T iMax, bool conj)
{
  T sign = T(1); if(conj) sign = T(-1);
  std::function<T(T, T)> fx, fy;
  fx = [&] (T re, T im) { return        real(f(std::complex<T>(re, im))); };
  fy = [&] (T re, T im) { return sign * imag(f(std::complex<T>(re, im))); };
  plotVectorField2D(fx, fy, Nr, rMin, rMax, Ni, iMin, iMax);
}
template void GNUPlotter::plotComplexVectorField(const std::function<std::complex<int>(std::complex<int>)>& f, int Nr, int rMin, int rMax, int Ni, int iMin, int iMax, bool conj);
template void GNUPlotter::plotComplexVectorField(const std::function<std::complex<float>(std::complex<float>)>& f, int Nr, float rMin, float rMax, int Ni, float iMin, float iMax, bool conj);
template void GNUPlotter::plotComplexVectorField(const std::function<std::complex<double>(std::complex<double>)>& f, int Nr, double rMin, double rMax, int Ni, double iMin, double iMax, bool conj);


//-------------------------------------------------------------------------------------------------

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
}

template <class T>
void GNUPlotter::plotBivariateFunction(int Nx, T xMin, T xMax, int Ny, T yMin, T yMax,
  T (*f)(T, T))
{
  addDataBivariateFunction(Nx, xMin, xMax, Ny, yMin, yMax, f);
  plot3D();
}
template void GNUPlotter::plotBivariateFunction(int Nx, double xMin, double xMax, int Ny,
  double yMin, double yMax, double (*f)(double, double));
template void GNUPlotter::plotBivariateFunction(int Nx, float xMin, float xMax, int Ny, float yMin,
  float yMax, float (*f)(float, float));
template void GNUPlotter::plotBivariateFunction(int Nx, int xMin, int xMax, int Ny, int yMin,
  int yMax, int (*f)(int, int));

template <class T>
void GNUPlotter::plotBivariateFunction(int Nx, T xMin, T xMax, int Ny, T yMin, T yMax,
  const std::function<T(T, T)>& f)
{
  addDataBivariateFunction(Nx, xMin, xMax, Ny, yMin, yMax, f);
  plot3D();
}
template void GNUPlotter::plotBivariateFunction(int Nx, double xMin, double xMax, int Ny,
  double yMin, double yMax, const std::function<double(double, double)>& f);
template void GNUPlotter::plotBivariateFunction(int Nx, float xMin, float xMax, int Ny,
  float yMin, float yMax, const std::function<float(float, float)>& f);
template void GNUPlotter::plotBivariateFunction(int Nx, int xMin, int xMax, int Ny,
  int yMin, int yMax, const std::function<int(int, int)>& f);

template<class T>
void GNUPlotter::plotContourMap(int Nx, T xMin, T xMax, int Ny, T yMin, T yMax,
  const std::function<T(T, T)>& f, int numContours, T zMin, T zMax)
{
  addDataBivariateFunction(Nx, xMin, xMax, Ny, yMin, yMax, f);
  std::vector<T> levels(numContours); 
  rangeLinear(&levels[0], numContours, zMin, zMax);
  setContourLevels(levels);

  // Use constant color fills between the contour lines if desired:
  bool useConstColors = true;  // make user parameter
  if(useConstColors)
  {
    std::string cmd = "set palette maxcolors " + std::to_string(levels.size() - 1);
    addCommand(cmd);
    size_t L = levels.size() - 1;        // last valid index
    std::string range = "[" + std::to_string(levels[0]) + ":" + std::to_string(levels[L]) + "]";
    addCommand("set zrange " + range);   // range for z values
    addCommand("set cbrange " + range);  // color bar range
  }

  // Plot:
  addCommand("set pm3d map impl");
  addCommand("set contour");
  addCommand("splot '" + dataPath + "' i 0 nonuniform matrix w pm3d notitle");
  //addCommand("set autoscale fix");
  invokeGNUPlot();

  // ToDo:
  // -When clicking on "Apply autoscale" on the GUI, the colors get messed up. I'm trying to fix
  //  this problem via "set autoscale fix" but that doesn't seem to help.
}
template void GNUPlotter::plotContourMap(
  int Nx, double xMin, double xMax, int Ny, double yMin, double yMax,
  const std::function<double(double, double)>& f, int numContours, double zMin, double zMax);
template void GNUPlotter::plotContourMap(
  int Nx, float xMin, float xMax, int Ny, float yMin, float yMax,
  const std::function<float(float, float)>& f, int numContours, float zMin, float zMax);
template void GNUPlotter::plotContourMap(
  int Nx, int xMin, int xMax, int Ny, int yMin, int yMax,
  const std::function<int(int, int)>& f, int numContours, int zMin, int zMax);

//-------------------------------------------------------------------------------------------------
// style setup:

void GNUPlotter::addCommand(std::string command)
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

void GNUPlotter::setToDarkMode()
{
  addCommand("set term wxt background rgb \"black\"");
  addCommand("set border lw 1 lc rgb \"white\"");
  addCommand("set grid ls 1 lw 1 lc rgb \"#404040\""); // old: "set grid lw 1 lc rgb \"white\""
  addCommand("set xtics textcolor rgb \"white\"");
  addCommand("set ytics textcolor rgb \"white\"");
  addCommand("set xlabel \"X\" textcolor rgb \"white\"");
  addCommand("set ylabel \"Y\" textcolor rgb \"white\"");
  addCommand("set key textcolor \"white\""); 
  addCommand("set style line 1 lt 1 lw 3 pt 3 linecolor rgb \"#F0F0F0\""); // What does this do?

  // Drawing All line white may be not so good:
  const char c[7] = "FFFFFF";
  setGraphColors(c, c, c, c, c, c, c, c, c, c);
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

void GNUPlotter::setLegends(CVR<std::string> legends)
{
  graphTitles = legends;
}

void GNUPlotter::setColorPalette(ColorPalette palette, bool inverted)
{
  std::string c;            // The command to be passed to plt
  using CP = ColorPalette;
  switch(palette)
  {
    // Linear:
  case CP::CB_Blues8:       c = "set palette defined (0 '#F7FBFF', 1 '#DEEBF7', 2 '#C6DBEF', 3 '#9ECAE1', 4 '#6BAED6', 5 '#4292C6', 6 '#2171B5', 7 '#084594')"; break;
  //case CP::CB_BuGn8:        c = "set palette defined (0 '#F7FCFD', 1 '#E5F5F9', 2 '#CCECE6', 3 '#99D8C9', 4 '#66C2A4', 5 '#41AE76', 6 '#238B45', 7 '#005824')"; break;
  //case CP::CB_BuGn8m:       c = "set palette defined (0 '#F7FCFD', 1 '#CCECE6', 2 '#99D8C9', 3 '#66C2A4', 4 '#41AE76', 5 '#238B45', 6 '#005824')"; break;
  case CP::CB_BuPu8:        c = "set palette defined (0 '#F7FCFD', 1 '#E0ECF4', 2 '#BFD3E6', 3 '#9EBCDA', 4 '#8C96C6', 5 '#8C6BB1', 6 '#88419D', 7 '#6E016B')"; break;
  case CP::CB_GnBu8:        c = "set palette defined (0 '#F7FCF0', 1 '#E0F3DB', 2 '#CCEBC5', 3 '#A8DDB5', 4 '#7BCCC4', 5 '#4EB3D3', 6 '#2B8CBE', 7 '#08589E')"; break;
  //case CP::CB_Greens8:      c = "set palette defined (0 '#F7FCF5', 1 '#E5F5E0', 2 '#C7E9C0', 3 '#A1D99B', 4 '#74C476', 5 '#41AB5D', 6 '#238B45', 7 '#005A32')"; break;
  //case CP::CB_Oranges8:     c = "set palette defined (0 '#FFF5EB', 1 '#FEE6CE', 2 '#FDD0A2', 3 '#FDAE6B', 4 '#FD8D3C', 5 '#F16913', 6 '#D94801', 7 '#8C2D04')"; break;
  case CP::CB_PuBu8:        c = "set palette defined (0 '#FFF7FB', 1 '#ECE7F2', 2 '#D0D1E6', 3 '#A6BDDB', 4 '#74A9CF', 5 '#3690C0', 6 '#0570B0', 7 '#034E7B')"; break;
  //case CP::CB_Purples8:     c = "set palette defined (0 '#FCFBFD', 1 '#EFEDF5', 2 '#DADAEB', 3 '#BCBDDC', 4 '#9E9AC8', 5 '#807DBA', 6 '#6A51A3', 7 '#4A1486')"; break;
  case CP::CB_RdPu8:        c = "set palette defined (0 '#FFF7F3', 1 '#FDE0DD', 2 '#FCC5C0', 3 '#FA9FB5', 4 '#F768A1', 5 '#DD3497', 6 '#AE017E', 7 '#7A0177')"; break;
  //case CP::CB_Reds8:        c = "set palette defined (0 '#FFF5F0', 1 '#FEE0D2', 2 '#FCBBA1', 3 '#FC9272', 4 '#FB6A4A', 5 '#EF3B2C', 6 '#CB181D', 7 '#99000D')"; break;
  case CP::CB_YlGn8:        c = "set palette defined (0 '#FFFFE5', 1 '#F7FCB9', 2 '#D9F0A3', 3 '#ADDD8E', 4 '#78C679', 5 '#41AB5D', 6 '#238443', 7 '#005A32')"; break;
  case CP::CB_YlGnBu8:      c = "set palette defined (0 '#FFFFD9', 1 '#EDF8B1', 2 '#C7E9B4', 3 '#7FCDBB', 4 '#41B6C4', 5 '#1D91C0', 6 '#225EA8', 7 '#0C2C84')"; break;
  case CP::CB_YlOrBr8:      c = "set palette defined (0 '#FFFFE5', 1 '#FFF7BC', 2 '#FEE391', 3 '#FEC44F', 4 '#FE9929', 5 '#EC7014', 6 '#CC4C02', 7 '#8C2D04')"; break;
  case CP::CB_YlOrRd8:      c = "set palette defined (0 '#FFFFCC', 1 '#FFEDA0', 2 '#FED976', 3 '#FEB24C', 4 '#FD8D3C', 5 '#FC4E2A', 6 '#E31A1C', 7 '#B10026')"; break;
  case CP::CB_YlGnBu9:      c = "set palette defined (0 '#ffffd9', 1 '#edf8b1', 2 '#c7e9b4', 3 '#7fcdbb', 4 '#41b6c4', 5 '#1d91c0', 6 '#225ea8', 7 '#253494', 8 '#081d58')"; break;
  case CP::CB_YlGnBu9m:     c = "set palette defined (0 '#ffffd9', 1 '#c7e9b4', 2 '#7fcdbb', 3 '#41b6c4', 4 '#1d91c0', 5 '#225ea8', 6 '#253494', 7 '#081d58')"; break;

  case CP::CB_YlOrBr9:      c = "set palette defined (0 '#ffffe5', 1 '#fff7bc', 2 '#fee391', 3 '#fec44f', 4 '#fe9929', 5 '#ec7014', 6 '#cc4c02', 7 '#993404', 8 '#662506')"; break;
  case CP::CB_YlOrRd9:      c = "set palette defined (0 '#ffffcc', 1 '#ffeda0', 2 '#fed976', 3 '#feb24c', 4 '#fd8d3c', 5 '#fc4e2a', 6 '#e31a1c', 7 '#bd0026', 8 '#800026')"; break;
  case CP::EF_Viridis:      c = "set palette defined (0 '#440154', 1 '#472c7a', 2 '#3b518b', 3 '#2c718e', 4 '#21908d', 5 '#27ad81', 6 '#5cc863', 7 '#aadc32', 8 '#fde725')"; break;
  case CP::GF_AfmHot:       c = "set palette rgbformulae 34,35,36"; break;
  case CP::GF_BkPuWt:       c = "set palette rgbformulae 3,23,21"; break;
  case CP::GF_Hot:          c = "set palette rgbformulae 21,22,23"; break;
  case CP::GF_Printable:    c = "set palette rgbformulae 30,31,32"; break;
  case CP::GF_TradPm3d:     c = "set palette rgbformulae  7, 5,15"; break;
  case CP::GP_Sand:         c = "set palette defined (0 '#604860', 1 '#784860', 2 '#a86060', 3 '#c07860', 4 '#f0a848', 5 '#f8ca8c', 6 '#feecae', 7 '#fff4c2', 8 '#fff7db', 9 '#fffcf6')"; break;
  case CP::ML_Parula:       c = "set palette defined (0 '#352a87', 1 '#0363e1', 2 '#1485d4', 3 '#06a7c6', 4 '#38b99e', 5 '#92bf73', 6 '#d9ba56', 7 '#fcce2e', 8 '#f9fb0e')"; break;
  case CP::RS_BkWt:         c = "set palette defined (0 '#000000', 1 '#ffffff')"; break;
  case CP::SW_Inferno:      c = "set palette defined (0 '#000004', 1 '#1f0c48', 2 '#550f6d', 3 '#88226a', 4 '#a83655', 5 '#e35933', 6 '#f9950a', 7 '#f8c932', 8 '#fcffa4')"; break;  
  case CP::SW_Magma:        c = "set palette defined (0 '#000004', 1 '#1c1044', 2 '#4f127b', 3 '#812581', 4 '#b5367a', 5 '#e55964', 6 '#fb8761', 7 '#fec287', 8 '#fbfdbf')"; break;
  case CP::SW_Plasma:       c = "set palette defined (0 '#0c0887', 1 '#4b03a1', 2 '#7d03a8', 3 '#a82296', 4 '#cb4679', 5 '#e56b5d', 6 '#f89441', 7 '#fdc328', 8 '#f0f921')"; break;
  //case CP::UA_YlRd:         c = "set palette defined (0 '#ffee00', 1 '#ff7000', 2 '#ee0000', 3 '#7f0000')"; break;
  case CP::CJ_YlRd9:        c = "set palette defined (0 '#ffffe0', 1 '#ffdfb8', 2 '#ffbc94', 3 '#ff9777', 4 '#ff6962', 5 '#ee4256', 6 '#d21f47', 7 '#b0062c', 8 '#8b0000')"; break;  

    // Diverging:
  case CP::AM_Turbo:        c = "set palette defined (0 '#30123b', 1 '#466be3', 2 '#28bceb', 3 '#32f298', 4 '#a4fc3c', 5 '#eecf3a', 6 '#fb7e21', 7 '#d02f05', 8 '#7a0403')"; break;
  
  case CP::CB_BrBG8:        c = "set palette defined (0 '#8C510A', 1 '#BF812D', 2 '#DFC27D', 3 '#F6E8C3', 4 '#C7EAE5', 5 '#80CDC1', 6 '#35978F', 7 '#01665E')"; break;
  case CP::CB_BrBG9:        c = "set palette defined (0 '#8c510a', 1 '#bf812d', 2 '#dfc27d', 3 '#f6e8c3', 4 '#f5f5f5', 5 '#c7eae5', 6 '#80cdc1', 7 '#35978f', 8 '#01665e')"; break;
  
  case CP::CB_PiYG8:        c = "set palette defined (0 '#C51B7D', 1 '#DE77AE', 2 '#F1B6DA', 3 '#FDE0EF', 4 '#E6F5D0', 5 '#B8E186', 6 '#7FBC41', 7 '#4D9221')"; break;
  case CP::CB_PRGn8:        c = "set palette defined (0 '#762A83', 1 '#9970AB', 2 '#C2A5CF', 3 '#E7D4E8', 4 '#D9F0D3', 5 '#A6DBA0', 6 '#5AAE61', 7 '#1B7837')"; break;
  case CP::CB_PuOr8:        c = "set palette defined (0 '#B35806', 1 '#E08214', 2 '#FDB863', 3 '#FEE0B6', 4 '#D8DAEB', 5 '#B2ABD2', 6 '#8073AC', 7 '#542788')"; break;
  case CP::CB_RdBu8:        c = "set palette defined (0 '#B2182B', 1 '#D6604D', 2 '#F4A582', 3 '#FDDBC7', 4 '#D1E5F0', 5 '#92C5DE', 6 '#4393C3', 7 '#2166AC')"; break;
  case CP::CB_RdYlBu8:      c = "set palette defined (0 '#D73027', 1 '#F46D43', 2 '#FDAE61', 3 '#FEE090', 4 '#E0F3F8', 5 '#ABD9E9', 6 '#74ADD1', 7 '#4575B4')"; break;
  case CP::CB_RdYlGn8:      c = "set palette defined (0 '#D73027', 1 '#F46D43', 2 '#FDAE61', 3 '#FEE08B', 4 '#D9EF8B', 5 '#A6D96A', 6 '#66BD63', 7 '#1A9850')"; break;
  case CP::CB_Spectral8:    c = "set palette defined (0 '#D53E4F', 1 '#F46D43', 2 '#FDAE61', 3 '#FEE08B', 4 '#E6F598', 5 '#ABDDA4', 6 '#66C2A5', 7 '#3288BD')"; break;

  case CP::CB_RdBu11:       c = "set palette defined (0 '#67001f', 1 '#b2182b', 2 '#d6604d', 3 '#f4a582', 4 '#fddbc7', 5 '#f7f7f7', 6 '#d1e5f0', 7 '#92c5de', 8 '#4393c3', 9 '#2166ac', 10 '#053061')"; break;  
  case CP::CB_RdYlBu11:     c = "set palette defined (0 '#a50026', 1 '#d73027', 2 '#f46d43', 3 '#fdae61', 4 '#fee090', 5 '#ffffbf', 6 '#e0f3f8', 7 '#abd9e9', 8 '#74add1', 9 '#4575b4', 10 '#313695')"; break;
  case CP::CB_RdYlGn11:     c = "set palette defined (0 '#a50026', 1 '#d73027', 2 '#f46d43', 3 '#fdae61', 4 '#fee08b', 5 '#ffffbf', 6 '#d9ef8b', 7 '#a6d96a', 8 '#66bd63', 9 '#1a9850', 10 '#006837')"; break;
  case CP::CB_PRGn11:       c = "set palette defined (0 '#40004b', 1 '#762a83', 2 '#9970ab', 3 '#c2a5cf', 4 '#e7d4e8', 5 '#f7f7f7', 6 '#d9f0d3', 7 '#a6dba0', 8 '#5aae61', 9 '#1b7837', 10 '#00441b')"; break;
  case CP::CB_Spectral11:   c = "set palette defined (0 '#9e0142', 1 '#d53e4f', 2 '#f46d43', 3 '#fdae61', 4 '#fee08b', 5 '#ffffbf', 6 '#e6f598', 7 '#abdda4', 8 '#66c2a5', 9 '#3288bd', 10 '#5e4fa2')"; break;
  case CP::CJ_BuYlRd11:     c = "set palette defined (0 '#00429d', 1 '#3d68aa', 2 '#6190b7', 3 '#86b8c4', 4 '#b6ded1', 5 '#ffffe0', 6 '#ffcab9', 7 '#fd9291', 8 '#e75d6f', 9 '#c52a52', 10 '#93003a')"; break;

  case CP::GF_PuGnRd:       c = "set palette rgbformulae 33,13,10"; break;
  case CP::KM_Moreland:     c = "set palette defined (0 '#3b4cc0', 1 '#688aef', 2 '#99baff', 3 '#c9d8ef', 4 '#edd1c2', 5 '#f7a789', 6 '#e36a53', 7 '#b40426')"; break;
  case CP::KM_BentCoolWarm: c = "set palette defined (0 '#5548c1', 1 '#7982d7', 2 '#abb8e7', 3 '#dde3ef', 4 '#ead3c6', 5 '#dba188', 6 '#ca6b55', 7 '#b10027')"; break;
  case CP::ML_Jet:          c = "set palette defined (0 '#000080', 1 '#0000ff', 2 '#0080ff', 3 '#00ffff', 4 '#80ff80', 5 '#ffff00', 6 '#ff8000', 7 '#ff0000', 8 '#800000')"; break;
  case CP::RS_RdGnBu:       c = "set palette defined (0 '#200000', 1 '#e0fff0', 2 '#000030')"; break;
  case CP::UA_GnPu:         c = "set palette defined (1 '#396353', 2 '#0db14b', 3 '#6dc067', 4 '#abd69b', 5 '#daeac1', 6 '#dfcce4', 7 '#c7b2d6', 8 '#9474b4', 9 '#754098', 10 '#504971')"; break;
  
    // Alternating:
  case CP::CB_Paired8:      c = "set palette defined (0 '#A6CEE3', 1 '#1F78B4', 2 '#B2DF8A', 3 '#33A02C', 4 '#FB9A99', 5 '#E31A1C', 6 '#FDBF6F', 7 '#FF7F00')"; break;
  case CP::CB_Paired10:   c = "set palette defined (0 '#a6cee3', 1 '#1f78b4', 2 '#b2df8a', 3 '#33a02c', 4 '#fb9a99', 5 '#e31a1c', 6 '#fdbf6f', 7 '#ff7f00', 8 '#cab2d6', 9 '#6a3d9a')"; break;
  // The ones with smaller count can actually be obtained by truncating the longer ones ...well...
  // not always, I think.
  
    // Cyclic:

    // Misc:
  //case CP::_test:   c = "set palette rgbformulae 8,9,7"; break;

  // ...more to come
  }
  
  // ToDo:
  // remove the common "set palette " from the cases and add it to the string like this:
  // c = "set palette " + c;

  if(inverted)
    c += " negative";

  addCommand(c);
}

void GNUPlotter::setGraphColors(CSR c0, CSR c1, CSR c2, CSR c3, CSR c4, CSR c5, CSR c6, CSR c7,
  CSR c8, CSR c9)
{
  std::vector<std::string> v;
  setStringVector(v, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
  setGraphColors(v);
}

void GNUPlotter::setGraphColors(CVR<std::string> c)
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

void GNUPlotter::setGraphStyles(CVR<std::string> styles)
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

void GNUPlotter::setLogScale(std::string axes, double /*base*/, bool shouldBeLogarithmic)
{
  std::string s;
  if( !shouldBeLogarithmic )
    s += "un";
  s += "set logscale " + axes + "\n";
  addCommand(s);
}
// todo: use base or get rid of the parameter

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

template<class T>
void GNUPlotter::setContourLevels(const std::vector<T>& levels)
{
  std::string str = "set cntrparam levels discrete ";
  str += std::to_string(levels[0]);
  for(size_t i = 1; i < levels.size(); i++)
    str += "," + std::to_string(levels[i]);
  addCommand(str);
}
template void GNUPlotter::setContourLevels(const std::vector<double>& levels);
template void GNUPlotter::setContourLevels(const std::vector<float>& levels);
template void GNUPlotter::setContourLevels(const std::vector<int>& levels);
// ToDo: 
// -Check, if we really need these instantiations. They could be generated automatically from
//  plotContourMap().

// data setup:

void GNUPlotter::setDataPrecision(unsigned int n)
{
  delete[] formatString;
  std::string str = "% 0" + s(n) + "." + s(n) + "le";
  formatString = toZeroTerminatedString(str);
}

template <class T>
void GNUPlotter::addDataBlockLineColumn(const std::vector<std::vector<std::vector<T>>>& d)
{
  std::ofstream out(dataPath, std::ofstream::app);
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
template void GNUPlotter::addDataBlockLineColumn(const std::vector<std::vector<std::vector<double>>>& d);
template void GNUPlotter::addDataBlockLineColumn(const std::vector<std::vector<std::vector<float>>>& d);
template void GNUPlotter::addDataBlockLineColumn(const std::vector<std::vector<std::vector<int>>>& d);

template <class T>
void GNUPlotter::addDataBlockColumnLine(const std::vector<std::vector<std::vector<T>>>& d)
{
  std::ofstream out(dataPath, std::ofstream::app);
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
template void GNUPlotter::addDataBlockColumnLine(const std::vector<std::vector<std::vector<double>>>& d);
template void GNUPlotter::addDataBlockColumnLine(const std::vector<std::vector<std::vector<float>>>& d);
template void GNUPlotter::addDataBlockColumnLine(const std::vector<std::vector<std::vector<int>>>& d);

template <class T>
void GNUPlotter::addData(int numBlocks, int *blockLengths, int numColumns, T **data)
{
  std::ofstream out(dataPath, std::ofstream::app);
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
void GNUPlotter::addDataComplex(const std::vector<std::complex<T>>& d)
{
  std::ofstream out(dataPath, std::ofstream::app);
  for(size_t i = 0; i < d.size(); i++)
  {
    out << sd(d[i].real()) + " " + sd(d[i].imag()) + " ";
    out << "\n";
  }
  out << "\n\n";
  out.close();
  dataInfo.push_back(DataInfo(1, 2));
}
template void GNUPlotter::addDataComplex(const std::vector<std::complex<int>>& d);
template void GNUPlotter::addDataComplex(const std::vector<std::complex<float>>& d);
template void GNUPlotter::addDataComplex(const std::vector<std::complex<double>>& d);

template <class T>
void GNUPlotter::addDataArrays(int N, T *x, int M, T **y)
{
  std::ofstream out(dataPath, std::ofstream::app);
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
  const std::vector<const T*> v = collectLeadingNonNullArguments(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9);
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
  std::vector<T(*)(T)> f = collectLeadingNonNullArguments(f0, f1, f2, f3, f4, f5, f6, f7, f8, f9);

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
  std::ofstream out(dataPath, std::ofstream::app);
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
  std::ofstream out(dataPath, std::ofstream::app);
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
void GNUPlotter::addDataMatrixFlat(int Nx, int Ny, T* x, T* y, T* z)
{
  T** rowPointers = new T*[Nx]; // int **array = new int*[10];
  for(int i = 0; i < Nx; i++)
    rowPointers[i] = (T*) &z[i*Ny]; // eww - casting away const? that's dirty!
  addDataMatrix(Nx, Ny, x, y, rowPointers);
  delete[] rowPointers;
}
// not yet tested

template <class T>
void GNUPlotter::addDataMatrixFlat(int Nx, int Ny, T* z)
{
  T* x = new T[Nx]; for(int i = 0; i < Nx; i++) x[i] = T(i);
  T* y = new T[Ny]; for(int i = 0; i < Ny; i++) y[i] = T(i);
  addDataMatrixFlat(Nx, Ny, x, y, z);
  delete[] x;
  delete[] y;
}
template void GNUPlotter::addDataMatrixFlat(int Nx, int Ny, int* z);
template void GNUPlotter::addDataMatrixFlat(int Nx, int Ny, float* z);
template void GNUPlotter::addDataMatrixFlat(int Nx, int Ny, double* z);

template <class T>
void GNUPlotter::addDataCurve2D(const std::function<T(T)>& fx, const std::function<T(T)>& fy,
  int Nt, T tMin, T tMax, bool writeT)
{
  std::vector<T> t(Nt), x(Nt), y(Nt);
  rangeLinear(&t[0], Nt, tMin, tMax);
  for(int i = 0; i < Nt; i++) {
    x[i] = fx(t[i]);
    y[i] = fy(t[i]);
  }
  if(writeT) addDataArrays(Nt, &t[0], &x[0], &y[0]); // write t-values into 1st column
  else       addDataArrays(Nt,        &x[0], &y[0]); // ...or don't

  // the t values may be useful for using color to indicate the t-value along the 
  // curve? see: http://www.gnuplotting.org/using-a-palette-as-line-color/
  // here for colormaps: https://github.com/Gnuplotting/gnuplot-palettes
  // ..or maybe use tick-marks along the curve - for example, if t in in 0..1 set markers
  // at 0.1, 0.2, 0.9, 1.0 - maybe these markers could be little arrows to also convey the 
  // orientation of the curve? in any case, the markers should not appear on each t - that would
  // be far too dense - maybe every 20th sample or so would be appropriate
  // how about drawing tangent vectors along the curve? make a function addDataCurveWithTangents2D
  // compute the actual tangents numerically
}



template <class T>
void GNUPlotter::addDataSurface(const std::function<T(T, T)>& fx, 
  const std::function<T(T, T)>& fy, const std::function<T(T, T)>& fz,
  int Nu, T uMin, T uMax, int Nv, T vMin, T vMax)
{
  // The outer index runs over the indices for parameter u, the middle index runs over v and the 
  // innermost vector index runs from 0...2 giving a 3-vector containing x, y, z coordinates for 
  // each point:
  std::vector<std::vector<std::vector<T>>> d;    // doubly nested vector of data
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
  for(i = 0; i < Nx; i++) {
    for(j = 0; j < Ny; j++)
      z[i][j] = f(x[i], y[j]); }
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
void GNUPlotter::addDataVectorField2D(const std::function<T(T, T)>& fx, 
  const std::function<T(T, T)>& fy, int Nx, T xMin, T xMax, int Ny, T yMin, T yMax)
{
  int Nv = Nx*Ny;                                // number of vectors to draw
  std::vector<T> 
    x(Nv), y(Nv), dx(Nv), dy(Nv), c(Nv);         // arrays to hold our data
  T xStep = (xMax-xMin) / T(Nx-1);               // step size for x
  T yStep = (yMax-yMin) / T(Ny-1);               // step size for y
  T arrowLength = std::min(xStep, yStep);        // length of arrows to draw
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
void GNUPlotter::addVectorField2D(const std::function<T(T, T)>& fx, 
  const std::function<T(T, T)>& fy, int Nx, T xMin, T xMax, int Ny, T yMin, T yMax)
{
  addDataVectorField2D(fx, fy, Nx, xMin, xMax, Ny, yMin, yMax);
  addGraph(std::string("index ") + std::to_string(dataInfo.size()-1) + 
    std::string(" using 1:2:3:4:5 with vectors head filled size 0.08,15 ls 2 lc palette notitle"));
}

template<class T>
void GNUPlotter::addFieldLine2D(const std::function<T(T, T)>& fx, const std::function<T(T, T)>& fy,
  T x0, T y0, T stepSize, int numPoints, int oversampling)
{
  addDataFieldLine2D(fx, fy, x0, y0, stepSize, numPoints, oversampling);
  addGraph("index " + std::to_string(dataInfo.size()-1) + " using 1:2 with lines lt 1 notitle");
}
template void GNUPlotter::addFieldLine2D(const std::function<float(float, float)>& fx, const std::function<float(float, float)>& fy, float x0, float y0, float stepSize, int numPoints, int oversampling);
template void GNUPlotter::addFieldLine2D(const std::function<double(double, double)>& fx, const std::function<double(double, double)>& fy, double x0, double y0, double stepSize, int numPoints, int oversampling);


//-------------------------------------------------------------------------------------------------
// Drawing:

void GNUPlotter::drawText(const std::string& attr,
  const std::string& text, double x, double y)
{
  addCommand("set label \"" + text + "\" at " + s(x) + "," + s(y) + " " + attr);
}
// https://stackoverflow.com/questions/16820963/how-to-add-text-to-an-arrow-in-gnuplot
// http://www.manpagez.com/info/gnuplot/gnuplot-4.4.3/gnuplot_259.php

void GNUPlotter::drawLine(const std::string& attributes,
  double x1, double y1, double x2, double y2)
{
  drawArrow("nohead " + attributes, x1, y1, x2, y2);
  // A line is just drawn as an arrow without head, i.e. the nohead attribute is added to the 
  // attributes that are already given.
}

void GNUPlotter::drawArrow(const std::string& attr,
  double x1, double y1, double x2, double y2)
{
  addCommand("set arrow from " + s(x1) + "," + s(y1) + " to "
    + s(x2) + "," + s(y2) + " " + attr);
}

void GNUPlotter::drawPolyLine(const std::string& attributes, 
  const std::vector<double>& x, const std::vector<double>& y)
{
  assume(x.size() == y.size(), "x and y must have the same size");
  for(int i = 0; i < (int)x.size() - 1; i++)
    drawLine(attributes, x[i], y[i], x[i+1], y[i+1]);
  // It's awkward that we have to draw a polyline as a bunch of lines which in turn are just arrows
  // without a head, but gnuplot doesn't provide line or polyline primitives. ...how weird is that?
}

void GNUPlotter::drawPolygon(const std::string& attributes,
  const std::vector<double>& x, const std::vector<double>& y)
{
  assume(x.size() == y.size(), "x and y must have the same size");
  if(x.size() < 3) return; // we don't draw degenerate polygons
  std::string cmd = "set object polygon from ";
  for(size_t i = 0; i < x.size(); i++)
    cmd += s(x[i]) + "," + s(y[i]) + " to ";
  cmd += s(x[0]) + "," + s(y[0]) + " " + attributes;
  addCommand(cmd);
  // The syntax for polygons, circles, ellipses, etc. is different than for arrows and (text) 
  // labels. It doesn't includes "object" - that's a bit inconsistent.
}

void GNUPlotter::drawCircle(const std::string& attr, double x, double y, double r)
{
  addCommand("set object circle at " + s(x) + "," + s(y) + " size " + s(r) + " " + attr);
}

void GNUPlotter::drawEllipse(const std::string& attr, 
  double x, double y, double w, double h, double a)  // angle is in degrees
{
  addCommand("set object ellipse center " + s(x) + "," + s(y) + " size " 
              + s(w) + "," + s(h) + " angle " + s(a) + " " + attr);
}
// http://gnuplot.sourceforge.net/demo/ellipse.html

// see here for what types of objects are supported:
// http://soc.if.usp.br/manual/gnuplot-doc/htmldocs/set_002dshow.html#set_002dshow

//-------------------------------------------------------------------------------------------------
// inquiry:

std::string GNUPlotter::getDataPath()
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

//-------------------------------------------------------------------------------------------------
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

void GNUPlotter::showMultiPlot(int numRows, int numCols, const std::string& how)
{
  addCommand("set multiplot");                             // init multiplot

  double h = 1.0 / numRows;                                // relative height of subplots
  double w = 1.0 / numCols;                                // relative width of subplots
  for(int i = 0; i < numRows; i++)                         // loop over the plot-rows
    for(int j = 0; j < numCols; j++)                       // loop over the plot-columns
      addSubPlot(j*w, 1-i*h-h, w, h, numCols*i+j, how);

  addCommand("unset multiplot");
  invokeGNUPlot();
}
// factor out initMultiPlot/finishMultiPlot...or maybe use show instead of finish

void GNUPlotter::addSubPlot(double x, double y, double w, double h, int datasetIndex, 
  const std::string& how)
{
  addCommand("set origin " + s(x) + "," + s(y));    // bottom-left corner of this subplot
  addCommand("set size "   + s(w) + "," + s(h));    // size of this subplot
  addCommand("plot '" + getDataPath() + "' i " + s((unsigned int)datasetIndex) + " " + how);
}
// todo: this should also allow to use splot instead of plot - have a boolean splot or plot3D 
// parameter

void GNUPlotter::clearCommandFile()
{
  initFile(commandPath);
}

void GNUPlotter::invokeGNUPlot()
{
  // create the callstring and invoke GNUPlot:
  //string callString = "\"" + gnuplotPath + "\" " + commandPath + " -";
  //string callString = "\"" + gnuplotPath + "\" " + commandPath;
  std::string callString = "\"" + gnuplotPath + "\" " + commandPath + " -persist";

  // wrapping gnuplotPath into quotes is required to handle installation paths with whitespaces,
  // but it doesn't seem to be possible to handle a commandPath with whitespaces (i tried
  // wrapping it into quotes too and wrapping the whole callString into quotes - none of that
  // works)
  // the minus at the end prevents gnuplot from immediately closing

  systemCall(callString);
  // is it possible to call GNUPlot in a separate process? - this is actually already the case
  // ...but we have to close it to call it again - i.e. we can't do multiple plots at once
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void GNUPlotter::initFile(const std::string &path)
{
  std::ofstream out(path);
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

void GNUPlotter::assume(bool condition, const char* errorMessage)
{
  if(!condition) {
    std::cout << errorMessage;
    std::terminate(); }
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
  return std::to_string(x);
}

// internal functions

void GNUPlotter::addPlotCommand(bool splot)
{
  // If dataInfo is empty (i.e. no data was added to the datafile by client code), we generate some
  // dummy data. This is needed, if the user just wants to draw geometric figures on an empty 
  // canvas. If there is no data at all, the code below crashes and even modifying it to avoid the 
  // crash doesn't work - then it doesn't crash GNUPlot doesn't open - so, we need at least one 
  // dummy dataset to plot:
  if(dataInfo.empty()) {
    double dummy = 0; addDataArrays(1, &dummy, &dummy, &dummy); } // may work for 2D and 3D

  // hmm - maybe in case of an empty dataset, we should just call replot?:
  // http://soc.if.usp.br/manual/gnuplot-doc/htmldocs/set_002dshow.html#set_002dshow



  // Auto-generate graph-descriptors, if user has not set them up manually:
  if(graphDescriptors.empty())
    generateGraphDescriptors(splot);

  // Initialize the plot command:
  int i;
  std::string pc;
  addCommand("\n# Plotting:");
  if( splot == true )
    pc = "splot \\\n";  // 3D plots
  else
    pc = "plot \\\n";   // 2D plots

  // Add the graph-descriptors to the plot command:
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
    for(i = 0; i < dataInfo.size(); i++) {
      if(dataInfo[i].type == "matrix") {
        addGraph("i " + s(i) + " nonuniform matrix" + getStyleString(k));
        k++; }
      else {
        if(dataInfo[i].numColumns == 1) { // this is new -> test...
          //addGraph("i " + s(i) + " u 1:1" + getStyleString(k));
          addGraph("i " + s(i) + getStyleString(k));
          k++; }
        // entered only, if condition above is false:
        for(j = 0; j < (unsigned int)dataInfo[i].numColumns-1; j++) {
          addGraph("i " + s(i) + " u 1:" + s(j+2) + getStyleString(k));
          k++;  }}}}
  else {
    for(i = 0; i < dataInfo.size(); i++) {
      if(dataInfo[i].type == "matrix")
        addGraph("i " + s(i) + " nonuniform matrix" + getStyleString(i));
      else
        addGraph("i " + s(i) + getStyleString(i)); }}
}
// refactor this - split out two functions generateGraphDescriptors2D/3D


// argument handling

template<class T>
T GNUPlotter::nullValue(T)
{
  return T(0);
}

std::string GNUPlotter::nullValue(std::string)
{
  return "";
}

template<class T>
const std::vector<T> GNUPlotter::collectLeadingNonNullArguments(const T a0, const T a1, const T a2, 
  const T a3, const T a4, const T a5, const T a6, const T a7, const T a8, const T a9)
{
  T null = nullValue(a0);
  std::vector<T> v;
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
void GNUPlotter::append(std::vector<T>& v, const std::vector<T>& appendix)
{
  v.reserve(v.size() + appendix.size());
  for(size_t i = 0; i < appendix.size(); i++)
    v.push_back(appendix[i]);
}

template<class T>
std::vector<std::vector<T>> GNUPlotter::wrapIntoVectors(int N, const T *a0, const T *a1, const T *a2, 
  const T *a3, const T *a4, const T *a5, const T *a6, const T *a7, const T *a8, const T *a9)
{
  std::vector<T*> pointers = collectLeadingNonNullArguments(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
  std::vector<std::vector<T>> v;
  v.resize(pointers.size());
  for(int i = 0; i < v.size(); i++)
  {
    v[i].resize(N);
    for(int j = 0; j < N; j++)
      v[i][j] = pointers[i][j];
  }
  return v;
}

void GNUPlotter::addToStringVector(std::vector<std::string>& v, CSR s0, CSR s1, CSR s2, CSR s3, 
  CSR s4, CSR s5, CSR s6, CSR s7, CSR s8, CSR s9)
{
  append(v, collectLeadingNonNullArguments(s0, s1, s2, s3, s4, s5, s6, s7, s8, s9));
}

void GNUPlotter::setStringVector(std::vector<std::string>& v, CSR s0, CSR s1, CSR s2, CSR s3, 
  CSR s4, CSR s5, CSR s6, CSR s7, CSR s8, CSR s9)
{
  v.clear();
  addToStringVector(v, s0, s1, s2, s3, s4, s5, s6, s7, s8, s9);
}

/*=================================================================================================

ToDo:

-maybe move the explicit template instantiations to another file...that would reduce clutter in 
 this implementation file - but would make the library harder to use - the user would have to deal
 with more files...so it's probably not such a good idea...simple use is more important than
 nice looking code

 there's an alternative way for doing this stuff - look into this:
 https://vijaypolimeru.github.io/research/plotting/programming/graphics/opensees/visual%20studio/2020/02/05/plotting-in-cpp.html

*/
