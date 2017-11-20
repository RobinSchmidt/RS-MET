//// code has become obsolete - delete soon
//
//#include "PlotterExperiments.h"
//
//void createDataFiveFunctions(int N, double xMin, double xMax, double *x, double *y1, double *y2,
//  double *y3, double *y4, double *y5)
//{
//  rsFillWithRangeLinear(x, N, xMin, xMax);
//  for(int n = 0; n < N; n++)
//  {
//    y1[n] = sin(x[n]);
//    y2[n] = cos(x[n]);
//    y3[n] = y1[n] + y2[n];
//    y4[n] = tanh(x[n]);
//    y5[n] = 1 / (1 + x[n]*x[n]);
//  }
//}
//
//void testDataPlot2()
//{
//  // a simple plot of 5 datasets representing functions using the default settings
//
//  // create data:
//  static const int N = 200;  // number of data points
//  double x[N], y1[N], y2[N], y3[N], y4[N], y5[N];
//  createDataFiveFunctions(N, -10, +10, x, y1, y2, y3, y4, y5);
//
//  // plot data:
//  GNUPlotter p;
//  p.addData(N, x, y1, y2, y3, y4, y5);
//  p.plot();
//    // write convenience function, so we can do this with one line
//
//  int dummy = 0;
//}
//
//void testIntVector()
//{
//  // create the Fibonacci sequence in y, write the index into x:
//  int N = 10;
//  vector<int> x, y;
//  x.resize(N);
//  y.resize(N);
//  x[0] = y[0] = 0;
//  x[1] = y[1] = 1;
//  for(int i = 2; i < N; i++)
//  {
//    x[i] = i;
//    y[i] = y[i-1] + y[i-2];
//  }
//
//  GNUPlotter p;
//  //p.addData(x, y);
//  p.addData(y);
//  p.setGraphStyles("points pt 7 ps 1");
//  p.setGrid();
//  p.plot();
//
//  int dummy = 0;
//}
//
//void testStyleSetup2()
//{
//  // create data:
//  static const int N = 35;  // number of data points per set
//  double x[N], y1[N], y2[N], y3[N], y4[N], y5[N];
//  createDataFiveFunctions(N, -5, +5, x, y1, y2, y3, y4, y5);
//
//  // create plotter object, pass the data: 
//  GNUPlotter p;
//  p.addData(N, x, y1, y2, y3, y4, y5);
//
//  // set up the plotter:
//  p.setLegends("sin(x)", "cos(x)", "sin(x)+cos(x)");  // graphs 1,2,3 have legends
//  p.setGraphColors("800000", "008000", "000080");     // 3 colors, used cyclically
//  p.setGraphStyles("lines lw 2", "points pt 6 ps 2"); // styles for graphs 1 and 2
//  p.addGraphStyles("impulses lw 2", "boxes lw 1.5");  // styles for graphs 3 and 4
//  p.addGraphStyles("linespoints pt 7 ps 1");          // style for graph 5
//  p.setSize(800, 400);                                // pixel size for plot
//  p.setRange(-4, +4, -1.6, +1.6);                     // range for x- and y-axis
//  p.setAxisLabels("x-axis", "y-axis");                // labels for the axes
//  p.setTitle("Plotter Demo");                         // caption for the plot
//  p.setGrid(true, false, false, true); 
//    // the fine grid does not yet work - maybe the tics must be specified before? via 
//    // set x2tics, etc. 
//  //p.setRange(0.1, 5, 0.1, +1.6);  p.setLogScale("xy", 2);  // test
//
//  // plot:
//  p.plot();
//}
//
//void testMultiColumn2()
//{
//  static const int N = 30;      // maximum number of datapoints
//  double d[17][N];              // data - columns of N values
//
//  double a      = 0.3;          // multplier for the sine's argument
//  double b      = 1.0;          // y-axis offset between datasets
//
//  // create the data for the datasets:
//  int i;
//  for(i = 0; i < N; i++)
//  {
//    double x = i;              
//    double y = sin(a*x);
//
//    d[0][i] = x;           // x for all datasets (except 1-column which has implicit x)
//
//    d[1][i] = y;           // y for 1-column dataset
//
//    y += b;
//    d[2][i] = y;           // y for 2-column dataset
//
//    y += b;
//    d[3][i] = y;           // y for 3-column dataset (for yerrorlines)
//    d[4][i] = 0.3;         // y-delta
//
//    y += b;
//    d[5][i] = y;           // y for 4-column dataset (for yerrorlines)
//    d[6][i] = y - 0.4;     // y min
//    d[7][i] = y + 0.2;     // y max
//
//    y += b;
//    d[8][i]  = y;          // box min for 5-column dataset (for candlesticks)
//    d[9][i]  = y - 0.2;    // whisker min
//    d[10][i] = y + 0.4;    // whisker max
//    d[11][i] = y + 0.2;    // box max
//
//    y += b;
//    d[12][i] = y;          // y for 6-column dataset (for xyerrorlines)
//    d[13][i] = x - 0.4;    // x min
//    d[14][i] = x + 0.5;    // x max
//    d[15][i] = y - 0.2;    // y min
//    d[16][i] = y + 0.3;    // y max
//  }
//
//  // Create plotter and add datasets. We add the 1, 2 and 3 column datasets at once with a single 
//  // call and then add the 4, 5 and6 column datasets each seperately with their own call. Note that
//  // the datasets need not have the same length (i.e. the same number of rows). This is 
//  // demonstrated by shortening the 3-column dataset by 5 rows:
//  GNUPlotter p;
//  //p.addDataMultiColumn(N, 1, d[1], N, 2, d[0], d[2], N-5, 3, d[0], d[3], d[4], 0); // 1,2,3 columns
//  //p.addDataMultiColumn(N, 4, d[0], d[5], d[6], d[7], 0);                           // 4 columns
//  //p.addDataMultiColumn(N, 5, d[0], d[8], d[9], d[10], d[11], 0);                   // 5 columns
//  //p.addDataMultiColumn(N, 6, d[0], d[12], d[13], d[14], d[15], d[16], 0);          // 6 columns
//
//  p.addDataMultiColumn(N, 1, d[1]);                                    // 1 column
//  p.addDataMultiColumn(N, 2, d[0], d[2]);                              // 2 columns
//  p.addDataMultiColumn(N, 3, d[0], d[3], d[4]);                        // 3 columns
//  p.addDataMultiColumn(N, 4, d[0], d[5], d[6], d[7]);                  // 4 columns
//  p.addDataMultiColumn(N, 5, d[0], d[8], d[9], d[10], d[11]);          // 5 columns
//  p.addDataMultiColumn(N, 6, d[0], d[12], d[13], d[14], d[15], d[16]); // 6 columns
//
//
//  // set up styles and plot:
//  p.setRange(0, 30, -2, 7);
//  //p.setGraphStyles("lines", "lines", "yerrorlines", "yerrorlines", "candlesticks", 
//  //  "xyerrorlines", 0);
//
//  p.setLegends("lines", "linespoints", "yerrorlines", "yerrorlines", "candlesticks", 
//    "xyerrorlines");
//  // try to place that bottom right
//
//  p.setGraphStyles("lines", "linespoints", "yerrorlines", "yerrorlines", "candlesticks", 
//    "xyerrorlines");
//  p.plot();
//}
//
//void testSurface2()
//{
//  static const int Nx = 41;    // number of x-values
//  static const int Ny = 41;    // number of y-values
//
//  double xMin = -10;
//  double xMax = +10;
//  double yMin = -10;
//  double yMax = +10;
//
//  double x[Nx];
//  double y[Ny];
//  double z1[Nx][Ny], z2[Nx][Ny];
//
//  int i, j;
//  double xi, yj;
//  double *z1p[Nx], *z2p[Nx];
//  for(i = 0; i < Nx; i++)
//  {
//    z1p[i] = &z1[i][0];
//    z2p[i] = &z2[i][0];
//  }
//
//  for(i = 0; i < Nx; i++) x[i] = xMin + i*(xMax-xMin)/(Nx-1);
//  for(j = 0; j < Ny; j++) y[j] = yMin + j*(yMax-yMin)/(Ny-1);
//  for(i = 0; i < Nx; i++)
//  {
//    xi = x[i];
//    for(j = 0; j < Ny; j++)
//    {
//      yj = y[j];
//      z1[i][j] = exp(-0.03*(xi*xi+yj*yj));       // 2D Gaussian bell
//      z2[i][j] = 2 + 1 / (1+0.1*(xi*xi+yj*yj));  // 2D Butterworth bell shifted
//    }
//  }
//
//  GNUPlotter p;
//  p.addDataSurface(Nx, Ny, x, y, z1p);
//  //p.addDataSurface(Nx, Ny, x, y, z2p); // works, but doesn't look nice
//
//  p.setGrid();
//  p.addCommand("set hidden3d"); // add function removeHiddenLines or setHiddenLineRemoval
//  p.plot(true);
//
//  // todo check out how to do contour lines (or maybe fills), etc.
//}
//
//
//// here is a good tutorial:
//// http://lowrank.net/gnuplot/index-e.html
//
//// demo for various line styles:
//// http://gnuplot.sourceforge.net/demo_4.6/dashcolor.html
//
//// for polar grids
//// http://stackoverflow.com/questions/6772135/how-to-get-a-radialpolar-plot-using-gnu-plot
//
//// for 3D plot, manual, pages 65,115,184
//// http://lowrank.net/gnuplot/datafile-e.html
//
//
//// for example plots: 
//// -Integer plot: fibonacci sequence (or prime numbers) - maybe we can have a special int plot
////  that writes the values above the datapoints (centered) - but maybe this is something for a 
////  subclass that deals specifically with integer plots
//// -Gaussian bells with diffent mu, sigma - also demonstrates how to use greek letters
//// -square and/or saw wave with linear, cubic and sinc interploation - use points and impulses for
////  the discrete time data, lines of different colors for the interpolants
//// -pole/zero plot of an elliptic bandpass in the z-domain
//// -for log/log plots: familiy of Butterworth magnitude responses
////  -maybe with phase-plots in the same plot (angles written on right y-axis, we need an expression 
////   for the (unwrapped) phase
//// -3D: 2D Gaussian, z- or s-domain pole/zero "landscapes", clouds of points (maybe using gaussian
////  distributions)
//// -spectrogram-plots, phasogram-plots
//
//
////=================================================================================================
//// old:
//
//std::vector<std::string> charPointersToStringVector(int numStrings, ...)
//{
//  std::vector<std::string> v;
//  va_list argList;
//  va_start(argList, numStrings);
//  for(int i = 0; i < numStrings; i++)
//    v.push_back(va_arg(argList, const char*));
//  va_end(argList);
//  return v;
//}
//
//std::vector<std::string> stringList(int numStrings, const char *string1, ...)
//{
//  std::vector<std::string> v;
//  v.push_back(string1);
//  va_list argList;
//  va_start(argList, numStrings);
//  for(int i = 1; i < numStrings; i++)
//    v.push_back(va_arg(argList, const char*));
//  va_end(argList);
//  return v;
//}
//
//void testVariableArgumentList()
//{
//  //std::vector<std::string> v = charPointersToStringVector(3, "s1", "s2", "s3");
//
//  std::vector<std::string> v = stringList(3, "s1", "s2", "s3");
//  int dummy = 0;
//}
//
//void testDataPlot1()
//{
//  // a simple plot of 5 datasets representing functions using the default settings
//
//  // create data:
//  static const int N = 1000;  // number of data points
//  double x[N], y1[N], y2[N], y3[N], y4[N], y5[N];
//  createDataFiveFunctions(N, -10, +10, x, y1, y2, y3, y4, y5);
//
//  // plot data:
//  Plotter::plotFunctionData(5, N, x, y1, y2, y3, y4, y5);
//}
//
//void testColorSetup()
//{
//  // a simple plot of 5 datasets representing functions using custom colors
//
//  // create data:
//  static const int N = 1000;  // number of data points
//  double x[N], y1[N], y2[N], y3[N], y4[N], y5[N];
//  createDataFiveFunctions(N, -10, +10, x, y1, y2, y3, y4, y5);
//
//  // create plotter object, pass the data: 
//  Plotter plt;
//  plt.setFunctionData(5, N, x, y1, y2, y3, y4, y5);
//
//  // set up the colors (we specify only 3 colors but have 5 datasets, so we will see how the 
//  // defined colors will be used cyclically):
//  std::vector<std::string> colors;
//  colors.push_back("770000");
//  colors.push_back("007700");
//  colors.push_back("000077");
//  plt.setGraphColors(colors);
//
//  // alternatively to the code above, we could use the more convenient form:
//  plt.setGraphColors(3, "770000", "007700", "000077");
//
//  // plot:
//  plt.plot();
//}
//
//void testStyleSetup()
//{
//  // rename function to demoPlotterSetup (or similar)
//
//  // create data:
//  static const int N = 35;  // number of data points per set
//  double x[N], y1[N], y2[N], y3[N], y4[N], y5[N];
//  createDataFiveFunctions(N, -5, +5, x, y1, y2, y3, y4, y5);
//
//  // create plotter object, pass the data: 
//  Plotter p;
//  p.setFunctionData(5, N, x, y1, y2, y3, y4, y5);
//
//  // set up the plotter:
//  p.setSize(800, 400);
//  p.setRange(-4, +4, -1.6, +1.6);
//  p.setGraphColors(3, "770000", "007700", "000077");
//  p.setGraphStyles(5, "lines lw 2", "points pt 6 ps 2", "impulses lw 2", "boxes lw 1.5", 
//    "linespoints pt 7 ps 1");
//  p.setGraphTitles(2, 1, "sin(x)", 3, "sin(x)+cos(x)");
//  p.setAxisLabels("x", "y");
//  p.setTitle("Plotter Demo");
//
//  // to do:
//  // p.addAnnotation("set arrow from 1,1 to 0.5,0.5");
//  // p.addAnnotation("set label \"Label 1\" at 2,1");
//
//
//  // plot:
//  p.plot();
//
//  // move to function testPointTypes
//  //plt.setGraphStyles(5, "points pt 1", "points pt 2", "points pt 3", "points pt 4", "points pt 5");
//  //plt.setGraphStyles(5, "points pt 6", "points pt 7", "points pt 8", "points pt 9", "points pt 10");
//
//  // pointtypes: 1: +, 2: crosses, 3: *, 4: squares, 5: filled squares, 6: circles, 
//  // 7: filled circles, 8: triangle up, 9: filled triangle up, 10: triangle down,
//
//  // move to function testLineStyles:
//  //plt.setGraphStyles(4, "lines", "lines lw 2", "points pt 6 ps 2", "linespoints pt 7 ps 1");
//  //plt.setGraphStyles(5, "lines", "points pt 2", "impulses", "boxes", "linespoints pt 7 ps 1");
//  //plt.setGraphStyles(4, "lines lt 1", "lines lt 2", "lines lt 3", "lines lt 3");
//    // should create dashed lines but doesn't work
//}
//
//void testMultiColumn()
//{
//  // We create datasets with different numbers of columns 1...6 and plot each dataset with a
//  // style that is appropriate for that number of columns. The datasets will be offsetted sine
//  // functions plus some "error tolerances" where applicable.
//
//  static const int numSets = 6;                         // number of data sets
//  int numPoints[numSets] = { 25, 20, 30, 23, 27, 28 };  // number of data points per set
//  double a      = 0.3;          // multplier for the sine's argument
//  double offset = 1.0;          // y-axis offset between datasets
//  int i, j, k;                  // loop indices
//  double x, y, z1, z2, z3, z4;  // coordinates
//
//  // create data vector and set up the dimensionalities:
//  std::vector<std::vector<std::vector<double>>> data;
//  data.resize(numSets);
//  for(i = 0; i < numSets; i++)
//  {
//    data[i].resize(numPoints[i]);
//    for(j = 0; j < data[i].size(); j++)
//      data[i][j].resize(i+1);
//  }
//
//  // create dataset 1 (for plot with lines):
//  for(j = 0; j < data[0].size(); j++)
//  {
//    x = j;                    // x (implicit - not stored in dataset)
//    y = sin(a*x);             // y
//    data[0][j][0] = y;
//  }
//
//  // create dataset 2 (for plot with lines):
//  for(j = 0; j < data[1].size(); j++)
//  {
//    x = j;                    // x
//    y = sin(a*x) + 1*offset;  // y 
//    data[1][j][0] = x;
//    data[1][j][1] = y;
//  }
//
//  // create dataset 3 (for plot with yerrorlines):
//  for(j = 0; j < data[2].size(); j++)
//  {
//    x  = j;                   // x
//    y  = sin(a*x) + 2*offset; // y
//    z1 = 0.2;                 // y delta
//    data[2][j][0] = x;
//    data[2][j][1] = y;
//    data[2][j][2] = z1;
//  }
//
//  // create dataset 4 (for plot with yerrorlines):
//  for(j = 0; j < data[3].size(); j++)
//  {
//    x = j;                   // x
//    y = sin(a*x) + 3*offset; // y
//    z1 = y - 0.3;            // y min
//    z2 = y + 0.2;            // y max
//    data[3][j][0] = x;
//    data[3][j][1] = y;
//    data[3][j][2] = z1;
//    data[3][j][3] = z2;
//  }
//
//  // create dataset 5 (for plot with candlesticks):
//  for(j = 0; j < data[4].size(); j++)
//  {
//    x = j;                   // x
//    y = sin(a*x) + 4*offset; // box min
//    z1 = y - 0.2;            // whisker min
//    z2 = y + 0.4;            // whisker max
//    z3 = y + 0.2;            // box max
//    data[4][j][0] = x; 
//    data[4][j][1] = y;
//    data[4][j][2] = z1;
//    data[4][j][3] = z2;
//    data[4][j][4] = z3; 
//  }
//
//  // create dataset 5 (for plot with xyerrorlines):
//  for(j = 0; j < data[5].size(); j++)
//  {
//    x = j;                   // x
//    y = sin(a*x) + 5*offset; // y
//    z1 = x - 0.4;            // x min
//    z2 = x + 0.5;            // x max
//    z3 = y - 0.2;            // y min
//    z4 = y + 0.3;            // y max
//    data[5][j][0] = x;
//    data[5][j][1] = y;
//    data[5][j][2] = z1;
//    data[5][j][3] = z2;
//    data[5][j][4] = z3;
//    data[5][j][5] = z4;
//  }
//
//  // create plotter, set up styles, pass the data and plot:
//  Plotter p;
//  p.setRange(0, 30, -2, 7);
//  p.setGraphStyles(6, "lines", "lines", "yerrorlines", "yerrorlines", "candlesticks", 
//    "xyerrorlines");
//  p.setData(data);
//  p.plot();
//}
//
//
//void testSurface()
//{
//  // does not yet work - see comment below
//
//  static const int Nx = 10;    // number of x-values
//  static const int Ny = 15;    // number of y-values
//  static const int N  = Nx*Ny; // number of datapoints
//
//  double xMin = -10;
//  double xMax = +10;
//  double yMin = -10;
//  double yMax = +10;
//
//  int i, j, k;
//  double x, y, z;
//
//  std::vector<std::vector<std::vector<double>>> data;
//  data.resize(1);
//  data[0].resize(N);
//  for(j = 0; j < data[0].size(); j++)
//    data[0][j].resize(3);
//
//  for(j = 0; j < Nx; j++)
//  {
//    x = j; // ...later scale
//
//    for(k = 0; k < Ny; k++)
//    {
//      y = k;
//      //z = 1 / (1+x*y*x*y);
//      z = 1 / (1+x*x+y*y);
//      i = j*Ny+k;         // index for x value
//      data[0][i][0] = x;
//      data[0][i][1] = y;
//      data[0][i][2] = z;
//
//      int dummy = 0;
//    }
//  }
//
//  int dummy = 0;
//
//  // for a 3D plot, we must seperate the blocks of data (for one particular x value) by single 
//  // blank lines
//  // introduce a blockSizes array for blocking the data inside each dataset
//
//  Plotter p;
//  p.setSurfaceMode(true);
//  //p.addCommand("set contour");
//  p.setData(data);
//  p.plot();
//}
//
//
//
