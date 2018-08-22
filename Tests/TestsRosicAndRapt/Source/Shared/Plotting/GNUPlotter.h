#ifndef GNUPLOTTER_H
#define GNUPLOTTER_H

#include <string>
#include <vector>
#include <complex>

/**

This class allows for plotting data and functions directly from C++ code by invoking GNUPlot.
You can pass in the data to be plotted and set up a couple of things that control how the data will
be presented. From this data and the setup, an object of this class will create a datafile and a 
batchfile with GNUPlot commands. Then, GNUPlot is invoked from the command line and told to read 
the batchfile, which will do the setup tasks and then tell GNUPlot to actually plot the data 
stored in the datafile.

For this to work, the following preliminaries must be met:
-GNUPlot must be installed in the path defined by the member variable gnuplotPath. Currently, this 
 is fixed to C:/Program Files/gnuplot/bin/gnuplot.exe, which is the default installation path on 
 windows
-the directory defined in the member variables datapath and commandPath must exist (currently 
 fixed to E:/Temp/) - that's where the temporary data- and command batchfile will be written

*/

typedef const std::string & CSR;     // type "CSR" is a const reference to a std::string
                                     // maybe rename to CRS
template<class T> 
using CVR = const std::vector<T> &;  // type CVR<T> is a const reference to a std::vector of type T
                                     // this is a template alias (C++11)

class GNUPlotter
{

public:

  /** \name Construction/Destruction */

  /** Constructor.*/
  GNUPlotter();

  /** Destructor. */
  virtual ~GNUPlotter();

  /** Initializes our data- and commandfiles. Called from constructor, but you may also call it 
  manually, if you want to re-initialize everything. */
  void initialize();

  /** \name Plotting */

  /** After the data has been set up and possibly a couple of formatting functions have been
  called, calling this function will actually invoke GNUPlot for plotting the data according to the
  desired settings. Used internally, but can be also called, if you set up the data manually. */
  void plot();

  /** Like plot, but is used for 3D plots (invoking GNUPlot's "splot" command instead of "plot"). 
  It is supposed that the appropriate data-setup functions have been called before, such that the 
  datasets may meaningfully be interpreted as 3D datasets.  Used internally, but can be also called, 
  if you set up the data manually.*/
  void plot3D();

  /** Plots the function values in the arrays y1, y2, ... against an abscissa given by the array x. */
  template <class T>
  void plotFunctionTables(int N, T *x, T *y1, T *y2 = nullptr, T *y3 = nullptr, T *y4 = nullptr, 
    T *y5 = nullptr, T *y6 = nullptr, T *y7 = nullptr, T *y8 = nullptr, T *y9 = nullptr);

  /** Plots the values in the arrays y1, y2, ... against an abscissa given by the array index. */
  template <class T>
  void plotArrays(int N, T *y1, T *y2 = nullptr, T *y3 = nullptr, T *y4 = nullptr, T *y5 = nullptr,
    T *y6 = nullptr, T *y7 = nullptr, T *y8 = nullptr, T *y9 = nullptr);

  /** Plots univariate functions. For example, assume you have an array of 100 x-values and you 
  want to plot the sine, cosine and tangent values of these x-values, using the x-array as 
  abscissa, you would call it like: plotFunctions(100, x, &sin, &cos, &tan). */
  template <class T>
  void plotFunctions(int N, T *x, T (*f0)(T), T (*f1)(T) = nullptr, T (*f2)(T) = nullptr,
    T (*f3)(T) = nullptr, T (*f4)(T) = nullptr, T (*f5)(T) = nullptr, T (*f6)(T) = nullptr, 
    T (*f7)(T) = nullptr, T (*f8)(T) = nullptr, T (*f9)(T) = nullptr);

  /** Like above, but you don't need to pass x-axis values. Instead, you pass a minimum and maximum
  value for the x-axis from which the axis values will be created internally. */
  template <class T>
  void plotFunctions(int N, T xMin, T xMax, T (*f0)(T), T (*f1)(T) = nullptr, T (*f2)(T) = nullptr,
    T (*f3)(T) = nullptr, T (*f4)(T) = nullptr, T (*f5)(T) = nullptr, T (*f6)(T) = nullptr, 
    T (*f7)(T) = nullptr, T (*f8)(T) = nullptr, T (*f9)(T) = nullptr);

  /** Plots a surface above the xy-plane, given by a two-dimensional function z = f(x,y). You pass 
  the arrays of x- and y values (which should be of length Nx and Ny respectively) and a matrix of 
  z-values such that z[i][j] = f(x[i], y[j]). */
  template <class T>
  void plotSurface(int Nx, int Ny, T *x, T *y, T **z);

  /** Plots the bivariate function f for x- and y-axis values given in x and y (of length Nx, Ny
  respectively). */
  template <class T>
  void plotBivariateFunction(int Nx, int Ny, T *x, T *y, T (*f)(T, T));

  /** Like above, but instead of giving the x- and y-axis values as arrays, you give minimum and
  maximum values. */
  template <class T>
  void plotBivariateFunction(int Nx, T xMin, T xMax, int Ny, T yMin, T yMax, T (*f)(T, T));
   // let these functions take multiple functions, use internally a function
   // addBivariateFunctionData

  // provide more functions for specialized plots
  // plotFunctionFamily, plotVectorField, plotHistogram, bodePlot, scatterPlot, plotComplexMapping


  /** \name Style Setup */

  /** Adds the given command to the batchfile. */
  void addCommand(std::string command);

  /** Adds the default commands to the command file. This function controls the default appearance
  of plots in cases, when the user doesn't add their own commands. */
  void addDefaultCommands();

  /** Sets the labels for the x-, y- and z-axis. */
  void setAxisLabels(std::string x, std::string y, std::string z = "");

  /** Sets legends for the graphs. You can at most set legends for 10 graphs with this function. If
  you need more than that, call it once with the 1st 10 legends and then add more using subsequent 
  calls of addLegends. */
  void setLegends(CSR l0 = "y = f(x)",  CSR l1 = "", CSR l2 = "", CSR l3 = "", CSR l4 = "", 
    CSR l5 = "", CSR l6 = "", CSR l7 = "", CSR l8 = "", CSR l9 = "");

  /** This function can be used, if more than 10 legends are needed. */
  void setLegends(CVR<std::string> legends);

  /** Sets the colors to be used for the datasets. Each color must be a zero terminated c-string of
  the form RRGGBB or AARRGGBB where RR is a hexadecimal value for the red component, GG for green, 
  BB for blue and the optional AA for an opacity ("alpha") value. If you have more datasets than
  defined colors, GNUPlot's default colors will be used for the additional sets.  */
  void setGraphColors(CSR c0 = "000000", CSR c1 = "", CSR c2 = "", CSR c3 = "", CSR c4 = "", 
    CSR c5 = "", CSR c6 = "", CSR c7 = "", CSR c8 = "", CSR c9 = "");

  /** This function can be used, if more than 10 graph colors are needed. */
  void setGraphColors(CVR<std::string> colors);

  /** Sets the color for the graph with given index (using the "set linetype" command). Indices 
  start at 1. */
  void setGraphColor(unsigned int index, CSR color);

  /** Sets the dash pattern for the graph with given index. The string to specify the pattern can
  be either of the form: (s1,e1,s2,e2,s3,e3,...) where the s-values are the solid lengths and 
  e-values the empty space lengths or a string containing the characters: '.', '-', '_', ' '. 
  Indices start at 1. */
  void setDashType(unsigned int index, CSR type);

  /** Sets the plotting styles to be used for the datasets. Each string for a style must follow
  the GNUPlot syntax for the specification of plot styles. The most important ones are: 
  "lines":       connects all datapoints with straight lines
  "points":      draws a point (of some selectable shape) for each datapoint
  "linespoints": draws points at the datapoints and connexts them with lines
  "impulses":    draws a vertical line from the x-axis to the y-value
  "boxes":       draws a box for each datapoint
  Some of the styles allow you specify further options, for example, the points style lets you 
  select a pointtype and pointsize via "pt" and "ps". For example "points pt 6 ps 2" means to use
  pointtype 6 (which typically is an empty circle) with a size of 2. The first 10 pointtypes are
  1: +, 2: crosses, 3: *, 4: squares, 5: filled squares, 6: circles, 7: filled circles, 
  8: triangle up, 9: filled triangle up, 10: triangle down. For more point types, plotting styles 
  and their respective options, refer to the GNUPlot manual. */
  void setGraphStyles(CSR s0 = "lines", CSR s1 = "", CSR s2 = "", CSR s3 = "", CSR s4 = "", 
    CSR s5 = "", CSR s6 = "", CSR s7 = "", CSR s8 = "", CSR s9 = "");

  /** This function can be used, if more than 10 graph styles are needed. */
  void setGraphStyles(CVR<std::string> styles);

  /** Sets up the grid. */
  void setGrid(bool x = true, bool y = true, bool x2 = false, bool y2 = false, bool z = false);

  /** Sets (or unsets) the axis or axes given in "axes" to be logarithmically scaled using the 
  given base. The "axes" can be any combination of x, x2, y, y2, z, cb, and r. If an empty string is 
  passed, all axes will be affected except r. */
  void setLogScale(std::string axes, double base = 10.0, bool shouldBeLogarithmic = true);
  // \todo find out why the base is important (does this determine the tics?). If not, it seems
  // unnecessary and might be removed - or at least, placed at the end of the argument list.

  /** Sets the plot range for x-, y- and z-axis. If you don't call this function, GNUPlot will 
  automatically choose an appropriate range. This is also true, if you pass invalid values, i.e. 
  values where the minium is greater or equal to the respective maximum. */
  void setRange(double xMin, double xMax, double yMin = 0.0, double yMax = 0.0, double zMin = 0.0,
    double zMax = 0.0);

  /** Sets the width and height of the plot in pixels. */
  void setPixelSize(unsigned int width, unsigned int height);

  /** Sets a title to appear above the plot. */
  void setTitle(std::string title);

  /** Adds an annotation to a 2D plot at the position given by the x- and y-coordinates. It 
  produces a "set label at (x,y)" command in the commandfile. The optional options can be given by
  a string of formatting options, for example "center" for centered text (see the set label 
  documentation in the GNUPlot manual). */
  void addAnnotation(double x, double y, CSR text, CSR options = "");


  /** \name Data setup */

  /** Sets the number of decimal digits after the dot with which the data will be written into the 
  datafile. This will have to be called before any of the setData functions is called, otherwise 
  it won't affect the respective call of setData(). Higher values will make the file larger. */
  void setDataPrecision(unsigned int numDecimals);

  /** This is a very general function for adding data from a nested vector. The outermost index
  runs over the blocks, the middle index runs over the lines in each block and the innermost index 
  over the columns in each line. It may be used, for example, for adding general 3-dimensional 
  datasets in which case the outermost index could run over a 1st parameter, the middle index over 
  a 2nd parameter and the innermost index runs from 0...2 (over the dimensions of the 3D vector).
  */
  template <class T>
  void addDataBlockLineColumn(const std::vector<std::vector<std::vector<T>>>& d);

  // for legacy - remove soon
  template <class T>
  void addData(const std::vector<std::vector<std::vector<T>>>& d)
  { addDataBlockLineColumn(d); }

  /** Similar to addDataBlockLineColumn but here the middle index runs over the columns and the 
  last index over the lines. In each block, all vectors for the same middle index must have the 
  same length. For example size(d[0][0]) == size(d[0][1]) == ... size[0][n] and also
  size(d[1][0]) == size(d[1][1]) == ... size(d[1][n]) etc. where n = size(d)-1. It can be used to
  add a function family with an arbitrary number of functions (or sets of such function 
  families). */
  template <class T>
  void addDataBlockColumnLine(const std::vector<std::vector<std::vector<T>>>& d);

  /** Adds a segmented dataset to our datafile. It's similar to 
  addData(int numRows, int numColumns, T **data) just that instead of passing a number of rows,
  you pass a number of blocks and an array of the blocklengths (the sum of the blocklengths gives 
  the number of rows). Each block (of rows) will be separated from the next one by a single blank
  line. */
  template <class T>
  void addData(int numBlocks, int *blockLengths, int numColumns, T **data);

  /** Adds the dataset given in the matrix "data" to our datafile. The 1st index is the column, the
  2nd index is the row in our datafile. This means, each row of the datafile will represent one 
  datapoint and the number of (space separated) values on that line is the dimensionality of that 
  datapoint. This means, the data[0] could represent an x-axis, data[1] could represent function 
  values y = f(x), data[2] could represent values for a second function y = g(x), etc. */
  template <class T>
  void addData(int numRows, int numColumns, T **data);

  /** Adds a dataset from a vector of complex numbers. The real parts will be the 1st column, the
  imaginary parts the 2nd column in the datafile. */
  template <class T>
  void addDataComplex(const std::vector<std::complex<T>>& d);

  /** Adds M+1 arrays of length N. The length-N x-array is the first column followed by M columns
  of y-arrays (all of length N) - so the y-matrix must be of the type: y[M][N]. */
  template <class T>
  void addDataArrays(int N, T *x, int M, T **y);

  /** Function to conveniently add a dataset from a bunch of arrays of the same length. For 
  example, the 1st array could contain the x-axis values and subsequent arrays could be values of
  a bunch of functions of x, like: c0=x, c1=f(x), c2=g(x), etc. */
  template <class T>
  void addDataArrays(int N, T *c0, T *c1 = nullptr, T *c2 = nullptr, T *c3 = nullptr, 
    T *c4 = nullptr, T *c5 = nullptr, T *c6 = nullptr, T *c7 = nullptr, T *c8 = nullptr, 
    T *c9 = nullptr);

  /** Adds data from univariate functions. */
  template <class T>
  void addDataFunctions(int N, T *x, T (*f0)(T), T (*f1)(T) = nullptr, T (*f2)(T) = nullptr,
    T (*f3)(T) = nullptr, T (*f4)(T) = nullptr, T (*f5)(T) = nullptr, T (*f6)(T) = nullptr, 
    T (*f7)(T) = nullptr, T (*f8)(T) = nullptr, T (*f9)(T) = nullptr);

  /** Like above, but you don't need to pass x-axis values. Instead, you pass a minimum and maximum
  value for the x-axis from which the axis values will be created internally */
  template <class T>
  void addDataFunctions(int N, T xMin, T xMax, T (*f0)(T), T (*f1)(T) = nullptr, 
    T (*f2)(T) = nullptr, T (*f3)(T) = nullptr, T (*f4)(T) = nullptr, T (*f5)(T) = nullptr, 
    T (*f6)(T) = nullptr, T (*f7)(T) = nullptr, T (*f8)(T) = nullptr, T (*f9)(T) = nullptr);

  /** Adds a dataset to our datafile that may represent a surface above the xy-plane, given - for 
  example - by a two-dimensional function z = f(x,y). You pass the arrays of x- and y values 
  (which should be of length Nx and Ny respectively) and a matrix of z-values such that 
  z[i][j] = f(x[i], y[j]). 
  Let N=Nx-1, M=Ny-1. Then the dataformat in the file is given by:
  x[0] y[0] z[0][0]    1st block
  x[0] y[1] z[0][1]
  ...  ...    ...
  x[0] y[M] z[0][M]

  x[1] y[0] z[1][0]    2nd block
  x[1] y[1] z[1][1]
  ...  ...    ...
  x[1] y[M] z[1][M]
  ...
  ...                  more blocks
  ...
  x[N] y[0] z[N][0]    N-th block
  x[N] y[1] z[N][1]
  ...  ...    ...
  x[N] y[M] z[N][M]   */
  template <class T>
  void addDataGrid(int Nx, int Ny, T *x, T *y, T **z);
    // maybe make it possible to take multiple z-pointers that produce a 4th, 5th etc. column
    // such additional columns may be needed for plotting direction-fields - but maybe not, let's
    // wait, if it is needed somewhere - if so, do it

  /** Adds data in the matrix-format to our datafile. This format is useful for more economical 
  storage of 3D data that is sampled at a grid of x-/y-values.
  If N=Nx-1, M=Ny-1, then the data format in the file is given by:
  Nx    x[0]    x[1]    ... x[N]
  y[0]  z[0][0] z[1][0] ... z[N][0]
  y[1]  z[0][1] z[1][1] ... z[N][1]
  ...     ...     ...   ...   ...
  y[M]  z[0][M] z[1][M] ... z[N][M] 
  the top-left Nx value is actually ignored by GNUPlot - it would be needed only for a binary 
  matrix format. */
  template <class T>
  void addDataMatrix(int Nx, int Ny, T *x, T *y, T **z);


  template <class T>
  void addDataBivariateFunction(int Nx, int Ny, T *x, T *y,  T (*f)(T,T));

  template <class T>
  void addDataBivariateFunction(int Nx, T xMin, T xMax, int Ny, T yMin, T yMax, T (*f)(T, T));
    // \todo: add 2 optional boolean parameters to let x- and/or y array be exponentially scaled 
    // instead of linear bool xLog = false, bool yLog = false


  /** Adds a graph to the plot. You must pass a string that describes how the data (from the 
  previously created datafile) should be used to produce one of the graphs in the plot. This is a 
  string that appears in the actual "plot" or "splot" command that will be created for GNUPlot 
  in the command file. For example, if you want to use column 1 and 2 of dataset 0 as x- and y-axis 
  for a function plot of a sine function using a black line and as legend "sin(x)", the string 
  could look like:  
  index 0 using 1:2 with lines lc rgb "#000000" title "sin(x)"
  If you don't call this function at all, this class will autogenerate these strings based on the 
  data you passed, thereby infering how you want to display the data. If you do call it, however, 
  you must call it for all graphs that you want to see in the plot (you can't mix auto-generation 
  and manual setup). Most of the times, you will probably rely on autogeneration and not have to 
  deal with this function, but it has been included (and exposed to client code) for greater 
  flexibility and generality in the use of this plotter class. */
  void addGraph(CSR descriptor);


  /** \name Inquiry */

  /** Returns the path of the datafile. */
  std::string getDataPath();


  /** \name Misc */

  /** Fills the array x of length N with equally spaced values from min...max. Useful for creating 
  an x-axis for a plot. */
  template <class T>
  static void rangeLinear(T *x, int N, T min, T max);

  /** Like rangeLinear, but the values are equally spaced on a logarithmic scale. Useful for 
  creating the array of values for a log scaled x-axis. */
  template <class T>
  static void rangeLogarithmic(T *x, int N, T min, T max);

  /** Executes GNUPlot with the appropriate commandline parameter to read the command file. */
  void invokeGNUPlot();




protected:

  /** Creates an initial empty file with the given path. */
  void initFile(const std::string &path);

  /** Converts the given string into a zero-terminated C-style string. The caller is responsible
  for deleting the returned pointer variable. */
  char* toZeroTerminatedString(const std::string &s);

  ///** Writes the passed data to the passed output stream object. The nested data vector is interpreted
  //as follows: the outermost vector runs over the blocks of data in the dataset (which are separated by 
  //single blank lines in the datafile), the middle vector runs over the data points (i.e. the number 
  //of lines of the current block) and the innermost vector runs over the dimensions of the datapoints
  //(i.e. the columns of one line in the datafile). After all blocks have been written, an additional 
  //blank line will be written (so we have two blank lines after the dataset, which indicates that the 
  //dataset is finished). */
  //template<class T>
  //void writeDataSet(const vector<vector<vector<T>>>& data, ostream& out);


  /** Calls the operating system to execute the command given by callString. */
  void systemCall(const std::string &callString);

  /** Adds the command for actually plotting the data to the commandfile. */
  void addPlotCommand(bool splot = false);

  /** Automaitcally creates the graph descriptors, in case the user didn't set them up manually. You 
  should tell it if GNUPlot is to be invoked with the regular 2D "plot" command or the 3D "splot" 
  command because the data will be interpreted differently in both cases. */
  void generateGraphDescriptors(bool splot);

  /** Returns a string to be used inside the plot command, that gives a graph with given index the 
  desired plotting style as defined in our member graphStyles and the desired color, as defined in 
  our graphColors member. */
  std::string getStyleString(unsigned int graphIndex);

  /** Returns the legend for the graph with given index. */
  std::string getGraphLegend(unsigned int graphIndex);



  // conversion of numbers to strings:
  std::string s(unsigned int x);   // conversion of unsigned integers for command file
  std::string s(double x);         // conversion of doubles for command file

  std::string sd(double x);        // conversion of doubles for data file
  // add one for floats
  std::string sd(int x);           // conversion of integers for data file


  // functions for the handling of variable argument lists

  template<class T> 
  T nullValue(T);

  std::string nullValue(std::string);

  template<class T>
  std::vector<T> collectLeadingNonNullArguments(T a0, T a1, T a2, T a3, T a4, T a5, T a6, T a7, T a8, 
    T a9);

  template<class T>
  void append(std::vector<T>& v, const std::vector<T>& appendix);

  template<class T>
  std::vector<std::vector<T>> wrapIntoVectors(int N, T *a0, T *a1, T *a2, T *a3, T *a4, T *a5, T *a6, 
    T *a7, T *a8, T *a9);

  void addToStringVector(std::vector<std::string>& v, CSR s0, CSR s1, CSR s2, CSR s3, CSR s4, 
    CSR s5, CSR s6, CSR s7, CSR s8, CSR s9);

  void setStringVector(std::vector<std::string>& v, CSR s0, CSR s1, CSR s2, CSR s3, CSR s4, CSR s5, 
    CSR s6, CSR s7, CSR s8, CSR s9);


  /** \name Data */

  // location for gnuplot and the directory for the temporary files:  
  std::string gnuplotPath;  // path, where gnuPlot is installed
  std::string dataPath;
  std::string commandPath;
  // make static

  // string arrays for styles and titles:
  std::vector<std::string> graphStyles;
  std::vector<std::string> graphTitles;


  std::vector<std::string> graphDescriptors;
    // an array of strings to be used in the plot command, like:
    // index 0 using 1:2 with lines lc rgb "#000000" title "sin(x)"
    // index 0 using 1:3 with lines lc rgb "#000000" title "cos(x)"
    // this vector can be set up from the user or - if left empty - will be created automatically 
    // by making some reasonable guess about how the user wants to plot the written data based on
    // out stored dataInfo array
    // maybe do: index 0 using 1:2, index 0 using 1:3, index 0 using 1:4, ...
    // then      index 1 using 1:2, index 1 using 1:3, index 1 using 1:4, ...
    //           index 2 ...
    // and so on, all with lines using colors and from our graphColors, graphStlyes arrays
    // hmm - but what if the style is one that uses several columns, have to think about it...
    //

  //char *formatString = "% 016.16le";
  char *formatString = nullptr;

  // for keeping track how much datasets we write and how many blocks and columns each dataset has:
  struct DataInfo
  {
    DataInfo(size_t numBlocks, size_t numColumns, std::string type = "")
    {
      this->numBlocks  = numBlocks;
      this->numColumns = numColumns;
      this->type       = type;
    }
    size_t numBlocks;
    size_t numColumns;
    std::string type;
  };
  std::vector<DataInfo> dataInfo; 


};

#endif
