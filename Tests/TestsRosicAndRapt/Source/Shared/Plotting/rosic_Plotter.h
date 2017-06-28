#ifndef rosic_Plotter_h
#define rosic_Plotter_h

//// rosic-indcludes:
//#include "../infrastructure/rosic_FileInputOutput.h"
//#include "../math/rosic_ElementaryFunctionsReal.h"
//#include "../math/rosic_Complex.h"

namespace rosic
{

  /**

  This class allows for plotting data directly from C++ code by invoking GnuPlot. For this to work,
  the following preliminaries must be met:
  -Gnuplot must be installed in the path defined by rosic::gnuplotPath
  -the directory defined in tmpDataDir must exist

  */

  // define the location for gnuplot and the directory for the temporary files:  
  static const char* gnuplotPath = "E:/SoftwareWindows/Development/GnuPlot/gp426win32/gnuplot/bin/wgnuplot";  
  static const char* tmpDataDir  = "D:/TmpData/";

  class Plotter
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    Plotter(); 

    /** Destructor. */
    ~Plotter();  

    //---------------------------------------------------------------------------------------------
    // parameter settings:


    //---------------------------------------------------------------------------------------------
    // inquiry:


    //---------------------------------------------------------------------------------------------
    // static functions:

    /** Plots the passed function for the specified range of x-values. The optional numPoints 
    parameter determines the number of datapoints that will be created. */
    static void plotFunction(double (*f) (double), double xMin, double xMax, int numPoints = 1000);

    /** Plots two functions in one plot. */
    static void plotTwoFunctions(double (*f1) (double), double (*f2) (double), 
      double xMin, double xMax, int numPoints = 1000);

    /** Plots a family of functions that is parametrized by some parameter p for 5 particular 
    values of the parameter. The function that defines the family must have a signature of:
    double f(double x, double p), where x denotes the primary function argument and p denotes the 
    parameter, like for example: double sn(double u, double m) - which could be the Jacobian 
    elliptic function sn with argument u and parameter m. */
    static void plotFunctionFamily(double (*f) (double, double), double xMin, double xMax, 
      double p1, double p2, double p3, double p4, double p5, int numPoints = 1000);

    /** Plots the passed y-arrays against the axis defined by the x-array. There should be at least
    one y-array, the others are optional. This function is intended for quick and dirty plots 
    without actually instantiating a Plotter object. If you want to customize your plot (ranges, 
    colors, etc.), you should instantiate an object and use the various set-functions to control
    the plot options. */
    static void plotData(int numValues, double *x, double *y1, double *y2 = NULL,
      double *y3 = NULL, double *y4 = NULL, double *y5 = NULL);


    static void plotComplexMagnitudes(int numValues, double *x, bool inDecibels, Complex *y1, 
      Complex *y2 = NULL, Complex *y3 = NULL, Complex *y4 = NULL, Complex *y5 = NULL);


    //static void plotImage(double *data, int pixelWidth, int pixelHeight, double xMin, double xMax,
    //  double yMin, double yMax);

    static void plotImage(double *data, int pixelWidth, int pixelHeight, double *xAxis, 
      double *yAxis);


    /** Plots the magnitude response of an analog filter with zeros, poles and gain given in "z", "p", "k". The order of the filter should
    be passed in "N". The "wl", "wu" parameters are the normalized radian frequencies between which the magnitude response is plotted and
    the "resolution" parameter gives the number of these frequecies. */
    static void plotAnalogMagnitudeResponse(Complex *z, Complex *p, double k, int N, double wl, double wu, int resolution = 200);

    /** Plots the phase response of an analog filter @see plotAnalogMagnitudeResponse. */
    static void plotAnalogPhaseResponse(Complex *z, Complex *p, double k, int N, double wl, double wu, int resolution = 200);



    //=============================================================================================

  protected:

  };

} // end namespace rosic

#endif // rosic_Plotter_h
