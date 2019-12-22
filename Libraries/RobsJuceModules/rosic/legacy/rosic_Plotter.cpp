#include "rosic_Plotter.h"
using namespace rosic;
//-------------------------------------------------------------------------------------------------
// construction/destruction:

Plotter::Plotter()
{

}

Plotter::~Plotter()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:


//-------------------------------------------------------------------------------------------------
// others:

void Plotter::plotTwoFunctions(double (*f1) (double), double (*f2) (double), 
                               double xMin, double xMax, int numPoints)
{
  double *x  = new double[numPoints];
  double *y1 = new double[numPoints];
  double *y2 = new double[numPoints];
  double tmp;
  for(int n=0; n<numPoints; n++)
  {
    tmp   = (double) n / (double) (numPoints-1);
    tmp   = RAPT::rsLinToLin(tmp, 0.0, 1.0, xMin, xMax);
    x[n]  = tmp;
    y1[n] = f1(tmp);
    y2[n] = f2(tmp);
  }
  plotData(numPoints, x, y1, y2);
  delete[] x;
  delete[] y1;
  delete[] y2;
}

void Plotter::plotFunctionFamily(double (*f) (double, double), double xMin, double xMax, double p1, 
                                 double p2, double p3, double p4, double p5, int numPoints)
{
  double *x  = new double[numPoints];
  double *y1 = new double[numPoints];
  double *y2 = new double[numPoints];
  double *y3 = new double[numPoints];
  double *y4 = new double[numPoints];
  double *y5 = new double[numPoints];
  double tmp;
  for(int n=0; n<numPoints; n++)
  {
    tmp   = (double) n / (double) (numPoints-1);
    tmp   = RAPT::rsLinToLin(tmp, 0.0, 1.0, xMin, xMax);
    x[n]  = tmp;
    y1[n] = f(tmp, p1);
    y2[n] = f(tmp, p2);
    y3[n] = f(tmp, p3);
    y4[n] = f(tmp, p4);
    y5[n] = f(tmp, p5);
  }
  plotData(numPoints, x, y1, y2, y3, y4, y5);
  delete[] x;
  delete[] y1;
  delete[] y2;
  delete[] y3;
  delete[] y4;
  delete[] y5;
}

void Plotter::plotFunction(double (*f) (double), double xMin, double xMax, int numPoints)
{
  double *x = new double[numPoints];
  double *y = new double[numPoints];
  double tmp;
  for(int n=0; n<numPoints; n++)
  {
    tmp  = (double) n / (double) (numPoints-1);
    tmp  = RAPT::rsLinToLin(tmp, 0.0, 1.0, xMin, xMax);
    x[n] = tmp;
    y[n] = f(tmp);
  }
  plotData(numPoints, x, y);
  delete[] x;
  delete[] y;
}

void Plotter::plotData(int numValues, double *x, double *y1, double *y2, double *y3, double *y4, 
                       double *y5)
{
  // store the data in a temporary file:
  char* tmpDataFile  = "tmpData.dat";
  char* dataPath = new char[40];
  strcpy(dataPath, tmpDataDir);
  strcat(dataPath, tmpDataFile);
  writeDataToFile(dataPath, numValues, x, y1, y2, y3, y4, y5);

  // create the string that represents the command-line invocation of gnuplot:
  char* commandFile = "gnuplotCommands.txt"; 
  char* commandPath = new char[40];
  strcpy(commandPath, tmpDataDir);
  strcat(commandPath, commandFile);

  char* callString = new char[200];
  strcpy(callString, "\"");
  strcat(callString, gnuplotPath);
  strcat(callString, "\" ");
  strcat(callString, " ");
  strcat(callString, commandPath);
  strcat(callString, " -");        // the minus avoids gnuplot to immediately exit 

  //---------------------------------------------------------------------------
  // create the file that contains the commands to be batch-processed by gnuplot:
  char* commands  = new char[1024];
  char* tmpString = new char[1024];

  // create commands for setting up the grid:
  strcpy(commands, "");
  strcat(commands, "set grid\n");

  // create commands for setting up the range for the y-axis:
  double yMax    = RAPT::rsArrayTools::maxValue(y1, numValues);
  double yMin    = RAPT::rsArrayTools::minValue(y1, numValues);
  if( y2 != NULL ) 
  {
    yMax = RAPT::rsMax( yMax, RAPT::rsArrayTools::maxValue(y2, numValues) );
    yMin = RAPT::rsMin( yMin, RAPT::rsArrayTools::minValue(y2, numValues) );
  }
  if( y3 != NULL ) 
  {
    yMax = RAPT::rsMax( yMax, RAPT::rsArrayTools::maxValue(y3, numValues) );
    yMin = RAPT::rsMin( yMin, RAPT::rsArrayTools::minValue(y3, numValues) );
  }
  if( y4 != NULL ) 
  {
    yMax = RAPT::rsMax( yMax, RAPT::rsArrayTools::maxValue(y4, numValues) );
    yMin = RAPT::rsMin( yMin, RAPT::rsArrayTools::minValue(y4, numValues) );
  }
  if( y5 != NULL ) 
  {
    yMax = RAPT::rsMax( yMax, RAPT::rsArrayTools::maxValue(y5, numValues) );
    yMin = RAPT::rsMin( yMin, RAPT::rsArrayTools::minValue(y5, numValues) );
  }
  double yRange = yMax-yMin;
  if( yRange < 1.e-8 )
    yRange = 1.0;
  double yMargin = 0.05*yRange;
  yMax += yMargin;
  yMin -= yMargin;
  sprintf(tmpString, "%s %.8f %s %.8f %s", 
    "set yrange [", yMin, ":", yMax, "]\n");
  strcat(commands, tmpString);

  // create commands for plotting the data:
  strcat(commands, "plot '");
  strcat(commands, dataPath);
  strcat(commands, "' index 0 lc rgbcolor \"#000000\" with lines notitle");
  if( y2 != NULL ) 
  {
    strcat(commands, ", '");
    strcat(commands, dataPath);
    strcat(commands, "' index 1 linecolor rgbcolor \"#0000FF\" with lines notitle");
  }
  if( y3 != NULL ) 
  {
    strcat(commands, ", '");
    strcat(commands, dataPath);
    strcat(commands, "' index 2 linecolor rgbcolor \"#007700\" with lines notitle");
  }
  if( y4 != NULL ) 
  {
    strcat(commands, ", '");
    strcat(commands, dataPath);
    strcat(commands, "' index 3 linecolor rgbcolor \"#BB0000\" with lines notitle");
  }
  if( y5 != NULL ) 
  {
    strcat(commands, ", '");
    strcat(commands, dataPath);
    strcat(commands, "' index 4 linecolor rgbcolor \"#9900FF\" with lines notitle");
  }
  strcat(commands, "\n");

  // write the string with the commands into a file to be batch processed by 
  // gnuplot:
  writeStringToFile(commandPath, commands);

  //---------------------------------------------------------------------------
  // batch command file created - now we are ready to actually invoke gnuplot:
  system(callString);

  // clean up memory:
  delete[] dataPath;
  delete[] commandPath;
  delete[] callString;
  delete[] commands;
  delete[] tmpString;
}

void Plotter::plotComplexMagnitudes(int numValues, double *x, bool inDecibels, 
                                    Complex *y1, Complex *y2, Complex *y3, 
                                    Complex *y4, Complex *y5)
{
  double *ym1 = new double[numValues];
  for(int i=0; i<numValues; i++)
  {
    if( inDecibels == true )
      ym1[i] = RAPT::rsAmpToDb(y1[i].getRadius());
    else
      ym1[i] = y1[i].getRadius();
  }

  double *ym2 = NULL; 
  double *ym3 = NULL; 
  double *ym4 = NULL; 
  double *ym5 = NULL; 
  if( y2 != NULL )
  {
    ym2 = new double[numValues];
    for(int i=0; i<numValues; i++)
    {
      if( inDecibels == true )
        ym2[i] = RAPT::rsAmpToDb(y2[i].getRadius());
      else
        ym2[i] = y2[i].getRadius();
    }
  }
  if( y3 != NULL )
  {
    ym3 = new double[numValues];
    for(int i=0; i<numValues; i++)
    {
      if( inDecibels == true )
        ym3[i] = RAPT::rsAmpToDb(y3[i].getRadius());
      else
        ym3[i] = y3[i].getRadius();
    }
  }
  if( y4 != NULL )
  {
    ym4 = new double[numValues];
    for(int i=0; i<numValues; i++)
    {
      if( inDecibels == true )
        ym4[i] = RAPT::rsAmpToDb(y4[i].getRadius());
      else
        ym4[i] = y4[i].getRadius();
    }
  }
  if( y5 != NULL )
  {
    ym5 = new double[numValues];
    for(int i=0; i<numValues; i++)
    {
      if( inDecibels == true )
        ym5[i] = RAPT::rsAmpToDb(y5[i].getRadius());
      else
        ym5[i] = y5[i].getRadius();
    }
  }

  plotData(numValues, x, ym1, ym2, ym3, ym4, ym5); 

  delete[] ym1;
  delete[] ym2;
  delete[] ym3;
  delete[] ym4;
  delete[] ym5;
}

/*
void Plotter::plotImage(double *data, int pixelWidth, int pixelHeight, double xMin, double xMax,
                        double yMin, double yMax)
{
  // store the data in a temporary file:
  char* tmpDataFile = "tmpData.dat";
  char* dataPath    = new char[40];
  strcpy(dataPath, tmpDataDir);
  strcat(dataPath, tmpDataFile);
  writeDataToFile(dataPath, data, pixelWidth, pixelHeight);

  // create the string that represents the command-line invocation of gnuplot:
  char* commandFile = "gnuplotCommands.txt"; 
  char* commandPath = new char[40];
  strcpy(commandPath, tmpDataDir);
  strcat(commandPath, commandFile);

  char* callString = new char[200];
  strcpy(callString, gnuplotPath);
  strcat(callString, " ");
  strcat(callString, commandPath);
  strcat(callString, " -");        // the minus avoids gnuplot to immediately exit 

  //---------------------------------------------------------------------------
  // create the file that contains the commands to be batch-processed by gnuplot:
  char* commands  = new char[1024];
  char* tmpString = new char[1024];

  // create commands for setting up the axes:
  strcpy(commands, "");



  double xRange  = xMax-xMin;
  double xMargin = 0.05*xRange;
  xMax += xMargin;
  xMin -= xMargin;
  sprintf(tmpString, "%s %d %s %d %s", "set xrange [", 0, ":", pixelWidth, "]\n");
  //sprintf(tmpString, "%s %.8f %s %.8f %s", "set xrange [", xMin, ":", xMax, "]\n");
  strcat(commands, tmpString);

  double yRange  = yMax-yMin;
  double yMargin = 0.05*yRange;
  yMax += yMargin;
  yMin -= yMargin;
  sprintf(tmpString, "%s %d %s %d %s", "set yrange [", 0, ":", pixelHeight, "]\n");
  //sprintf(tmpString, "%s %.8f %s %.8f %s", "set xrange [", xMin, ":", xMax, "]\n");
  strcat(commands, tmpString);


  // create commands for plotting the data:
  strcat(commands, "plot '");
  strcat(commands, dataPath);
  //strcat(commands, "' binary array=16x16 flipy format='%uchar' with rgbimage");
  strcat(commands, "' binary array=");
  sprintf(tmpString, "%d", pixelWidth);
  strcat(commands, tmpString);
  strcat(commands, "x");
  sprintf(tmpString, "%d", pixelHeight);
  strcat(commands, tmpString);
  //strcat(commands, " flipy format='%uchar' with rgbimage \n");
  strcat(commands, " format='%uchar' with rgbimage \n");

  // plot 'blutux.rgb' binary array=128x128 flipy format='%uchar' with rgbimage

  // write the string with the commands into a file to be batch processed by 
  // gnuplot:
  writeStringToFile(commandPath, commands);

  //---------------------------------------------------------------------------
  // batch command file created - now we are ready to actually invoke gnuplot:
  system(callString);

  // clean up memory:
  delete[] dataPath;
  delete[] commandPath;
  delete[] callString;
  delete[] commands;
  delete[] tmpString;
}
*/



void Plotter::plotImage(double *data, int pixelWidth, int pixelHeight, double *xAxis, 
                        double *yAxis)
{
  // store the data in a temporary file:
  char* tmpDataFile = "tmpData.dat";
  char* dataPath    = new char[40];
  strcpy(dataPath, tmpDataDir);
  strcat(dataPath, tmpDataFile);
  writeDataToFile(dataPath, data, pixelWidth, pixelHeight, xAxis, yAxis);

  // create the string that represents the command-line invocation of gnuplot:
  char* commandFile = "gnuplotCommands.txt"; 
  char* commandPath = new char[40];
  strcpy(commandPath, tmpDataDir);
  strcat(commandPath, commandFile);

  char* callString = new char[200];
  strcpy(callString, gnuplotPath);
  strcat(callString, " ");
  strcat(callString, commandPath);
  strcat(callString, " -");        // the minus avoids gnuplot to immediately exit 

  //---------------------------------------------------------------------------
  // create the file that contains the commands to be batch-processed by gnuplot:
  char* commands  = new char[1024];
  char* tmpString = new char[1024];

  // create commands for setting up the axes:
  strcpy(commands, "");


  double xMax    = RAPT::rsArrayTools::maxValue(xAxis, pixelWidth);
  double xMin    = RAPT::rsArrayTools::minValue(xAxis, pixelWidth);
  double xRange  = xMax-xMin;
  double xMargin = 0.05*xRange;
  xMax += xMargin;
  xMin -= xMargin;
  //sprintf(tmpString, "%s %d %s %d %s", "set xrange [", 0, ":", pixelWidth-1, "]\n");
  sprintf(tmpString, "%s %.8f %s %.8f %s", "set xrange [", xMin, ":", xMax, "]\n");
  strcat(commands, tmpString);

  double yMax    = RAPT::rsArrayTools::maxValue(yAxis, pixelHeight);
  double yMin    = RAPT::rsArrayTools::minValue(yAxis, pixelHeight);
  double yRange  = yMax-yMin;
  double yMargin = 0.05*yRange;
  yMax += yMargin;
  yMin -= yMargin;
  //sprintf(tmpString, "%s %d %s %d %s", "set yrange [", 0, ":", pixelHeight-1, "]\n");
  sprintf(tmpString, "%s %.8f %s %.8f %s", "set yrange [", yMin, ":", yMax, "]\n");
  strcat(commands, tmpString);


  double zMax    = RAPT::rsArrayTools::maxValue(data, pixelWidth*pixelHeight);
  double zMin    = RAPT::rsArrayTools::minValue(data, pixelWidth*pixelHeight);
  double zRange  = zMax-zMin;
  double zMargin = 0.05*zRange;
  zMax += zMargin;
  zMin -= zMargin;
  //sprintf(tmpString, "%s %d %s %d %s", "set yrange [", 0, ":", pixelHeight-1, "]\n");
  sprintf(tmpString, "%s %.8f %s %.8f %s", "set zrange [", zMin, ":", zMax, "]\n");
  strcat(commands, tmpString);

  // create commands for plotting the data:
  //strcat(commands, "plot '");
  //strcat(commands, dataPath);
  //strcat(commands, "' binary with rgbimage \n");

  strcat(commands, "plot '");
  strcat(commands, dataPath);
  strcat(commands, "' binary with image \n");


  // plot 'blutux.rgb' binary array=128x128 flipy format='%uchar' with rgbimage

  // write the string with the commands into a file to be batch processed by 
  // gnuplot:
  writeStringToFile(commandPath, commands);

  //---------------------------------------------------------------------------
  // batch command file created - now we are ready to actually invoke gnuplot:
  system(callString);

  // clean up memory:
  delete[] dataPath;
  delete[] commandPath;
  delete[] callString;
  delete[] commands;
  delete[] tmpString;
}

void Plotter::plotAnalogMagnitudeResponse(Complex *z, Complex *p, double k, int N, double wl, double wu, int resolution)
{
  double *w = new double[resolution];
  double *m = new double[resolution];
  RAPT::rsArrayTools::fillWithRangeLinear(w, resolution, wl, wu);
  rosic::rsFilterAnalyzerD::getAnalogMagnitudeResponse(
    rsCastPointer(z), rsCastPointer(p), k, N, w, m, resolution);
  Plotter::plotData(resolution, w, m);
  delete[] w;
  delete[] m;
}

void Plotter::plotAnalogPhaseResponse(Complex *z, Complex *p, double k, int N, double wl, double wu, int resolution)
{
  double *w = new double[resolution];
  double *m = new double[resolution];
  RAPT::rsArrayTools::fillWithRangeLinear(w, resolution, wl, wu);
  rosic::rsFilterAnalyzerD::getAnalogPhaseResponse(
    rsCastPointer(z), rsCastPointer(p), k, N, w, m, resolution);
  Plotter::plotData(resolution, w, m);
  delete[] w;
  delete[] m;
}
