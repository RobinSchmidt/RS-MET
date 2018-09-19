#ifndef RAPT_REALTIMESPLINE_H_INCLUDED
#define RAPT_REALTIMESPLINE_H_INCLUDED

/** A class to compute dot-locations and brightness-values for realtime scope visualization of
incoming sample points. The class computes the dot location by interpolating the incoming 
data points in various ways. In the simplest case, it just produces dots that sit along a line that
connects two successive input sample points, but it can do more sophisticated things, like 
connecting the incoming data-points with spline segments and it can also ensure that the dots along
the spline segments are equally spaced. */

template<class TCor, class TWgt>
class rsRealTimeSpline
{

public:

  /** Standard constructor. Initializes to line drawing mode. */
  rsRealTimeSpline();



  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  enum drawModes
  {
    LINEAR = 0,
    CUBIC_HERMITE,
    // BRESENHAM,
    // WU,
    NUM_DRAW_MODES
  };

  /** Sets the drawing mode as one of the value in enum drawModes. */
  void setDrawMode(int newMode) { drawMode = newMode; }

  /** Sets the overall brightness. This parameter, together with the sample rate, determines the
  weight by which new dot are added in. */
  void setBrightness(TWgt newBrightness) { brightness = newBrightness; }
  // maybe needs to update the "insertFactor"?

  /** When we draw line segments, we normally scale the brightness by the reciprocal of the length
  so as to always add the same total amount of brightness into the screen. However, this may lead
  to abrupt color changes on line segment joints - so with this function, you can switch into a
  mode where a length-wise color gradient is used instead of a sudden switch. */
  void setUseColorGradient(bool shouldUseGradient) { useGradient = shouldUseGradient; }



  void setNormalizeDotDensity(bool shouldNormalize) { normalizeDensity = shouldNormalize; }

  /** Sets up the density of the lines the connect our actual datapoints. It actually determines the
  number of artificial datapoints that are inserted (by linear interpolation) between our actual
  incoming datapoints. If set to zero, it will just draw the datapoints as dots. */
  void setDensity(TCor newDensity) { density = newDensity; }

  /** Sets a limit of the number of artificial datapoints that may be inserted per drawing
  operation. This is important to keep the cpu usage under control. */
  void setMaxNumDotsPerSegment(int newMaxNumDots) { maxNumDots = newMaxNumDots; }




  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** This is the main processing function to be called from client code. You supply a new pair of 
  input coordinates newX, newY and the function gives you back the locations and brightness values
  for the dots to draw. ... The return value is the number of dots produced. */
  int getDotsForInputPoint(TCor newX, TCor newY, TCor* dotsX, TCor* dotsY, TWgt *weights, 
    int xywLength);
  // maybe don't pass in the buffers each time, instead, have a function 
  // setDotBuffers(TCor *bufX, TCor *bufY, TWgt *bufWeights, int bufferLengths) that client code
  // calls once on intitialization

  /** Resets the input point buffers to the given coodinate values. You will typically want to pass
  incoordinates that correspond to the origin of you coordinate system, either in pixel coordinates
  or in normalized device coordinates. */
  void reset(TCor x = TCor(0), TCor y = TCor(0));
  // todo: facilitate use with pixel-coordinates (as required by rsPhaseScopeBuffer) or normalized 
  // device coordinates (as required for use with OpenGL)


protected:


  int dotsLinear(TCor* dotsX, TCor* dotsY, TWgt weights, int xywLength);

  int dotsCubicHermite(TCor* dotsX, TCor* dotsY, TWgt weights, int xywLength);

  // getDotsQuadratic, ...


  int  drawMode = LINEAR;
  bool useGradient = true;     // use color gradient to seamlessly join line segments

  bool normalizeDensity = false;  // or maybe use an int for precision of normalization (samples in 
                                  // numeric integration routine)

  int maxNumDots = -1;         // -1 is code for: no limit
  TCor density = TCor(1);      // dot-density: 1: full. 0: dots only at input points
  TWgt brightness;             // determines weight by which dots are added in
  TWgt insertFactor;           // factor by which are pixels "inserted" (applied at sampleRate)
  TCor x[4], y[4];             // input signal buffers
  TWgt cOld;                   // old line end color



};
// maybe rename to RealTimeSplineGenerator

#endif