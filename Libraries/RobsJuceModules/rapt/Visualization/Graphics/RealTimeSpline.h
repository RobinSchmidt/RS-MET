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
    QUADRATIC,      // not yet optimized - probably not so useful anyway
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

  /** When drawing splines, the resulting dots are not naturally equidistant which may result in
  artifacts like color discontinuities at the spline joints. This function switches some processing
  on that normalizes the dot-density along the spline segment - which means that the dots become
  equidistant. The computations for that are rather expensive, though - i think, much more 
  expensive that the spline drawing itself (but didn't measure yet). */
  void setNormalizeDotDensity(bool shouldNormalize) { normalizeDensity = shouldNormalize; }

  /** Sets up the density of the lines the connect our actual datapoints. It actually determines the
  number of artificial datapoints that are inserted (by linear interpolation) between our actual
  incoming datapoints. If set to zero, it will just draw the datapoints as dots. */
  void setDensity(TCor newDensity) { density = newDensity; }

  /** Sets a limit of the number of artificial datapoints that may be inserted per drawing
  operation. This is important to keep the CPU usage under control. */
  void setMaxNumDotsPerSegment(int newMaxNumDots) { maxNumDots = newMaxNumDots; }

  /** Sets the buffers to be used for the produced points (x- and y-coordinates as well as weights
  (i.e. color/brightness values). These are the buffers that will be filled when you call
  updateDotBuffers. */
  void setDotBuffers(TCor* bufX, TCor* bufY, TWgt* bufWeights, int bufLengths);

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** This is the main processing function to be called from client code. You supply a new pair of 
  input coordinates newX, newY and the function gives you back the locations and brightness values
  for the dots to draw. ... The return value is the number of dots produced. 
  The dot coordinates and weights will be written into the buffers that you have previously passed 
  via setDotBuffers.  */
  int updateDotBuffers(TCor newX, TCor newY);

  /** Adds a new incoming datapoint to the front of our point buffers and discards the last. Called 
  internally from updateDotBuffers, but there are situations where you may also need to call it 
  yourself.  */
  void updatePointBuffers(TCor newX, TCor newY);

  /** Shifts all the stored points in our point buffers by the specified amount. This can be 
  necessary in a scope in 1D mode when the ray jumps from the right border of the screen back
  to the left, i.e. wraps around. Then, all point locations must be shifted left by the screen 
  width. */
  void shiftPointBuffers(TCor shiftX, TCor shiftY);

  /** Resets the input point buffers to the given coodinate values. You will typically want to pass
  incoordinates that correspond to the origin of you coordinate system, either in pixel coordinates
  or in normalized device coordinates. */
  void reset(TCor x = TCor(0), TCor y = TCor(0));
  // todo: facilitate use with pixel-coordinates (as required by rsPhaseScopeBuffer) or normalized 
  // device coordinates (as required for use with OpenGL)


protected:


  /** Computes the number of dots to be produced for the current segment of given length. The 
  length of the segment depends on the interpolation mode (straight lines produce shorter segments 
  than spline arcs), so it must be passed as parameter. */
  int numDotsForSegment(TCor segmentLength);

  /** Computes the target weights/brightness values at the start and end of the current segment to 
  be used for the linear color gradient along the segment. */
  void getStartAndEndWeights(int numDots, TWgt* wStart, TWgt* wEnd);

  // the actual workhorse functions to update the dot-buffers using various interpolation modes:
  int dotsLinear();
  int dotsCubicHermite();
  int dotsQuadratic();
  // , dotsCubic2D, ...

  /** After the polynomial coefficient arrays a, b have been computed, this function fills the 
  array for the parameter t with values at which we evaluate the spline segment given by the 
  parametric equations:
  x(t) = a0 + a1*t + a2*t^2 + a3*t
  y(t) = b0 + b1*t + b2*t^2 + b3*t  
  The return value is the number of evaluation points/samples. */
  int prepareParameterArray();

  /** Evaluates the two cubic polynomials x(t),y(t) for the coordinates at t-values determined by 
  our t-array... Called from dotsCubicHermite and dotsQuadratic (the quadratic spline can just the
  same function by setting a3=b3=0 in the cubic spline equation. */
  void dotsCubic(TWgt w1, TWgt w2, int numDots);


  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  int  drawMode = LINEAR;
  bool useGradient = true;       // use color gradient to seamlessly join line segments

  bool scaleWeightsByNumDots = true;
  bool normalizeDensity = false; // or maybe use an int for precision of normalization (samples in 
                                 // numeric integration routine) ...expensive feature

  int maxNumDots = -1;         // -1 is code for: no limit
  TCor density = TCor(1);      // dot-density: 1: full. 0: dots only at input points
  TWgt brightness;             // determines weight by which dots are added in
  TWgt insertFactor;           // factor by which are pixels "inserted" (applied at sampleRate)
  TCor x[4], y[4];             // input signal buffers
  TWgt wOld;                   // old line end weight (brightness/color)

  // buffers to write the produced dots into:
  TCor *dotsX, *dotsY;
  TWgt *dotsW;
  int  dotBufferLength;

  // polynomial coefficients for for x(t), y(t):
  TCor a[4], b[4];

  // buffers for the density compensation and spline arc-length computations:
  std::vector<TCor> r, s, t, u; 

};
// maybe rename to RealTimeSplineGenerator or Interpolator

#endif