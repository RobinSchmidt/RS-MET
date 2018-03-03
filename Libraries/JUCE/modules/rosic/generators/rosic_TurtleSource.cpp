TurtleSource::TurtleSource()
{
  // init to square:
  turnAngle = 90;
  turtle.setAngle(turnAngle);
  setTurtleCommands("F+F+F+F+");
  updateWaveTable();

  // experimental:
  turtleLowpass.setSampleRate(1.1); // later try 1.0
  turtleLowpass.setFrequency(0.5); // check out, if 0.5 leads to bypass coeffs
  //turtleLowpass.setApproximationMethod(RAPT::rsPrototypeDesigner<double>::BESSEL);
  turtleLowpass.setApproximationMethod(RAPT::rsPrototypeDesigner<double>::ELLIPTIC);
  turtleLowpass.setPrototypeOrder(4);
  turtleLowpass.setMode(RAPT::rsInfiniteImpulseResponseDesigner<double>::LOWPASS);
}

void TurtleSource::setTurtleCommands(const std::string& commands)
{
  turtleCommands = commands;
  numLines = turtle.getNumberOfLines(turtleCommands); 
  updateMeanAndNormalizer();
  updateIncrement();
  reset();
}

void TurtleSource::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  incUpToDate = false;
}

void TurtleSource::setFrequency(double newFrequency)
{
  frequency = newFrequency;
  incUpToDate = false;
}

void TurtleSource::setRotation(double newRotation)
{
  rotator.setAngle((PI/180) * newRotation);
}

void TurtleSource::setResetAfterCycles(int numCycles) 
{ 
  cyclicReset = numCycles; 
}

void TurtleSource::setResetAfterLines(double numLines)  // doule take a double
{
  lineCountReset      = numLines;
  lineCountResetFloor = (int) lineCountReset;
  lineCountResetFrac  = lineCountReset - lineCountResetFloor; 

  // maybe factor out:
  if(lineCountResetErr > 0.5) {
    lineCountResetAlt  = lineCountResetFloor+1;
    lineCountResetErr -= 1.0; }
  else
    lineCountResetAlt  = lineCountResetFloor;

  incUpToDate = false;
}

void TurtleSource::setAngle(double newAngle) 
{
  turnAngle = newAngle;
  turtle.setAngle(turnAngle);
  updateMeanAndNormalizer();
  tableUpToDate = false; 
}

int TurtleSource::getNumTurtleLines()
{
  return numLines;
}

void TurtleSource::updateWaveTable()
{
  TurtleGraphics tmpTurtle;
  tmpTurtle.setAngle(turnAngle);
  tmpTurtle.translate(turtleCommands, tableX, tableY);
  tableUpToDate = true;
}

void TurtleSource::reset()
{
  pos = 0;
  cycleCount = 0;
  resetTurtle();
  updateXY();      // rename to drawNextLineToBuffer
  lineIndex = 0;
}

void TurtleSource::resetTurtle()
{
  commandIndex = 0;
  turtle.init(0, 0, 1, 0);
}

void TurtleSource::goToLineSegment(int targetLineIndex)
{
  if(numLines == 0) { x[0] = y[0] = x[1] = y[1] = 0; return; }
  if(useTable)
  {
    lineIndex = targetLineIndex;
    x[0] = tableX[lineIndex];
    y[0] = tableY[lineIndex];
    x[1] = tableX[lineIndex+1];
    y[1] = tableY[lineIndex+1];
  }
  else
  {
    while(lineIndex != targetLineIndex)
      goToNextLineSegment(); // increments lineIndex with wrap around
  }

  // later: update interpolator coeffs here according to x,y buffers
}

void TurtleSource::goToNextLineSegment()
{
  updateXY();

  if(lineCountReset != 0) 
  {
    lineCount++;

    if(lineCount >= lineCountResetAlt) { // >= not ==, bcs we may get beyond when user adjusts it
      lineCount = 0;
      resetTurtle();
      updateXY(); 

      lineCountResetErr += lineCountResetFrac;

      // maybe factor out:
      if(lineCountResetErr > 0.5) {
        lineCountResetAlt  = lineCountResetFloor+1;
        lineCountResetErr -= 1.0;
      }
      else
        lineCountResetAlt  = lineCountResetFloor;
    }
  }



  lineIndex++;
  if(lineIndex == numLines) {
    lineIndex = 0;

    cycleCount++;
    if(cyclicReset == 0)
      cycleCount = 0;
    else if(cycleCount >= cyclicReset) {
      cycleCount = 0;
      resetTurtle();
      updateXY();
    }


  }
}

// make reset of the turtle optional, or reset after a certain number of cycles, 0: never
// reset. the turtle is free running in this case. when the angle is slightly off the perfect
// value, the picture rotates - but we need to make sure to define axioms in a way that they
// head into the same direction as in the beginning after completing the cycle, for example
// "F--F--F--" instead of "F--F--F" for the initial triangle for the koch snowflake
// maybe have a soft-reset parameter that resets only partially (i.e. interpolates between
// current and initial state )

void TurtleSource::updateXY()  // // rename to drawNextLineToBuffer or updateLineBuffer
{
  if(turtleCommands.size() == 0)
    return;
  bool xyUpdated = false;
  while(xyUpdated == false) {
    bool draw = turtle.interpretCharacter(turtleCommands[commandIndex]);
    commandIndex++;
    if(commandIndex == turtleCommands.size())
      commandIndex = 0;
    if(draw) {

      if(!antiAlias)
      {
        x[0] = turtle.getStartX();
        y[0] = turtle.getStartY();
        x[1] = turtle.getEndX();
        y[1] = turtle.getEndY();
      }
      else
      {
        //// experimental anti aliasing:
        //double k = 0.9*inc / (numLines+1); // or maybe make the factor 0.9 a parameter
        //k = 0;
        //x[0] = x[1];
        //y[0] = y[1];
        //x[1] = (1-k)*turtle.getEndX() + k*x[0];
        //y[1] = (1-k)*turtle.getEndY() + k*y[0];

        x[0] = x[1];
        y[0] = y[1];
        x[1] = turtle.getEndX();
        y[1] = turtle.getEndY();
        turtleLowpass.getSampleFrameStereo(&x[1], &y[1]);
        // hmm - with FourArmAlien, it sounds worse with anti-aliasing - we need to create an 
        // experiment for this - also compare it to using oversampling
        // ...with HexaGrid, the anti aliasing improves the sound

        // ...but it deosn't seem to work well anyway ...maybe the turtle should apply a lowpass
        // internally - in each step, instead of setting x += dx, 
        // set x = lowpass.getSample(x+dx) ...the turtle could use a bessel lowpass ...make a class
        // FilteredTurtle...do things like turtle.setSmoothing/Lowpass
        // when inc>1, we are actually stepping over multiple trurtle-genertated line segments in 
        // each sample - so the turtle produces an oversampled output with respect to our sample 
        // rate ..instead of the simple averaging above, we could do:
        // x[1] = turtleLowpass.getSample(turtle.getEndX());
      }

      xyUpdated = true;
    }
  }
}

void TurtleSource::updateIncrement()
{
  // old:
  //inc = numLines * frequency / sampleRate; // avoid division

  // new:
  if(lineCountReset == 0)
    inc = numLines * frequency / sampleRate;
  else
    inc = rmin(lineCountReset, double(numLines)) * frequency / sampleRate;




  turtleLowpass.setSampleRate(rmin(1/inc, 1.1)); // use 1.0 later
  incUpToDate = true;
}

void TurtleSource::updateMeanAndNormalizer()
{
  // run through all turtle commands once, thereby keep track of the mean and min/max values of the
  // generated curve, then set our meanX, meanY, normalizer members accordingly

  double minX = 0, maxX = 0, minY = 0, maxY = 0, sumX = 0, sumY = 0;
  TurtleGraphics tmpTurtle;
  tmpTurtle.setAngle(turnAngle);
  tmpTurtle.init();
  int N = (int) turtleCommands.size();
  for(int i = 0; i < N; i++) {
    bool lineDrawn = tmpTurtle.interpretCharacter(turtleCommands[i]);
    if( lineDrawn == true )  {
      double x = tmpTurtle.getX();
      double y = tmpTurtle.getY();
      minX  = rmin(minX, x);
      maxX  = rmax(maxX, x);
      minY  = rmin(minY, y);
      maxY  = rmax(maxY, y);
      sumX += x;
      sumY += y;
    }
  }
  double tmp = 1.0 / (numLines+1); // maybe have a member for that (optimization)
  meanX = sumX * tmp;
  meanY = sumY * tmp;
  minX -= meanX; maxX -= meanX;
  minY -= meanY; maxY -= meanY;
  maxX  = rmax(fabs(minX), fabs(maxX));
  maxY  = rmax(fabs(minY), fabs(maxY));
  normalizer = 1.0 / rmax(maxX, maxY);

  // in free running mode, this so computed mean is sometimes (often) wrong because one cycle is 
  // only part (for example half) of a larger shape - in this case, we would need the mean over two 
  // cycles ...maybe we can somehow deduce the symmetry in order to compute a better mean? maybe we
  // should offer different normalization modes...maybe have a highpass and leveller running on it
  // in realtime (but not oversampled)
}

/*
Ideas:

-make TurningAngle modulatable
-optionally compensate for changes in turning-angle by a rotation of the whole image, i.e. rotate
 image by -TurningAngle ..or maybe -compensationAmount*TurningAngle
-requires to not pre-render a wavetable but only render the L-system output string and interpret 
 it on the fly (the string should be stripped from meaningless characters, retaining only F,+,- 
 before)
-may require to compute shift and scale (to remove DC and normalize) when a new angle is set up
-maybe these values should be precomputed for various angles, for example 0,5,10,15,.. and 
 interpolation be used (to avoid to expensively compute the actually desired value whenever the
 angle changes)
-allow a turtle syntax like +30 or -45 to mean: turn 30° left or 45° right (instead of whatever 
 value the turtle drawer is set to)
-make rotation angle quantizable to a set of numbers that the user can enter, for example
 120,90,72,60,45,36,30,24,22.5,20,15,12,10
-have different loop modes: forward, backward, alternating (avoids jumps for non-closed curves)
-maybe in backward passes, (optionally) interpret angles as their negative (necessarry to actually
 "go back")
-maybe have an additional string of turtle commands that can be applied after each cycle has passed
 (empty by default)...or maybe even different command strings to be applied after different numbers 
 of cycles ..this defines the "turn/wrap around behavior"
-interpolation modes: left, right, nearest, linear, cubic
-allow a traversal speed (or better: duration/length) to be associated with each point (i.e. the 
 line segment that goes toward that point)...like a normal 'F' would assume 1, an 'f' would 
 assume '0', but we could also have F2 to let it take 2 time units to traverse the segment or 
 1/3 to take one thrid of the time unit, the increment computation would then not use "numLines" 
 as multiplier but sum-of-traversal-duration
-optionally "upsample" the turtle commands like F+F-FF -> x3 -> FFF+FFF-FFFFFF, this givens a finer
 resulution of lines. by itself, it's useless but now we apply additional bending inside ths line
 (like turning by a few degrees after each F)
-randomize turn angles, user can set seed and distribution (uniform, bell, bimodal, etc.)
-randomize step lengths
-maybe let the turle have scaleX, scaleY members that are applied in each step - these can be 
 modulated and/or randomized
-have a fine turning angle delta (maybe 0 to 10%) with keytracking
-DC blocker and leveller may replace or complement normalization (it doesn't always work right, 
 especially in free running mode and turning angle modulation is prohibitively expensive with it)
 ...maybe we could handle it in ModulationTarget by defining an update-interval (by default, 1 
 sample but could also be something like 100 or 1000 samples - greatly reduces modulation costs)
-maybe, as (pseudo) anti-alias strategy, we should use a fixed internal sample-rate that is some 
 multiple of the signal frequency. technically, the signal would still produce alias frequencies, 
 but they would be harmonically related to the signal frequency and therefore much less annoying, 
 maybe even beneficial
 -maybe this is a good opportunity to write a resampler class for arbitrary ratios
-the turtle could interpret '*' as "multiply speed by factor" (which is a (realtime) parameter) and 
'/' as "divide speed by factor" (or maybe the speed should increase decreas linearly?)
-in free-running mode, a turn angle that is slightly off some ideal value leads to rotating image 
-maybe this can be simulated in wavetable mode, by modulating the increments for x,y in a particular 
 way - maybe modulate incX in any way and modulate incY according to the condition that 
 incX^2 + incY^2 = const = 2*inc^2 -> incY = sqrt(2*inc^2 - incX^2)...why?...because it 
 renormalizes the implied "turtle" step size...i think...maybe try it..or maybe not modulate it 
 but just set it to that value?
-provide soft-reset (interpolate between current turtle state and initial state)
-provide independent resets for turtle position and direction
-maybe the line count reset can be made continuous by alternating between floor and ceil of the 
 given number...or maybe by just comparing "pos" to lineCountReset in getSampleFrame?
 -maybe then it can be expressed also as fraction of the cycle-Length?
 -maybe then the cyclic reset can be continuous? too
 -maybe we can then have an arbitrary number of independent reset-counters?
-the reset intervals need keytracking because otherwise the speed of the "modulation" depends on
 the note (and is unpleasantly fast for higher notes)
-ideal would be a parameter that the user can set up in terms of the modulation speed 
 -that seems to depend on the difference between (a multiple of) numLines and lineCountReset 
  if lineCountReset is close to a multiple of numLines, modulation is slow - but the multiple
  can also be 1.5 ...but also, the higher the multiple, the slower and the higher the played
  note, the faster, so 
  speed = k * (numLines-resetInterval) * noteFreq?
  that k depends in some way on the ratio resetInterval/numLines - maybe gcd/lcm is involved?



*/