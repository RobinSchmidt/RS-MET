ResetCounter::ResetCounter()
{
  setParametersToOffMode();
  reset(); 
}

void ResetCounter::setInterval(double t)
{
  if(t >= maxInterval)
    setParametersToOffMode();
  else {
    interval = t; 
    intPart  = (int)t;
    fracPart = t - intPart; 
  }
  updateJitteringLimit();
}

bool ResetCounter::tick()
{
  if(intPart == maxInterval)
    return false; // counter is off 
  counter++;
  if(counter >= jitteringLimit) { // >= not ==, bcs we may get beyond when user adjusts it
    counter = 0;
    errAccu += fracPart;
    updateJitteringLimit();
    return true;
  }
  return false;
}

void ResetCounter::setParametersToOffMode()
{
  interval = INF;
  intPart  = maxInterval;
  fracPart = 0;
}

void ResetCounter::reset()
{
  counter = 0;
  errAccu = 0;
  jitteringLimit = 0; 
}

void ResetCounter::updateJitteringLimit()
{
  if(errAccu > 0.5) {
    jitteringLimit = intPart+1;
    errAccu -= 1.0; }
  else
    jitteringLimit = intPart;
}

//=================================================================================================

TurtleSource::TurtleSource()
{
  // init to square:
  turnAngle = 90;
  turtle.setAngle(turnAngle);
  setTurtleCommands("F+F+F+F+");
  updateWaveTable();

  setResetRatio( 0, 1);   // first resetter resets after each cycle by default
  setResetOffset(0, 0);
  for(int i = 1; i < numResetters; i++) {
    setResetRatio( i, 0); // other resetters are off by default
    setResetOffset(i, 0);
  }

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
  updateLineCommandIndices();
  updateMeanAndNormalizer();
  updateResetters();
  updateIncrement();
  reset();
}

void TurtleSource::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  updateResetters();
  incUpToDate = false;
}

void TurtleSource::setFrequency(double newFrequency)
{
  frequency = newFrequency;
  updateResetters();
  incUpToDate = false;
}

void TurtleSource::setFrequencyScaler(double newScaler)
{
  freqScaler = newScaler;
  updateResetters();
  incUpToDate = false;
}

void TurtleSource::setRotation(double newRotation)
{
  rotator.setAngle((PI/180) * newRotation);
}

void TurtleSource::setStartPosition(double newPosition)
{
  startPos = newPosition;
  startLineIndex = (int) (startPos*numLines);
  if(startLineIndex >= numLines)
    startLineIndex -= numLines;
  startCommandIndex = lineCommandIndices[startLineIndex];
}

void TurtleSource::setResetRatio(int i, double newRatio)
{
  resetRatios[i] = newRatio;
  updateResetter(i);
}

void TurtleSource::setResetOffset(int i, double newOffset)
{
  resetOffsets[i] = newOffset;
  updateResetter(i);
}

void TurtleSource::setTurnAngle(double newAngle) 
{
  turnAngle = newAngle;
  turtle.setAngle(turnAngle);
  updateMeanAndNormalizer();
  tableUpToDate = false; 
}

void TurtleSource::setSkew(double newAngle)
{
  skew = newAngle;
  turtle.setAngleDelta(skew);
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

void TurtleSource::updateLineCommandIndices()
{
  lineCommandIndices.resize(numLines);
  int j = 0;
  for(int i = 0; i < turtleCommands.size(); i++) {
    if(turtle.isLineCommand(turtleCommands[i])) {
      lineCommandIndices[j] = i;
      j++;
    }
  }
}

void TurtleSource::reset()
{
  pos = 0;
  resetTurtle();
  resetCounters();
  updateLineBuffer();

  // old:
  //lineIndex = 0; // shouldn't this be done in resetTurtle? ...but this gives a hang

  // new:
  lineIndex = startLineIndex;
}

bool TurtleSource::checkIndexConsistency()
{
  // find target command index:
  int tmp = lineCommandIndices[lineIndex];
  if(reverse) tmp -= 1;  // in reverse mode, it should be before the 'F' for current line
  else        tmp += 1;  // and in normal mode directly after the 'F'

  // ...with wrap around:
  int size = (int)turtleCommands.size();
  if(tmp >= size)  tmp -= size;
  if(tmp < 0)      tmp += size;

  // check and return result:
  bool result = commandIndex == tmp;
  rsAssert(result);
  return result;
}

void TurtleSource::resetTurtle()
{
  //// old:
  //commandIndex = 0;
  //turtle.init(0, 0, 1, 0);

  // new:
  if(!tableUpToDate)
    updateWaveTable();
  int i = startLineIndex;
  commandIndex = lineCommandIndices[i];
  double x0 = tableX[i];
  double y0 = tableY[i];
  double x1 = tableX[i+1];
  double y1 = tableY[i+1];
  turtle.init(x0, y0, x1-x0, y1-y0);

  // todo: use startPos, read out the table at the position corresponding to startPos (with linear
  // interpolation...or maybe rounding is better?), set dx,dy to the differences between the two
  // neighbouring points...oh...and also find the command-index that corresponds to the line 
  // segment...hmm...but how? maybe we need to keep a table of command-indices (indexed by 
  // line-number)
}

void TurtleSource::resetCounters()
{
  for(int i = 0; i < numResetters; i++)
    resetters[i].reset();
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
    // if we want anti-aliasing, we would actually also have to go through the intermediate 
    // segments and filter while reading out the table
  }
  else
  {
    while(lineIndex != targetLineIndex)
      goToNextLineSegment(); // increments lineIndex with wrap around

    // maybe do it like this:
    // lineIndex = targetLineIndex;
    // goToCommand(lineCommandIndices[lineIndex]);
    // but from where do i invoke the resetters then? maybe from getSampleFrame? may make more
    // sense anyway - reset instants are then not restricted to occur after a number of lines
    // has been drawn
  }

  // later: update interpolator coeffs here according to x,y buffers (but maybe only if inc > 1 in
  // which case the same set of coeffs may be used more than once)
}

void TurtleSource::goToNextLineSegment()
{
  updateLineBuffer();
  bool reset = false;
  for(int i = 0; i < numResetters; i++)
    reset |= resetters[i].tick();
  if(reset) {
    resetTurtle();
    updateLineBuffer();
  }

  if(reverse) {
    lineIndex--;
    if(lineIndex < 0)
      lineIndex = numLines-1; 
  }
  else {
    lineIndex++;
    if(lineIndex >= numLines)
      lineIndex = 0; 
  }
}

void TurtleSource::updateLineBuffer()
{
  if(turtleCommands.size() == 0)
    return;
  bool xyUpdated = false;
  while(xyUpdated == false) {
    bool draw = turtle.interpretCharacter(turtleCommands[commandIndex]);

    // increment or decrement command index:
    if(reverse) {
      commandIndex--;
      if(commandIndex < 0)
        commandIndex = (int)turtleCommands.size()-1;
    }
    else {
      commandIndex++;
      if(commandIndex == turtleCommands.size())
        commandIndex = 0;
    }

    if(draw) {

      if(!antiAlias)
      {
        if(reverse) // do we have to reverse the line buffer in reverse mode?
        {
          x[1] = turtle.getStartX();
          y[1] = turtle.getStartY();
          x[0] = turtle.getEndX();
          y[0] = turtle.getEndY();
        }
        else
        {
          x[0] = turtle.getStartX();
          y[0] = turtle.getStartY();
          x[1] = turtle.getEndX();
          y[1] = turtle.getEndY();
        }
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

void TurtleSource::updateResetters()
{
  for(int i = 0; i < numResetters; i++)
    updateResetter(i);
}

void TurtleSource::updateResetter(int i)
{
  double interval = numLines * (1/resetRatios[i] + resetOffsets[i]/(freqScaler*frequency));
  resetters[i].setInterval(interval);
  incUpToDate = false; // because computing the inc uses min(numLines, minResetInterval)
}

void TurtleSource::updateIncrement()
{
  double minLength = double(numLines);
  for(int i = 0; i < numResetters; i++)
    minLength = rmin(minLength, resetters[i].getInterval());
  inc = minLength * freqScaler * frequency / sampleRate;

  // hmm...the pitch goes down when ratio is below 1 (and the old, cyclic reset is off)
  // ...ahh - i think that's ok and expected - we should use a direction fix in the axiom to 
  // avoid this

  reverse = inc < 0.0;
  turtle.setReverseMode(reverse);
  turtleLowpass.setSampleRate(rmin(1/fabs(inc), 1.1)); // use 1.0 later
  incUpToDate = true;
}

void TurtleSource::updateMeanAndNormalizer()
{
  if(!tableUpToDate)
    updateWaveTable();

  double minX = 0, maxX = 0, minY = 0, maxY = 0, sumX = 0, sumY = 0;
  int N = (int) tableX.size();
  for(int i = 0; i < N; i++) {
    double x = tableX[i];
    double y = tableY[i];
    minX  = rmin(minX, x);
    maxX  = rmax(maxX, x);
    minY  = rmin(minY, y);
    maxY  = rmax(maxY, y);
    sumX += x;
    sumY += y;
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

  // maybe don't use the "mean" but the "center" defined as (min+max)/2 -> avoids computations and 
  // is probably just as good (especially, when a DC blocking filter is used later anyway)
}


/*
// old version - avoids creating the table:
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
*/

/*
BUGS:
-reverse seems to work even though we don't switch the turtle into reverse mode - why? ...maybe we 
 need a more asymmetric shape (like - a pentagon, and see, if it switches upside down?) - pentagon
 works too, but when the turn angle is 74 (and reset is off), it doesn't work anymore
-in non-table mode, there's a click at the beginning (try with InitSquare patch)
-BuzzingTriangles patch is different in table-mode vs non-table-mode - aahhh - i think, it is 
 because in non-table mode, and 'f' leads to a jump in position and in table-mode, the point is not
 being drawn...needs fix - in table-mode, we also need a jump...but is that possible?
-the pitch is resetting mode is not the same as in free-runing mode for closed curves - test
 with patch BuzzingTriangles (maybe add ++ to axiom, if necessarry) ...still true?
-in resetting mode, there's a click at the beginning of a note, in free-running mode, there's no 
 such click
-when start-position is > 0 and the number of iterations is reduced, we may get access violations
 (i think, it tries to reset to a line-number that has become invalid)
-for the AlienFace patch, start-point modulation gives chaotic behavior. May this have to do with
 the state-stack? If so, what can be done - maybe just disable start-point modulation when the 
 rules and/or axiom contains '[', ']'? ...may this mess up reversal behavior, too?
 ...maybe disable the feature for now...doesn't seem to be very useful anyway


Ideas:

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
-maybe reset or don't reset when the counter has reached its limit, depending the value of a 
 binary sequence ("reset sequence"). see PGCA for generation methods for binary sequences. also,
 L-system output can be used, interpreting '+' as '1' and '-' as '0'...or a dedicated L-system,
 producing 0s and 1s directly can be set up
 -maybe, it can produce more instructions than just "reset" and "don't reset", like "half reset"
  "set to" or whatever (then, the sequence wouldn't be binary anymore)
-allow for reversal of readout direction - allow negative frequencies/increments
 -introduce backward-mode to TurtleGraphics: subtracts dx,dy (do the turns also need to be 
  reversed? ...maybe, maybe not - try it: just draw F+F and go back - see, if we end up at (0,0))
-have multiple resetters and reversers that all work independently from each other - the ratios of
 them will probably be useful parameters for sound design
-have a start-phase/pos parameter (modulatable) - modulating it and reset interval at the same time
 (but with different modulators) will give interesting sounds
-Let the turtle interpret '|': turn around (180°), g: forward without line, B,b: backward with and 
 without line
-maybe dc-subtraction, normalization, rotation can be applied already to the content of the
 line buffer - that would make sounds where inc < 1 faster to compute and sounds where inc > 1 
 slower...hmm...the latter part seems bad, since they are already slow to compute...
 so maybe it's not a good idea


*/