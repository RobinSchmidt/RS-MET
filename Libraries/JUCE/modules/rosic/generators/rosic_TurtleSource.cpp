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
  resetLineIndex(0);
  resetCounters();
}

bool TurtleSource::checkIndexConsistency()
{
  return commandIndex == lineCommandIndices[lineIndex];
}

bool TurtleSource::isInInitialState()
{
  bool r = true;

  r &= lineIndex == 0;
  r &= commandIndex == lineCommandIndices[lineIndex];

  // check line buffer and turtle state:
  double xs, ys, xe, ye;
  xs = tableX[0];
  ys = tableY[0];
  xe = tableX[1];
  ye = tableY[1];
  if(reverse && !useTable) {
    swap(xs, xe);
    swap(ys, ye); }
  r &= x[0] == xs;
  r &= y[0] == ys;
  r &= x[1] == xe;
  r &= y[1] == ye;
  r &= turtle.getStartX() == xs;
  r &= turtle.getStartY() == ys;
  r &= turtle.getEndX()   == xe;
  r &= turtle.getEndY()   == ye;

  // actually, it seems wrong to me that that the entries of the line-buffers should be swapped in
  // reverse mode - after all, we could switch between forward and reverse mode multiple times 
  // within one line segment in which case we would just go back and forth in the same line buffer

  return r;
}

void TurtleSource::resetLineIndex(double fracPos)
{
  pos = fracPos;
  lineIndex = startLineIndex;
  if(useTable)
    updateLineBufferFromTable();
  else
  {
    resetTurtle();
    goToCommand(lineCommandIndices[lineIndex]);
    updateLineBufferFromTurtle();
  }
}

void TurtleSource::resetTurtle()
{
  //// test:
  //turtle.init(0, 0, 1, 0); 
  //return; 
  //// when doing this plot 3 is right and 1,2 are wrong, when not doing it, it's the other
  //// way around

  // new since introducing startLineIndex:
  if(!tableUpToDate)
    updateWaveTable();
  int i = startLineIndex;
  commandIndex = lineCommandIndices[i];
  lineIndex = i;

  double x0 = tableX[i];
  double y0 = tableY[i];
  double x1 = tableX[i+1];
  double y1 = tableY[i+1];
  turtle.init(x0, y0, x1-x0, y1-y0);

  // initialization in reverse mode may still by buggy
  // this seems to magically work - at least with the test project, but i don't know why:
  if(reverse)
    return;

  turtle.interpretCharacter(turtleCommands[commandIndex]);
  updateLineBufferFromTurtle();
  

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
    updateLineBufferFromTable();
    // if we want anti-aliasing, we would actually also have to go through the intermediate 
    // segments and filter while reading out the table and filter during this process
  }
  else
  {
    lineIndex = targetLineIndex;

    // new:
    goToCommand(lineCommandIndices[lineIndex]); 

    //// old:
    //bool reset = false;
    //for(int i = 0; i < numResetters; i++)
    //  reset |= resetters[i].tick();
    //if(reset)
    //  resetTurtle();
    //else
    //  goToCommand(lineCommandIndices[lineIndex]);
    // maybe invoke resetting from getSampleFrame - may make more sense anyway - reset instants are
    // then not restricted to occur after a number of lines has been drawn (we need to adapt the 
    // reset-intervals then)..also, we may want to have resetting in table-mode, too
    // actually, it's wrong here anyway - we would have to do as mayn reset as the difference 
    // between targetLineIndex and lineIndex
  }

  // later: update interpolator coeffs here according to x,y buffers (but maybe only if inc > 1 in
  // which case the same set of coeffs may be used more than once)
}

void TurtleSource::goToCommand(int targetCommandIndex)
{
  if(turtleCommands.size() == 0)
    return;

  while(commandIndex != targetCommandIndex)
  {
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

    bool lineDrawn = turtle.interpretCharacter(turtleCommands[commandIndex]);
    if(lineDrawn)
      updateLineBufferFromTurtle();
    // todo: take anti-alias flag into account - if it's false, it may actually suffice to update
    // the line buffers after leaving the loop...hmm...maybe we can generally fill them after the 
    // loop keeping track only of some temporary variables here

  }
}

void TurtleSource::updateLineBufferFromTurtle()
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
  

  /*
  // test - don't reverse line-buffer - nope, seems wrong:
  x[0] = turtle.getStartX();
  y[0] = turtle.getStartY();
  x[1] = turtle.getEndX();
  y[1] = turtle.getEndY();
  */
}

void TurtleSource::updateLineBufferFromTable()
{
  int i = lineIndex;
  x[0] = tableX[i];
  y[0] = tableY[i];
  x[1] = tableX[i+1];
  y[1] = tableY[i+1];
}

void TurtleSource::updateResetters()
{
  for(int i = 0; i < numResetters; i++)
    updateResetter(i);
}

void TurtleSource::updateResetter(int i)
{
  double interval = numLines * (1/resetRatios[i] + resetOffsets[i]/(freqScaler*frequency));
  //double interval = (1/resetRatios[i] + resetOffsets[i]/(freqScaler*frequency));
  resetters[i].setInterval(interval);
  incUpToDate = false; // because computing the inc uses min(numLines, minResetInterval)
}

void TurtleSource::updateIncrement()
{
  double oldInc = inc;

  // old (before dragging reset into getSample):
  double minLength = double(numLines);
  for(int i = 0; i < numResetters; i++)
    minLength = rmin(minLength, resetters[i].getInterval());
  inc = minLength * freqScaler * frequency / sampleRate;

  // new:
  for(int i = 0; i < numResetters; i++)
    resetters[i].setIncrement(inc);


  reverse = inc < 0.0;

  turtle.setReverseMode(reverse);
  if(inc*oldInc < 0)    // a direction change occured...
    turtle.backtrack(); // ...turtle needs to set its position back to the start of current line


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
-in reverse mode, the picture is shifted when a new note is started (something wrong with reset?)
 -it always shifts to the left, it seems to shift less when there are more lines
 -seems to be fixed - but i don't know, why it works
-the resetting does not seem to work anymore (wrong interval?)...maybe move it to getSample anyway
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