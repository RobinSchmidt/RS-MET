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
  tableUpToDate = false;
  incUpToDate = false;

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

  // reversers are all off by default:
  for(int i = 0; i < numReversers; i++) {
    setReverseRatio( i, 0);
    setReverseOffset(i, 0);
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

void TurtleSource::setReverseRatio(int i, double newRatio)
{
  reverseRatios[i] = newRatio;
  updateReverser(i);
}

void TurtleSource::setReverseOffset(int i, double newOffset)
{
  reverseOffsets[i] = newOffset;
  updateReverser(i);
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
  RS_ASSERT(tableX.size() == tableY.size());           // tables for x and y must have same size
  while(tableX.size() < 2) {
    tableX.push_back(0);
    tableY.push_back(0);
  }
  RS_ASSERT(tableX.size() >= 2 && tableY.size() >= 2); // tables need to have at least 2 values
  tableUpToDate = true;
}

void TurtleSource::updateLineCommandIndices()
{
  //lineCommandIndices.resize(numLines);
  lineCommandIndices.resize(RAPT::rsMax(numLines,1));
  int j = 0;
  for(int i = 0; i < (int)turtleCommands.size(); i++) {
    if(turtle.isLineCommand(turtleCommands[i])) {
      lineCommandIndices[j] = i;
      j++;
    }
  }
}

void TurtleSource::reset()
{
  resetPhase();
  resetCounters();
  reverseFlipFlop = false;
  updateDirection();
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
  r &= reverseFlipFlop == false;

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

void TurtleSource::resetPhase(double newPhase)  
{
  // rename newPhase to phaseAdvance or get rid - i wasn't there befroe and was inetroduced to 
  // handle the andvance when receiving resets in getSample - but that doesn't seem to work well. 
  // we may need a spearate resetPhaseWithAdvance function in the subclass

  // old:
  if(newPhase == 0.0) {
    pos = 0;
    lineIndex = 0; }
  else {
    pos = newPhase;
    lineIndex = floorInt(getLinePosition(pos)); }

  // new:
  //pos = phaseOffset + newPhase;
  //lineIndex = floorInt(getLinePosition(pos));

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
  int i = 0;
  commandIndex = lineCommandIndices[i];
  lineIndex = i;

  //rsAssert(tableX.size() >= 2 && tableY.size() >= 2);
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
}

void TurtleSource::resetCounters()
{
  for(int i = 0; i < numResetters; i++)
    resetters[i].reset();
  for(int i = 0; i < numReversers; i++)
    reversers[i].reset();
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
    goToCommand(lineCommandIndices[lineIndex]);
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

double TurtleSource::computeInterval(double ratio, double offset)
{
  double freq = freqScaler*frequency;
  if(ratio == 0 || freq == 0)
    return RS_INF(double);       // avoid Nan in potential inf + inf computation
  return 1/ratio + offset/freq;
}

void TurtleSource::updateResetters()
{
  for(int i = 0; i < numResetters; i++)
    updateResetter(i);
}

void TurtleSource::updateResetter(int i)
{
  //double interval = (1/resetRatios[i] + resetOffsets[i]/(freqScaler*frequency));
  double interval = computeInterval(resetRatios[i], resetOffsets[i]);
  resetters[i].setInterval(interval);
  incUpToDate = false; // because computing the inc uses min(numLines, minResetInterval)
}

void TurtleSource::updateReversers()
{
  for(int i = 0; i < numReversers; i++)
    updateReverser(i);
}

void TurtleSource::updateReverser(int i)
{
  //double interval = 0.5*(1/reverseRatios[i] + reverseOffsets[i]/(freqScaler*frequency));
  double interval = 0.5*computeInterval(reverseRatios[i], reverseOffsets[i]);
  reversers[i].setInterval(interval);
  incUpToDate = false;

  // 0.5 * because the resulting periodicity is twice the interval of reversals
}

void TurtleSource::reverseDirection()
{
  reverseFlipFlop = !reverseFlipFlop;
  updateDirection();
}

bool rsXor(bool a, bool b) // move somewhere else
{
  return (a || b) && !(a && b); // can this be optimized?
}

void TurtleSource::updateDirection()
{
  bool oldReverse = reverse;
  reverse = rsXor(reverseFlipFlop, (freqScaler*frequency < 0));
  if(reverse && inc > 0 || !reverse && inc < 0) // can this be optimized?
    inc = -inc;
  turtle.setReverseMode(reverse);
  if(rsXor(reverse, oldReverse))  // a direction change occured - turtle needs to set its position
    turtle.backtrack();           // back to the start of current line
}

void TurtleSource::updateIncrement()
{
  double oldInc = inc;
  double minLength = 1;
  for(int i = 0; i < numResetters; i++)
    minLength = RAPT::rsMin(minLength, resetters[i].getInterval());

  for(int i = 0; i < numReversers; i++)
    minLength = RAPT::rsMin(minLength, 2*reversers[i].getInterval());



  inc = minLength * freqScaler * frequency / sampleRate;

  for(int i = 0; i < numResetters; i++)
    resetters[i].setIncrement(fabs(inc));

  for(int i = 0; i < numReversers; i++)
    reversers[i].setIncrement(fabs(inc));


  // new:
  updateDirection();

  /*
  // old:
  reverse = inc < 0.0;
  turtle.setReverseMode(reverse);
  if(inc*oldInc < 0)    // a direction change occured...
    turtle.backtrack(); // ...turtle needs to set its position back to the start of current line
  */


  turtleLowpass.setSampleRate(RAPT::rsMin(1/fabs(inc), 1.1)); // use 1.0 later
  incUpToDate = true;
}

void TurtleSource::updateMeanAndNormalizer()
{
  if(!tableUpToDate)
    updateWaveTable();

  double minX = 0, maxX = 0, minY = 0, maxY = 0;
  int N = (int) tableX.size();
  for(int i = 0; i < N; i++) {
    minX  = RAPT::rsMin(minX, tableX[i]);
    maxX  = RAPT::rsMax(maxX, tableX[i]);
    minY  = RAPT::rsMin(minY, tableY[i]);
    maxY  = RAPT::rsMax(maxY, tableY[i]);
  }

  centerX = 0.5*(minX+maxX);
  centerY = 0.5*(minY+maxY);

  minX -= centerX; maxX -= centerX;
  minY -= centerY; maxY -= centerY;
  maxX  = RAPT::rsMax(fabs(minX), fabs(maxX));
  maxY  = RAPT::rsMax(fabs(minY), fabs(maxY));
  normalizer = 1.0 / RAPT::rsMax(maxX, maxY);

  // in free running mode, this so computed mean is sometimes (often) wrong because one cycle is
  // only part (for example half) of a larger shape - in this case, we would need the mean over two
  // cycles ...maybe we can somehow deduce the symmetry in order to compute a better mean? maybe we
  // should offer different normalization modes...maybe have a highpass and leveller running on it
  // in realtime (but not oversampled)

  // maybe don't use the "mean" but the "center" defined as (min+max)/2 -> avoids computations and
  // is probably just as good (especially, when a DC blocking filter is used later anyway)
}

//=================================================================================================

void TurtleSourceAntiAliased::getSampleFrameStereoAA(double* outL, double* outR)
{
  // some checks (optimize - have a single readyToPlay flag so we only need one check here):
  if(numLines < 1)                return;
  if(!tableUpToDate && useTable)  updateWaveTable();
  if(!incUpToDate)                updateIncrement(); // must be done before goToLineSegment
  updatePosition();

  // handle periodic resetting:
  double slopeChangeX, slopeChangeY;
  double blepTime;
  for(int i = 0; i < numResetters; i++)
  {
    bool shouldReset = resetters[i].tick();
    if(shouldReset)
    {
      double newPhase = resetters[i].getPosition();
      blepTime = newPhase / inc;
      double stepX, stepY;
      resetPhase(newPhase, &stepX, &stepY, &slopeChangeX, &slopeChangeY);
      xBlep.prepareForStep(  blepTime, stepX);
      yBlep.prepareForStep(  blepTime, stepY);
      xBlep.prepareForCorner(blepTime, slopeChangeX); // commented for debuging the steps
      yBlep.prepareForCorner(blepTime, slopeChangeY);
    }
  }
  // this seems to work if only one resetter is active, but for more than one, maybe we need to 
  // first compute the data for all (both) resetters, then sort them by order of occurrence and 
  // then execute them? or will the first reset obviate all others? ...that would be convenient
  // ...but we may still need to figure out, which one comes first
  // ...or maybe just use a single resetter for the time being
 

  /*
  // handle periodic direction reversal:
  bool shouldReverse = false;
  for(int i = 0; i < numReversers; i++)
    shouldReverse |= reversers[i].tick();
  if(shouldReverse)
    reverseDirection();
    */
  // reversers should be handled similarly to resetters, they should invert the slope, so a blamp
  // should be inserted


  // integer and fractional part of position:
  double linePos = getLinePosition(pos); // linePos = 0...numLines
  // maybe this should just use linePos = pos*numLines

  int iPos = floorInt(linePos);          // iPos = 0...numLines...maybe subtract 1 in case == numLines?
  double fPos = linePos - iPos;

  // jump to appropriate line segment, thereby prepare bleps for the slope change:
  if(iPos != lineIndex)
  {
    goToLineSegment(iPos, &slopeChangeX, &slopeChangeY);
    blepTime = fPos / (numLines*inc);      // OPTIMIZE: precompute 1/(numLines*inc)
    xBlep.prepareForCorner(blepTime, slopeChangeX);
    yBlep.prepareForCorner(blepTime, slopeChangeY);
  }

  // read out buffered line segment (not yet anti-aliased):
  double x, y;
  interpolate(&x, &y, fPos);

  // apply anti-aliasing:
  x = xBlep.getSample(x);
  y = yBlep.getSample(y);

  // apply some final scaling and rotation:
  double a = normalizer * amplitude;  // maybe precompute as finalAmplitude
  *outL = a * (x - centerX);
  *outR = a * (y - centerY);
  rotator.apply(outL, outR);
}

void TurtleSourceAntiAliased::goToLineSegment(int targetLineIndex,
  double* slopeChangeX, double* slopeChangeY)
{
  double s = inc * numLines;
  double oldSlopeX = x[1] - x[0];
  double oldSlopeY = y[1] - y[0];
  Base::goToLineSegment(targetLineIndex);
  double newSlopeX = x[1] - x[0];   // get rid
  double newSlopeY = y[1] - y[0];   // get rid
  *slopeChangeX = s * (newSlopeX - oldSlopeX);
  *slopeChangeY = s * (newSlopeY - oldSlopeY);
}

void TurtleSourceAntiAliased::resetPhase(double targetPhase, double* stepX, double* stepY, 
  double* slopeChangeX, double* slopeChangeY)
{
  double linePos, fPos, oldX, oldY, /*newX, newY,*/ oldSlopeX, oldSlopeY, s;
  int iPos;

  linePos = getLinePosition(pos);
  iPos    = floorInt(linePos);
  fPos    = linePos - iPos;
  interpolate(&oldX, &oldY, fPos);
  oldSlopeX = x[1] - x[0];
  oldSlopeY = y[1] - y[0];

  // oldX, oldY already contain an additional full step by a full increment! We must actually 
  // subtract targetPhase * numLines * oldSlope as correction to set it back to a partial step:

  s = -targetPhase * numLines;

  //// test to fix problem when resetRatio = 0.25:
  //double tol = 1.e-11;
  //if(targetPhase <= tol)  s = 1;  // is this correct?

  oldX += s * oldSlopeX;
  oldY += s * oldSlopeY;

  if(fPos == 0)   // edge case, occurs when linePos == floor(linePos)
  //if(fPos <= tol || fPos >= 1-tol)   // hmmm ... doesn't seem right
  { 
    oldX = x[1]; 
    oldY = y[1]; 
  } 

  Base::resetPhase(targetPhase);


  // experimental:
  double linePos2 = getLinePosition(pos);
  int    iPos2    = floorInt(linePos);
  if(iPos == iPos2)
  {
    *stepX = *stepY = 0; 
    //*slopeChangeX = *slopeChangeY = 0;
    //return;
  }
  else
  {  
    *stepX = x[0] - oldX;
    *stepY = y[0] - oldY;
  }
  // this seems to fix the spurious spikes that occur when a wrap-around occurs immediately before
  // a reset...but why is this needed? should the result of the computaions below yield zero in 
  // this case? after all, when we are in the same segment pre and post Base::resetPhase, 
  // x[0], x[1], y[0], y[1] should be the same before and after...but they do not seem to be
  // ...soo - why does resetPhase change the x,y buffers in this case
  // oh: but i think the s = -targetPhase * numLines; applies only if we really had a 
  // reset that actually did an additional wraparound - otherwise, it makes no sense
  // -this fix feels like a dirty workaround, but maybe it's indeed the way to go - we'll see
  // -it seems to create dirt when the start-phase is not zero





  s = inc * numLines;
  *slopeChangeX = s * (x[1] - x[0] - oldSlopeX);
  *slopeChangeY = s * (y[1] - y[0] - oldSlopeY);

  int dummy = 0;
}
// -needs more verification with different resetRatios


//=================================================================================================

/*
Features to do:
-anti-aliasing
-continuous number of iterations (do it in Snowflake)
 -we need TurtleSourceMulti that
  either: runs all turtles in parallel to ensure that they all stay in sync (expensive)
  or: only run two turtles at a time and figure out a way to meaningfully initialize the state when
  we have to switch on a new turtle (like from the pair 2/3 to 3/4, we have to switch on 4)
  -for initialization of turtle 4 we may use its stored table and/or the state of "neighbour"
   turtle 3


BUGS:
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
