//#include "rosic_PiecewiseFunction.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

PiecewiseFunction::PiecewiseFunction()
{
  minY                  = -10.0;
  maxY                  = +10.0;
  minBreakpointDistance = 0.001;
  interpolationMethod   = LINEAR;
  a = b = c = d         = NULL;
  numSegments           = 0;
  initialize();
}

PiecewiseFunction::~PiecewiseFunction()
{
  delete[] a;
  delete[] b;
  delete[] c;
  delete[] d;
}

void PiecewiseFunction::copyDataFrom(PiecewiseFunction &source)
{
  // copy everything except the mutex-member:
  mutex.lock();
  source.mutex.lock();

  breakpoints           = source.breakpoints;
  minY                  = source.minY;
  maxY                  = source.maxY;
  minBreakpointDistance = source.minBreakpointDistance;
  interpolationMethod   = source.interpolationMethod;

  if( source.numSegments != this->numSegments )
  {
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
    numSegments = source.numSegments;
    a = new double[numSegments+1]; // +1 for a_N = f_N
    b = new double[numSegments];
    c = new double[numSegments+1]; // +1 for c_N = 0
    d = new double[numSegments];
  }
  for(int s=0; s<numSegments; s++)
  {
    a[s] = source.a[s];
    b[s] = source.b[s];
    c[s] = source.c[s];
    d[s] = source.d[s];
  }
  a[numSegments] = source.a[numSegments];
  c[numSegments] = source.c[numSegments];

  mutex.unlock();
  source.mutex.unlock();
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void PiecewiseFunction::setInterpolationMethod(int newInterpolationMethod)
{
  if( newInterpolationMethod >= 0 && newInterpolationMethod < NUM_INTERPOLATION_METHODS )
    interpolationMethod = newInterpolationMethod;
}

void PiecewiseFunction::setMinimumAllowedY(double newMinY)
{
  if( newMinY < maxY-0.001 )
    minY = newMinY;
  for(int i=0; i <= lastBreakpointIndex(); i++)
  {
    if( breakpoints[i].y < minY )
      breakpoints[i].y = minY;
  }
  calculateCubicCoefficients();
}

void PiecewiseFunction::setMaximumAllowedY(double newMaxY)
{
  if( newMaxY > minY+0.001 )
    maxY = newMaxY;

  for(int i=0; i <= lastBreakpointIndex(); i++)
  {
    if( breakpoints[i].y > maxY )
      breakpoints[i].y = maxY;
  }
  calculateCubicCoefficients();
}

int PiecewiseFunction::insertBreakpoint(double newX, double newY)
{
  if( newX >= getMinX() )
  {
    // Loop through the existing breakpoints and inspect their x-values. As soon as this x-value
    // gets (strictly) larger than the x--value of the new breakpoint to be inserted, we know that
    // this is the index, right before which the new breakpoint has to be inserted:
    int    bpIndex = 1;  // omit the first breakpoint(index 0)
    double bpX     = 0.0;
    while( bpIndex < (int) breakpoints.size() )
    {
      bpX = breakpoints[bpIndex].x;
      if( bpX > newX )
      {
        // make sure that the new  breakpoint to be inserted is not too close to
        // another existing breakpoint:
        if( newX - breakpoints[bpIndex-1].x < minBreakpointDistance
          || bpX - newX                     < minBreakpointDistance )
        {
          return -1;
        }

        // create a new breakpoint for insertion right before the breakpoint with the current
        // index:
        FunctionBreakpoint newBreakpoint;
        newBreakpoint.x = newX;
        newBreakpoint.y = clip(newY, minY, maxY);

        // insert the new breakpoint right before the breakpoint with the current index (this must
        // be done inside a mutex-lock):
        mutex.lock();
        breakpoints.insert(breakpoints.begin()+bpIndex, 1, newBreakpoint);
        calculateCubicCoefficients();
        mutex.unlock();

        // nothing more to do here. jump out of the function:
        return bpIndex;
      } // end of  if( bpTime >= newTimeStamp )
      else
        bpIndex++; // lets look at the next breakpoint...

    } // end of while( bpIndex < (int) breakpoints.size() )

    return -1; // this command should actually never be executed due to the
    // logic of the stuff above, but we need it to suppress a compiler-warning

  } // end of  if( newX >= 0.0 )
  else
    return -1; // breakpoints with newX < 0 are illegal  ....get rid of this !!!!
}

bool PiecewiseFunction::removeBreakpoint(int index)
{
  // do not allow removal of the first and last breakpoint:
  if( index > 0 && index < lastBreakpointIndex() )
  {
    mutex.lock();
    breakpoints.erase(breakpoints.begin()+index);
    calculateCubicCoefficients();
    mutex.unlock();
    return true;
  }
  else
    return false;
}

bool PiecewiseFunction::modifyBreakpoint(int index, double newX, double newY)
{
  if( index >= 0 && index < (int) breakpoints.size()  )
  {
    if( index == lastBreakpointIndex() )
    {
      // make sure that the breakpoint does not come too close to its neighbours:
      if( newX - breakpoints[index-1].x < minBreakpointDistance )
        newX = breakpoints[index-1].x + minBreakpointDistance;

      breakpoints[index].x = newX;
      breakpoints[index].y = clip(newY, minY, maxY);
      calculateCubicCoefficients();

      return true;
    }
    else // we deal with an intermediate breakpoint
    {

      // make sure that the breakpoint does not come too close to its neighbours:
      if( index > 0 )
      {
        if( newX - breakpoints[index-1].x < minBreakpointDistance )
          newX = breakpoints[index-1].x + minBreakpointDistance;
      }
      if( breakpoints[index+1].x - newX < minBreakpointDistance )
        newX = breakpoints[index+1].x - minBreakpointDistance;

      breakpoints[index].x = newX;
      breakpoints[index].y = clip(newY, minY, maxY);
      calculateCubicCoefficients();

      return true;
    } // end of if(index==0), else if(index==lastBreakpointIndex()), else

  } // end of if( index >= 0 && index < (int) data->breakpoints.size()  )
  else
    return false;
}

//-------------------------------------------------------------------------------------------------
// inquiry:

int PiecewiseFunction::getNumBreakpoints() const
{
  return (int) breakpoints.size();
}

int PiecewiseFunction::lastBreakpointIndex() const
{
  return (int) breakpoints.size()-1;
}

double PiecewiseFunction::getMinX()
{
  mutex.lock();
  double min = breakpoints[0].x;
  mutex.unlock();
  return min;
}

double PiecewiseFunction::getMaxX()
{
  mutex.lock();
  double max = breakpoints[lastBreakpointIndex()].x;
  mutex.unlock();
  return max;
}

double PiecewiseFunction::getMinY()
{
  mutex.lock();
  double min = breakpoints[0].y;
  for(int p=0; p < (int) breakpoints.size(); p++)
  {
    if( breakpoints[p].y < min )
      min = breakpoints[p].y;
  }
  mutex.unlock();
  return min;
}

double PiecewiseFunction::getMaxY()
{
  mutex.lock();
  double max = breakpoints[0].y;
  for(int p=0; p < (int) breakpoints.size(); p++)
  {
    if( breakpoints[p].y > max )
      max = breakpoints[p].y;
  }
  mutex.unlock();
  return max;
}

double PiecewiseFunction::getBreakpointX(int index) const
{
  if( index < 0 || index > (int) breakpoints.size()-1 )
    return NAN;
  else
    return breakpoints[index].x;
}

double PiecewiseFunction::getBreakpointMinX(int index) const
{
  if( index < 0 || index > (int) breakpoints.size()-1 )
    return NAN;
  else
  {
    if( index == 0 )
      return breakpoints[index].x - 1.0;
    else
      return breakpoints[index-1].x + minBreakpointDistance;
  }
}

double PiecewiseFunction::getBreakpointMaxX(int index) const
{
  if( index < 0 || index > (int) breakpoints.size()-1 )
    return NAN;
  else
  {
    if( index == lastBreakpointIndex() )
      return 30.0; // this is somewhat arbitrary
    else
      return breakpoints[index+1].x - minBreakpointDistance;
  }
}

bool PiecewiseFunction::setBreakpointX(int index, double newX)
{
  // check if the first and last data->breakpoints are being modified and impose
  // some restrictions on them....
  if( index == 0 )
  {
    if( breakpoints[index+1].x-newX < minBreakpointDistance )
      newX = breakpoints[index+1].x - minBreakpointDistance;
    breakpoints[index].x = newX;
    calculateCubicCoefficients();
    return true;
  }
  else if( index == lastBreakpointIndex() )
  {
    if( index<1 ) return false; // should never happen -> would result in an access violation!

    // make sure that the breakpoint does not come too close to its predecessor:
    if( newX - breakpoints[index-1].x < minBreakpointDistance )
      newX = breakpoints[index-1].x + minBreakpointDistance;
    breakpoints[index].x = newX;
    calculateCubicCoefficients();
    return true;
  }
  else if( index > 0 && index < lastBreakpointIndex() )  // we deal with an intermediate breakpoint
  {
    // make sure that the breakpoint does not come too close to its neighbours:
    if( newX - breakpoints[index-1].x < minBreakpointDistance )
      newX = breakpoints[index-1].x + minBreakpointDistance;
    if( breakpoints[index+1].x - newX < minBreakpointDistance )
      newX = breakpoints[index+1].x - minBreakpointDistance;
    breakpoints[index].x = newX;
    calculateCubicCoefficients();
    return true;
  } // end of if(index==0), else if(index==lastBreakpointIndex()), else
  else
    return false;
}

double PiecewiseFunction::getBreakpointY(int index) const
{
  if( index < 0 || index > (int) breakpoints.size()-1 )
    return 0.0;
  else
    return breakpoints[index].y;
}

bool PiecewiseFunction::setBreakpointY(int index, double newY)
{
  if( index >= 0 && index <= lastBreakpointIndex() )
  {
    breakpoints[index].y = clip(newY, minY, maxY);
    calculateCubicCoefficients();
    return true;
  }
  else
    return false;
}

//-------------------------------------------------------------------------------------------------
// others:

int PiecewiseFunction::getNextNonCoincidingIndex(int startIndex)
{
  //data->mutex.lock();
  // mutex is not necesarry - this function is called only internally from functions which
  // already aquire the lock

  // when there are no breakpoints at the same time instant, this will be the default value:
  int foundIndex = startIndex + 1;

  // make sure to not access out-of-range indices:
  if( foundIndex > lastBreakpointIndex() )
    return 0;

  while( breakpoints[startIndex].x >= breakpoints[foundIndex].x )
  {
    // timeStamp at foundIndex is not strictly larger than at startIndex - skip to next breakpoint:
    foundIndex++;

    // make sure to not access out-of-range indices:
    if( foundIndex > lastBreakpointIndex() )
      return 0;
  }

  // O.K. we now have skipped to the first index for which the timeStamp is actually larger than
  // at the startIndex (and returned 0, if that would be out-of-range) - now we can return our
  // finding:
  return foundIndex;
}

int PiecewiseFunction::getLastCoincidingIndex(int index)
{
  // make sure to not access out-of-range indices:
  if( index+1 > lastBreakpointIndex() )
    return lastBreakpointIndex();

  while( breakpoints[index].x >= breakpoints[index+1].x )
  {
    // timeStamp at foundIndex is not strictly larger than at startIndex - skip to next breakpoint:
    index++;

    // make sure to not access out-of-range indices:
    if( index+1 > lastBreakpointIndex() )
      return lastBreakpointIndex();
  }

  return index;
}

void PiecewiseFunction::initialize()
{
  mutex.lock();

  // initialize the breakpoint-vector with two entries, these two will always be there (their data
  // can be modified, though), additional entries can be inserted and removed at will in between:
  breakpoints.clear();
  FunctionBreakpoint newBreakpoint;

  newBreakpoint.x = -1.0;
  newBreakpoint.y = -1.0;
  breakpoints.push_back(newBreakpoint);

  newBreakpoint.x = +1.0;
  newBreakpoint.y = +1.0;
  breakpoints.push_back(newBreakpoint);

  calculateCubicCoefficients();

  mutex.unlock();
}

void PiecewiseFunction::calculateCubicCoefficients()
{
  mutex.lock();

  if( (int) breakpoints.size() != numSegments-1 )
  {
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
    numSegments = (int) breakpoints.size()-1;
    a = new double[numSegments+1];  // +1 for a_N = f_N
    b = new double[numSegments];
    c = new double[numSegments+1];  // +1 for c_N = 0
    d = new double[numSegments];
  }

  double *h = new double[numSegments];  // segment lengths
  for(int i=0; i<numSegments; i++)
  {
    a[i] = breakpoints[i].y;                        // constant for segment i (Eq. 19.235)
    h[i] = breakpoints[i+1].x - breakpoints[i].x;   // length of segment i
  }
  a[numSegments] = breakpoints[numSegments].y;

  // compute the right-hand-sides for the tridiagonal equation system:
  double rhs[10];
  for(int i=1; i<numSegments; i++)
    rhs[i] = 3 * ( (a[i+1]-a[i])/h[i] - (a[i]-a[i-1])/h[i-1] );  // Eq. 19.239

  // establish the main-diagonal for the tridiagonal system (Eq.19.238):
  double md[10];
  for(int i=0; i<numSegments-1; i++)
    md[i] = 2 * (h[i]+h[i+1]);

  // solve the system (resulting in our c-coefficients) and set c[0] and c[N-1] zero for
  // a 'natural' cubic spline (with zero second derivative at the endpoints):
  solveTriDiagonalSystem(&h[1], md, &h[1], &rhs[1], &c[1], (int) breakpoints.size()-2); // Eq.19.240
  c[0]                    = 0;
  c[breakpoints.size()-1] = 0;

  // use the c-coefficients to solve for the b- and d-coefficients
  for(int i=1; i<=numSegments; i++)
    b[i-1] = (a[i]-a[i-1])/h[i-1] - (2*c[i-1]+c[i]) * h[i-1] / 3;   // Eq.19.238
  for(int i=1; i<=numSegments; i++)
    d[i-1] = (c[i]-c[i-1]) / (3*h[i-1]);                            // Eq.19.237

  delete[] h;

  mutex.unlock();
}

