//-------------------------------------------------------------------------------------------------
// positive range bell functions:

double rsPositiveBellFunctions::linear(double x)
{
  if(x > 1.0)
    return 0.0;
  else
    return 1.0 - x;
}

double rsPositiveBellFunctions::cubic(double x)
{
  if(x > 1.0)
    return 0.0;
  else
    return 1 + (2*x - 3) * x*x;
}

double rsPositiveBellFunctions::quintic(double x)
{
  if(x > 1.0)
    return 0.0;
  else
  {
    double x2 = x*x;
    return 1 + (-10 + 15*x - 6*x2) * x*x2; // 1 - 10*x^3 + 15*x^4 - 6*x^5
  }
}

double rsPositiveBellFunctions::heptic(double x)
{
  if(x > 1.0)
    return 0.0;
  else
  {
    double x2 = x*x;
    return 1 + (-35 + 84*x - 70*x2 + 20*x2*x) * x2*x2; // 1 - 35*x^4 + 84*x^5 - 70*x^6 + 20*x^7
  }
}

//-------------------------------------------------------------------------------------------------
// class rsParametricBellFunction:

rsParametricBellFunction::rsParametricBellFunction()
{
  bell   = &rsPositiveBellFunctions::cubic;
  center = 0.0; 
  setWidth(2.0);
  setFlatTopWidth(0.0);
}

void rsParametricBellFunction::setCenter(double newCenter)
{
  center = newCenter; 
}

void rsParametricBellFunction::setWidth(double newWidth)
{
  a = 2 / newWidth;
}

void rsParametricBellFunction::setFlatTopWidth(double newWidth)
{
  flat = newWidth;
  b = 1 / (1 - flat);
}

void rsParametricBellFunction::setPrototypeBell(double (*newFunction)(double))
{
  bell = newFunction;
}
