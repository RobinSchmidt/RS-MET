#ifndef RSLib_FunctionIterators_cpp  // do we need this???
#define RSLib_FunctionIterators_cpp

using namespace RSLib;

rsComplexExponentialIterator::rsComplexExponentialIterator(rsComplexDbl a, rsComplexDbl z)
{
  this->w = a;
  this->z = z;
}

void rsComplexExponentialIterator::resetValue(rsComplexDbl initialValue)
{ 
  w = initialValue;
}

void rsComplexExponentialIterator::setZ(rsComplexDbl newZ)
{
  z = newZ;
}


    
rsSineIterator::rsSineIterator()
{
  a1 =  1.0806046117362795;
  s1 = -0.84147098480789650;
  s2 = -0.90929742682568171;
    // calling setup(1, 0, 1) would compute these values, but that would be more costly.
}

rsSineIterator::rsSineIterator(double w, double p, double a)
{
  setup(w, p, a);
}

void rsSineIterator::setup(double w, double p, double a)
{
  a1 = 2.0*cos(w);
  s1 = a*sin(p-    w);
  s2 = a*sin(p-2.0*w);
}




#endif
