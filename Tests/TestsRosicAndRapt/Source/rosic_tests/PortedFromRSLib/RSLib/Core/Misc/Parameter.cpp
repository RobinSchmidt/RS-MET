#include "Parameter.h" // remove later

// construction/destruction:

rsNumericParameter::rsNumericParameter()
{
  toString        = &rsDoubleToString2;
  value           = 0.0;
  defaultValue    = 0.0;
  normalizedValue = 0.0;
  mapper          = new rsRealLinearMap;
}

rsNumericParameter::rsNumericParameter(const rsString& name, double minValue, double maxValue, 
  double defaultValue, DoubleToStringFunctionPointer toStringFunction, int scaling, 
  const rsString& unit)
{
  this->name         = name;
  this->unit         = unit;
  this->toString     = toStringFunction;
  this->defaultValue = defaultValue;
  switch( scaling )
  {
  case EXPONENTIAL: mapper = new rsRealLinToExpMap(0.0, minValue, 1.0, maxValue); break;
  default:          mapper = new rsRealLinearMap(  0.0, minValue, 1.0, maxValue);
  }
  setValue(defaultValue);
}

rsNumericParameter::~rsNumericParameter()
{
  deleteAllCallbacks();
  delete mapper;
}

// setup:

void rsNumericParameter::addCallback(rsCallbackBase1<void, double> *callbackToAdd)
{
  callbacks.push_back(callbackToAdd);
}

void rsNumericParameter::deleteAllCallbacks()
{  
  for(unsigned int i = 0; i < callbacks.size(); i++)
    delete[] callbacks[i];
  callbacks.clear();
}

void rsNumericParameter::setNormalizedValue(double newNormalizedValue)
{
  normalizedValue = rsClipToRange(newNormalizedValue, 0.0, 1.0);
  value = mapper->evaluate(normalizedValue);
  invokeCallbacks();
}

void rsNumericParameter::setValue(double newValue)
{
  setNormalizedValue(mapper->evaluateInverse(newValue));
}

void rsNumericParameter::setMapper(rsRealFunctionInvertible *newMapper)
{
  rsAssert( newMapper != mapper );
  delete mapper;
  mapper = newMapper;
}

// misc:

void rsNumericParameter::invokeCallbacks()
{
  for(unsigned int i = 0; i < callbacks.size(); i++)
    callbacks[i]->call(value);
}

