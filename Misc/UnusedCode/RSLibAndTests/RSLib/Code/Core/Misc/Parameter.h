#ifndef RS_PARAMETER_H
#define RS_PARAMETER_H

namespace RSLib
{

  /** Abstract baseclass for representing functions of the real numbers. */

  class RSLib_API rsRealFunction
  {

  public:

    virtual double evaluate(double x) = 0;

  };

  /** Abstract baseclass for representing invertible functions of the real numbers. */

  class RSLib_API rsRealFunctionInvertible : public rsRealFunction
  {

  public:

    virtual double evaluateInverse(double y) = 0;

  };

  /** A function object that represents a linear mapping from the real numbers to the real numbers 
  of the general form: y = a*x + b. */

  class RSLib_API rsRealLinearMap : public rsRealFunctionInvertible
  {

  public:

    rsRealLinearMap(double a = 1.0, double b = 0.0)
    {
      this->a = a;
      this->b = b;
    }

    rsRealLinearMap(double x1, double y1, double x2, double y2)
    {
      setParametersFromTwoPoints(x1, y1, x2, y2);
    }

    void setParameters(double newA, double newB)
    {
      a = newA;
      b = newB;
    }

    void setParametersFromTwoPoints(double x1, double y1, double x2, double y2)
    {
      a = (y2-y1) / (x2-x1);
      b = y1 - a*x1;
    }

    virtual double evaluate(double x) 
    {
      return a*x + b;
    }

    virtual double evaluateInverse(double y)
    {
      return (y-b) / a;
    }

  protected:

    double a, b;

  };

  /** A function object that represents an exponential mapping from the real numbers to the real 
  numbers of the general form: y = a * exp(b*x). */

  class RSLib_API rsRealLinToExpMap : public rsRealFunctionInvertible
  {

  public:

    rsRealLinToExpMap(double a = 1.0, double b = 1.0)
    {
      this->a = a;
      this->b = b;
    }

    rsRealLinToExpMap(double x1, double y1, double x2, double y2)
    {
      setParametersFromTwoPoints(x1, y1, x2, y2);
    }

    void setParameters(double newA, double newB)
    {
      a = newA;
      b = newB;
    }

    void setParametersFromTwoPoints(double x1, double y1, double x2, double y2)
    {
      b = log(y2/y1) / (x2-x1);
      a = y1 / exp(b*x1);
    }

    virtual double evaluate(double x) 
    {
      return a * exp(b*x);
    }

    virtual double evaluateInverse(double y)
    {
      return log(y/a) / b;
    }

  protected:

    double a, b;

  };

  // \todo move these rsRealFunction classes into the Math module and the rsNumericParameter class
  // into a Control module


  /**

  This is a class for representing numeric parameters for applications and plugins. The class 
  facilitates mapping between normalized and real-world values as well as invocation of callbacks 
  that shall be called whenever the parameter changes its value.

  \todo: maybe move this into a "Control" module
  \todo: make classes rsChoiceParameter, rsBoolParameter, rsStringParameter
  \todo: maybe rename this to NumericCallbackParameter and provide an alternative implementation
         that uses an observer pattern

  */

  class RSLib_API rsNumericParameter
  {

  public:

    /** Enumeration of the available scaling modes. */
    enum scalings
    {  
      LINEAR = 0,
      EXPONENTIAL
    };


    /** \name Construction/Destruction */

    /** Constructor. */
    rsNumericParameter();

    /** Constructor with intializations. */
    rsNumericParameter(const rsString& name, double minValue, double maxValue, double defaultValue,
      DoubleToStringFunctionPointer toStringFunction, int scaling, const rsString& unit);


    /** Destructor. */
    ~rsNumericParameter();


    /** \name Setup */

    /** Adds a new callback that will be invoked whenever the value of this parameter changes. 
    This object will take over responsibility to eventually delete the passed object. */
    void addCallback(rsCallbackBase1<void, double> *callbackToAdd);

    /** Deletes all callback objects in our array and clears the array. */
    void deleteAllCallbacks();

    /** Sets the value from a normalized (0...1) value (and updates the real-world value 
    accordingly and invokes the callbacks). */
    void setNormalizedValue(double newNormalizedValue);

    /** Sets the value, updates the normalized value (0...1) accordingly and invokes the 
    callbacks. */
    void setValue(double newValue);

    /** Assigns a new parameter mapper (to/from normalized/actual values). This object will take 
    over responsibility to eventually delete the passed object. */
    void setMapper(rsRealFunctionInvertible *newMapper);


    /** \name Inquiry */

    /** Returns the value normalized to the range 0...1. */
    double getNormalizedValue() const { return normalizedValue; }

    /** Returns the real-wrold value. */
    double getValue() const           { return value; }
  
    /** Returns the name of the parameter. */
    rsString getName() const          { return name;}


  protected:

    /** \name Internal */

    /** Iterates over our callback array, calling each of the callbacks once with our current value
    as argument. */
    void invokeCallbacks();


    /** \name Data */

    double value, normalizedValue, defaultValue;

    rsRealFunctionInvertible *mapper;

    std::vector< rsCallbackBase1<void, double>* > callbacks;

    rsString name, unit;

    DoubleToStringFunctionPointer toString; // function to convert the value to a string

  };

}

#endif
