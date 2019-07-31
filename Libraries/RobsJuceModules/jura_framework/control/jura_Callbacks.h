#ifndef jura_Callbacks_h
#define jura_Callbacks_h

//#include <string>  // do we need this? if not - get rid

//=================================================================================================
// class GenericMemberFunctionCallback1:

/** This template can be used to create pointers to callback objects to member-functions of 
as-of-yet-unspecified classes. This is useful, for example, to store an array with pointers to 
callbacks where each array-element may call to an instance of a different class. You can only store
pointer-arrays because this class here is actually abstract so you can't instantiate it. You can, 
however instantiate objects of the subclass SpecificMemberFunctionCallback1 for some specific 
object and assign a baseclass pointer to point to the subclass object. For a single generic 
pointer, it may look like this:

\code
GenericMemberFunctionCallback1<MyReturnType, MyParameterType> *myCallbackPointer =
  new SpecificMemberFunctionCallback1<MyClass, MyReturnType, 
    MyParameterType>(&myInstance, &MyClass::myMemberFunction);
\endcode

When you later invoke the callback via myCallbackPointer->call(myArgument), polymorphism will take 
care of calling the right member function of the concrete object to be called back. */

template<class ReturnType, class ArgumentType>
class JUCE_API GenericMemberFunctionCallback1
{

public:

  /** Virtual destructor */
  virtual ~GenericMemberFunctionCallback1() {}


  /** This is the function to actually invoke the callback. Inside your client-code, the invocation 
  code should look like this:
  \code
  myReturnValue = myCallback.call(myArgument);
  \endcode
  This will then call the assigned member-function on the assigned object with argument 'myArgument' 
  and store the return-value of the invoked member-function in 'myReturnValue'. It is still purely 
  virtual here in this basclass and overriden in the subclass-template 
  SpecificMemberFunctionCallback1 to actually invoke some member-function of some specific class. */
  virtual ReturnType call(const ArgumentType argument) = 0;

  /** Implement function call operator for allowing alternative syntax:
  \code
  myReturnValue = myCallback(myArgument);
  \endcode
  in place of:
  \code
  myReturnValue = myCallback.call(myArgument);
  \endcode  
  It makes it also possible to assign the callback object to a std::function. */
  ReturnType operator()(ArgumentType argument) { return call(argument); }

};

//=================================================================================================
// class SpecificMemberFunctionCallback1:

/** This template can be used to create callback objects that can be assigned to a (non-static) 
member function for some instance of a class that has to be specified via the first template 
argument. This template here is made for member functions that have one parameter and a 
return-value (which may be void). If you have a class 'MyClass' with a member-function 
'myMemberFunction' that takes a parameter of type 'MyParameterType' and returns a value of type 
'MyReturnType', you can instantiate the template like this:

\code
MemberFunctionCallback1<MyClass, MyReturnType, MyParameterType> 
  myCallback(&myInstance, &MyClass::myMemberFunction);
\endcode

where 'myInstance' is the instance (and '&myInstance' therefore is a pointer to it) on which the 
member-function 'myMemberFunction' will be invoked whenever you call myCallback.call(myArgument).

Remark: the '()' operator is deliberately not implemented so you are forced to use the '.call()' 
syntax - mainly because i believe that callback invocations should be distinguished from regular 
function calls in client code to avoid confusion. */

template<class CalleeObjectType, class ReturnType, class ArgumentType>
class JUCE_API SpecificMemberFunctionCallback1  // get rid of "Specific"
  : public GenericMemberFunctionCallback1<ReturnType, ArgumentType>
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction:

  /** Constructor. syntax: 
  MemberFunctionCallback1<MyReturnType, MyParameterType, MyClass> 
    myCallback(&myInstance, &MyClass::myMemberFunction); */
  SpecificMemberFunctionCallback1(CalleeObjectType *calleeObject, 
    ReturnType (CalleeObjectType::*memberToCall) (ArgumentType))
  {
    this->calleeObject = calleeObject;
    this->memberToCall = memberToCall;
  }

  //-----------------------------------------------------------------------------------------------
  // callback invocation:

  /** Overriden from GenericMemberFunctionCallback1::call to actually invoke the desired 
  member-function on the callee object. */
  ReturnType call(const ArgumentType argument)
  {
    return (*calleeObject.*memberToCall)(argument);
  }

  //-----------------------------------------------------------------------------------------------
  // operators:

  /** Compares two Callbacks for equality - two callbacks are considered equal when they point to 
  the same member-function and to the same instance. */
  bool operator==(const SpecificMemberFunctionCallback1& other) const
  {
    if(calleeObject == other.calleeObject && memberToCall == other.memberToCall)
      return true;
    else
      return false;
  }

  /** Compares two Callbacks for inequality - as usual, this is just defined as the negation of 
  equality. */
  bool operator!=(const SpecificMemberFunctionCallback1& other) const
  {
    return !(*this == other);
  }

protected:

  ReturnType (CalleeObjectType::*memberToCall) (ArgumentType);
  CalleeObjectType *calleeObject;

};

#endif 
