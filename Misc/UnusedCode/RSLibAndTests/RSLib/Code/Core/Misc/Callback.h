#ifndef RS_CALLBACK_H
#define RS_CALLBACK_H

namespace RSLib
{

  /**

  This template can be used to create pointers to callback objects to member-functions of
  as-of-yet-unspecified classes. This is useful, for example, to store an array with pointers to
  callbacks where each array-element may call to an instance of a different class. You can only
  store pointer-arrays because this class here is actually abstract so you can't instantiate it.
  You can, however instantiate objects of the subclass rsCallback1 for some specific object and
  assign a baseclass pointer to point to the subclass object. For a single generic pointer, it may
  look like this:

  \code
  rsCallbackBase1<MyReturnType, MyParameterType> *myCallbackPointer =
    new rsCallback1<MyClass, MyReturnType, MyParameterType>
    (&myInstance, &MyClass::myMemberFunction);
  \endcode

  When you later invoke the callback via myCallbackPointer->call(myArgument), polymorphism will
  take care of calling the right member function of the concrete object to be called back.

  */

  template<class ReturnType, class ArgumentType>
  class rsCallbackBase1
  {

  public:

    /** An abstract class should have a virtual destructor. */
    //virtual ~rsCallbackBase1() {}
      // leads to segmentation fault error in the "testParameter" unit test with gcc

    /** This is the function to actually invoke the callback. Inside your client-code, the
    invocation code should look like this:
    \code
    myReturnValue = myCallback.call(myArgument);
    \endcode
    This will then call the assigned member-function on the assigned object with argument
    'myArgument' and store the return-value of the invoked member-function in 'myReturnValue'.
    It is still purely virtual here in this basclass and overriden in the subclass-template
    Callback1 to actually invoke some member-function of some specific class. */
    virtual ReturnType call(const ArgumentType argument) = 0;

  };

  //===============================================================================================
  // class rsCallback1:

  /**

  This template can be used to create callback objects that can be assigned to a (non-static)
  member function for some instance of a class that has has to specified via the first template
  arguments. This template here is made for member functions that have one parameter and a
  return-value (which may be void). If you have a class 'MyClass' with a member-function
  'myMemberFunction' that takes a parameter of type 'MyParameterType' and returns a value of type
  'MyReturnType', you can instantiate the template like this:

  \code
  rsCallback1<MyClass, MyReturnType, MyParameterType>
    myCallback(&myInstance, &MyClass::myMemberFunction);
  \endcode

  where 'myInstance' is the instance (and '&myInstance' therefore is a pointer to it) on which the
  member-function 'myMemberFunction' will be invoked whenever you call myCallback.call(myArgument).

  Remark: the '()' operator is deliberately not implemented so you are forced to use the '.call()'
  syntax - mainly because i believe that callback invocations should be distinguished from regular
  function calls in client code to avoid confusion.

  */

  template<class CalleeObjectType, class ReturnType, class ArgumentType>
  class rsCallback1 : public rsCallbackBase1<ReturnType, ArgumentType>
  {

  public:

    /** Constructor. Syntax:
    rsCallback1<MyReturnType, MyParameterType, MyClass>
      myCallback(&myInstance, &MyClass::myMemberFunction);  */
    rsCallback1(CalleeObjectType *calleeObject,
      ReturnType (CalleeObjectType::*memberToCall) (ArgumentType))
    {
      this->calleeObject = calleeObject;
      this->memberToCall = memberToCall;
    }

    /** Destructor.  */
    //virtual ~rsCallback1() {}

    /** Overriden from rsGenericMemberCallback1::call to actually invoke the desired
    member-function on the callee object. */
    ReturnType call(const ArgumentType argument)
    {
      return (*calleeObject.*memberToCall)(argument);
    }

    /** Compares two callback objacts for equality - two callbacks are considered equal when they
    point to the same member-function and to the same instance. */
    bool operator==(const rsCallback1& other) const
    {
      if( calleeObject == other.calleeObject && memberToCall == other.memberToCall )
        return true;
      else
        return false;
    }

    /** Compares two callback objects for inequality - as usual, this is just defined as the
    negation of equality. */
    bool operator!=(const rsCallback1& other) const
    {
      return !(*this == other);
    }

  protected:

    ReturnType (CalleeObjectType::*memberToCall) (ArgumentType);
    CalleeObjectType *calleeObject;

  };


  /**

  \todo create boiler-plate versions for callbacks with other numbers of parameters (0,2,3,4, etc.)
  when needed - though, as soon as C++ supports variadic templates, they should be replaced by
  that

  \todo optimizations:
  try the following variants:
  -replace the purely virtual 'call' function in the baseclass with a function pointer
   ("manual polymorphism - see paper by Jakubik)
  -maybe try __fastcall calling convention
  -maybe use void-pointers in the baseclass and 'type-destroyer is type-restorer' (see paper by
   Jakubik) - may replace virtual declaration and even allow for enhanced by-value semantics
   (comparability of generic callbacks)

  */

}

#endif
