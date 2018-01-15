#ifndef jura_UnitTestView_h
#define jura_UnitTestView_h  

/** A component to perform unit tests for jura classes, print results, etc. */

class JUCE_API UnitTestsView : public Editor
{

public:

  UnitTestsView() {}


protected:

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(UnitTestsView)
};

#endif