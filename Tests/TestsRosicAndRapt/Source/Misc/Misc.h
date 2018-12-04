#pragma once

#include "../rosic_tests/PortedFromRSLib/ExamplesRSLib.h"


void sampleTailExtenderTest();



class MemLeakTest
{
  //std::vector<ModuleTypeInfo*> typeInfos; // having this as member causes the memleak
  //std::vector<int> test; // this also
  // so it seems that when i have a global object of a class that has a std::vector member,
  // i'll get a memleak


  //std::vector<int> testVector;         // triggers pre-exit-of-main memleak detection
  std::vector<int>* testPointerToVector; // doesn't trigger it

  // there is a way to prevent the debugger to trigger false memory leaks, explained here:
  // http://www.cplusplus.com/forum/general/9614/ it says:
  // So the way to get around the problem is this:
  // struct MemCheck {
  //   ~MemCheck() {
  //     _CrtDumpMemoryLeaks();
  //   }
  // };
  // And then instantiate MemCheck as the first variable in the block/function so that it gets 
  // destroyed last, after everything else.

  // so, using a class/struct that calls _CrtDumpMemoryLeaks(); in its destructor and then 
  // instantiating a variable of that class, we ensure that the destructor of that variable is 
  // called *after* the exit point of our main function
  // maybe have two checks pre-exit memory-leaks and post-exit memory leaks
  // maybe name the struct PostExitMemLeakChecker

};
extern MemLeakTest memLeakTest;
// ok - something like this has been added to romos.h/cpp needs to be cleaned up - move this to 
// somewhere else - maybe make a test-project solely for demonstrating the behavior of the memleak 
// checker - how it triggers false positives and how to avoid it