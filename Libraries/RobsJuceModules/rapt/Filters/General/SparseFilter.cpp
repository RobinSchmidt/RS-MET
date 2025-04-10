

/*

ToDo:

- Move the code for rsSparseDigitalTransferFunction into a dedicated pair of .h/.cpp files. I moved
  it temporarily over from SparseFilter.h/cpp into the SparseRationalTransferFunction.h/cpp pair 
  because the class needs to be defined before DelayLine.h is included because clas rsDelay has a 
  function that returns such a transfer function and on the other hand, the class rsSparseFiler 
  needs the definition of rsDelay. Moving it into SparseRationalTransferFunction solves these 
  problems but conceptually, the class feels a bit misplaced in that file. It needs its own file
  to sove the include order problem and have it placed into a file that makes sense. Then include
  SparseDigitalTransferFunction.h before DelayLine.h before SparseFilter.h


*/