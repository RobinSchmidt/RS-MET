#ifndef romos_GenerateDesiredOutput_h
#define romos_GenerateDesiredOutput_h

//#include "../Framework/romos_Module.h"
//#include "../Framework/romos_NoteEvent.h"
//#include "../Algorithms/romos_FilterDesign.h"

namespace rsTestRomos
{

/** This class can generate the desired correct output signals that the romos-modules should match.

Why do i have a class for this? Wouldn't it be better to have the code directly in the respective 
unit test? */

class GenerateDesiredOutput
{

public:

  // make these work also for non-distinct in/out buffers...

  static void forImpulse(int N, double  *x);

  static void forWhiteNoiseUniform(int N, double  *d, unsigned long seed);

  static void forUnitDelay(int N, double  *x, double  *d);
  static void forSummedDiffs(int N, double **x, double **d);
  static void forMovingAverage(int N, double  *x, double *b0, double *b1, double *d);
  static void forLeakyIntegrator(int N, double  *x, double *c, double *d);
  static void forLeakyIntegratorDoubleDelay(int N, double  *x, double *c, double *d);
  static void forBiquad(int N, double  *x, double *b0, double *b1, double *b2, double *a1, double *a2, double *d);
  static void forBiquadWithFixedCoeffs(int N, double  *x, double  b0, double  b1, double  b2, double  a1, double  a2, double *d);


  static void forFormula1In1Out(int N, double *x, double *d);

  static void forTestFilter1(int N, double  *x, double *b0, double *b1, double *c,
    double *dSum, double *dDiff, double *dProd);




  static void forFilterBlip(int N, double frequency, double q, double *desiredOutput);

  static void forGatedNoteFrequencies(int N, std::vector<romos::NoteEvent> *events, double ***d,
    bool containerIsPolyphonic, bool noteFreqModuleIsPolyphonic);

};

}

#endif 
