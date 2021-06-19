using namespace rotes;

#include "rosic/rosic.h"
using namespace rosic;

bool rotes::testMinimumPhaseReconstruction()
{
  bool   ok  = true;
  double tol = 1.e-13;  // tolerance

  // Even length case:
  double xe[10] = {1, 2, 1, 0, -2, -1, 1, 0, 2, -1};
  double xem[10];
  minimumPhaseReconstruction(xe, 10, xem);
  ok &= RAPT::rsIsCloseTo(xem[0],  3.2687459271567021,  tol);
  ok &= RAPT::rsIsCloseTo(xem[1],  0.62968035674446710, tol);
  ok &= RAPT::rsIsCloseTo(xem[2],  0.87862532246369252, tol);
  ok &= RAPT::rsIsCloseTo(xem[3], -1.5376097533484183,  tol);
  ok &= RAPT::rsIsCloseTo(xem[4], -1.3261247574468302,  tol);
  ok &= RAPT::rsIsCloseTo(xem[5],  0.45208166454780080, tol);
  ok &= RAPT::rsIsCloseTo(xem[6],  0.32384035922187748, tol);
  ok &= RAPT::rsIsCloseTo(xem[7],  0.77095722011738244, tol);
  ok &= RAPT::rsIsCloseTo(xem[8], -0.14508685139545574, tol);
  ok &= RAPT::rsIsCloseTo(xem[9], -0.31510948806123357, tol);

  // Odd length case:
  double xo[11] = {1, 2, 1, 0, -2, -1, 1, 0, 2, -1, 1};
  double xom[11];
  minimumPhaseReconstruction(xo, 11, xom);
  ok &= RAPT::rsIsCloseTo(xom[0],   2.4491173144733538,    tol);
  ok &= RAPT::rsIsCloseTo(xom[1],   1.0524544833893181,    tol);
  ok &= RAPT::rsIsCloseTo(xom[2],   2.2884160106887808,    tol);
  ok &= RAPT::rsIsCloseTo(xom[3],  -1.1060807731566638,    tol);
  ok &= RAPT::rsIsCloseTo(xom[4],  -1.2870178891607584,    tol);
  ok &= RAPT::rsIsCloseTo(xom[5],  -0.27009549127078175,   tol);
  ok &= RAPT::rsIsCloseTo(xom[6],  -0.49593357262191812,   tol);
  ok &= RAPT::rsIsCloseTo(xom[7],   1.5494684926807107,    tol);
  ok &= RAPT::rsIsCloseTo(xom[8],  -0.0041995463282778469, tol);
  ok &= RAPT::rsIsCloseTo(xom[9],   0.057200299639294377,  tol);
  ok &= RAPT::rsIsCloseTo(xom[10], -0.23332932833305497,   tol);

  return ok;

  // Where are these reference numbers from? matlab/octave? I can't remember.
  // https://de.mathworks.com/help/signal/ref/rceps.html
  //
  // Try:
  // x = [1, 2, 1, 0, -2, -1, 1, 0, 2, -1]
  // [y,yr] = rceps(x)
  // yr
  //
  // If that produces these numbers, document that!
}

