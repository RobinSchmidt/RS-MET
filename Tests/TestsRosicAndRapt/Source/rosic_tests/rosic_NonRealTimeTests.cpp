using namespace rotes;

#include "rosic/rosic.h"
using namespace rosic;

bool rotes::testMinimumPhaseReconstruction()
{
  bool   result = true;
  double tol    = 1.e-14;  // tolerance

  // even length case:
  double xe[10] = {1, 2, 1, 0, -2, -1, 1, 0, 2, -1};
  double xem[10];
  minimumPhaseReconstruction(xe, 10, xem);
  result &= RAPT::rsIsCloseTo(xem[0],  3.2687459271567021,  tol);
  result &= RAPT::rsIsCloseTo(xem[1],  0.62968035674446710, tol);
  result &= RAPT::rsIsCloseTo(xem[2],  0.87862532246369252, tol);
  result &= RAPT::rsIsCloseTo(xem[3], -1.5376097533484183,  tol);
  result &= RAPT::rsIsCloseTo(xem[4], -1.3261247574468302,  tol);
  result &= RAPT::rsIsCloseTo(xem[5],  0.45208166454780080, tol);
  result &= RAPT::rsIsCloseTo(xem[6],  0.32384035922187748, tol);
  result &= RAPT::rsIsCloseTo(xem[7],  0.77095722011738244, tol);
  result &= RAPT::rsIsCloseTo(xem[8], -0.14508685139545574, tol);
  result &= RAPT::rsIsCloseTo(xem[9], -0.31510948806123357, tol);
  // Where are these numbers from? matlab/octave?


  // odd length case:
  double xo[11] = {1, 2, 1, 0, -2, -1, 1, 0, 2, -1, 1};
  double xom[11];
  minimumPhaseReconstruction(xo, 11, xom);
  result &= RAPT::rsIsCloseTo(xom[0],   2.4491173144733538,    tol);
  result &= RAPT::rsIsCloseTo(xom[1],   1.0524544833893181,    tol);
  result &= RAPT::rsIsCloseTo(xom[2],   2.2884160106887808,    tol);
  result &= RAPT::rsIsCloseTo(xom[3],  -1.1060807731566638,    tol);
  result &= RAPT::rsIsCloseTo(xom[4],  -1.2870178891607584,    tol);
  result &= RAPT::rsIsCloseTo(xom[5],  -0.27009549127078175,   tol);
  result &= RAPT::rsIsCloseTo(xom[6],  -0.49593357262191812,   tol);
  result &= RAPT::rsIsCloseTo(xom[7],   1.5494684926807107,    tol);
  result &= RAPT::rsIsCloseTo(xom[8],  -0.0041995463282778469, tol);
  result &= RAPT::rsIsCloseTo(xom[9],   0.057200299639294377,  tol);
  result &= RAPT::rsIsCloseTo(xom[10], -0.23332932833305497,   tol);

  return result;
  // fails! maybe we need higher tolerance?
}

