// ToDo: Figure out if we still need this file anywhere. If not, delete it. If yes, try to change
// the places where we need it such that it's not needed there anymore - and then delete it.
// Or: Document why this file is still here. Maybe it should serve as a reminder for how a
// library include file could look like if we don't need to adhere to the JUCEs module definition
// format? That may make sense. But it should then be documented.

/**

This is the main include file for the RoSiC (Robin's Signal Processing 
Classes) library.

*/

#ifndef rosic_h
#define rosic_h

#include "analysis/rosic_Analysis.h"
#include "basics/rosic_Basics.h"
#include "datastructures/rosic_DataStructures.h"
#include "delaylines/rosic_DelayLines.h"
#include "dynamics/rosic_Dynamics.h"
#include "effects/rosic_Effects.h"
#include "filters/rosic_Filters.h"
#include "generators/rosic_Generators.h"
#include "infrastructure/rosic_Infrastructure.h"
#include "instruments/rosic_Instruments.h"
#include "math/rosic_Math.h"
#include "modulators/rosic_Modulators.h"
#include "neural/rosic_Neural.h"
#include "numerical/rosic_Numerical.h"
#include "others/rosic_Others.h"
#include "rendering/rosic_Rendering.h"
#include "scripting/rosic_Scripting.h"
//#include "plugins/rosic_PlugIns.h"
#include "transforms/rosic_Transforms.h"

#endif // end of #ifndef rosic_h







