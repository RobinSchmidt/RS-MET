#ifndef RAPT_TYPEDEFINSTANTIATIONS_H
#define RAPT_TYPEDEFINSTANTIATIONS_H

#include "../../../Modules/RAPT.h"
using namespace RAPT;

// Does not yet work

// We create explicit instantiations of the RAPT template classes here and assign a typename to 
// such instances using the following convention: we append one or more letters to template class
// name that indicate the template parameters. For example, PhaseScopeBuffer<float, float, double>
// becomes: PhaseScopeBufferFFD

// Visualization:
typedef RAPT::Image<float> ImageF;
typedef RAPT::AlphaMask<float> AlphaMaskF;
typedef RAPT::ImagePainter<float, float, float> ImagePainterFFF;
typedef RAPT::PhaseScopeBuffer<float, float, double> PhaseScopeBufferFFD;


#endif