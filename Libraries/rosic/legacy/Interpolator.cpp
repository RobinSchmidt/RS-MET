#include "Interpolator.h"

Interpolator::Interpolator()
{
 previousOutput = 0;
}

Interpolator::~Interpolator()
{

}

void Interpolator::reset()
{
 previousOutput = 0;
}