#include "DoubleNumberPair.h"

DoubleNumberPair::DoubleNumberPair()
{
 theNumberPair = _mm_set1_pd(0.0);
}

DoubleNumberPair::DoubleNumberPair(double values)
{
 theNumberPair = _mm_set1_pd(values);
}

DoubleNumberPair::DoubleNumberPair(__m128d newPair)
{
 theNumberPair = newPair;
}

DoubleNumberPair::DoubleNumberPair(double leftValue, double rightValue)
{
 theNumberPair = _mm_set_pd(leftValue, rightValue);
}

DoubleNumberPair::~DoubleNumberPair()
{

}

void DoubleNumberPair::setLeft(double newLeft)
{
 theNumberPair = _mm_set_pd(newLeft, getRight());
}

double DoubleNumberPair::getLeft()
{
 return * (double*) &theNumberPair;
}

void DoubleNumberPair::setRight(double newRight)
{
 theNumberPair = _mm_set_pd(getLeft(), newRight);
}

double DoubleNumberPair::getRight()
{
 return * (((double*) &theNumberPair)+1) ;
}
