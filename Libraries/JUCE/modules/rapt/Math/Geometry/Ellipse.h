#ifndef RAPT_ELLIPSE_H_INCLUDED
#define RAPT_ELLIPSE_H_INCLUDED

/** This is a class for dealing with ellipses. */

template<class T>
class Ellipse : public ConicSection<T>
{

public:


protected:

  T scale = 1, ratio = 1, angle = 0, centerX = 0, centerY = 0;

};

#endif
