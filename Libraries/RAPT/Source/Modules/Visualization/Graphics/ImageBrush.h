#ifndef RAPT_IMAGEBRUSH_H_INCLUDED
#define RAPT_IMAGEBRUSH_H_INCLUDED

/** This is a subclass of Image intended to be used as a kind of prototype dot */

template<class TPix, class TPar>  // pixel, parameter types
class ImageBrush : public Image<TPix>
{

public:


  /** \name Setup */

  void setSize(TPar newSize);

  void setShape(int newShape);
    // circle, rectangle, etc


protected:


};

#endif