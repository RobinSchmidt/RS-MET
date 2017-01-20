#ifndef RAPT_ALPHAMASK_H_INCLUDED
#define RAPT_ALPHAMASK_H_INCLUDED

/** This is a subclass of Image intended to be used as a kind of prototype dot.

Maybe rename to ImageMask

*/

template<class TPix>  // pixel type
class AlphaMask : public ImageResizable<TPix>
//class AlphaMask : private Image<TPix>
{

public:


  /** \name Setup */

  //void setMaxPixelSize(int newMaxWidth, int newMaxHeight);

  void setSize(double newSize);

  void setShape(int newShape);
    // circle, rectangle, etc


protected:


};

#endif