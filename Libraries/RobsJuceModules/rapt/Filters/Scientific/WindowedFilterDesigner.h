#ifndef RAPT_WINDOWEDFILTERDESIGNER_H_INCLUDED
#define RAPT_WINDOWEDFILTERDESIGNER_H_INCLUDED


/** This class has functions to design various filters using the windowing method. The impulse 
response will be written into the array "h" which must be "numTaps" long. You also need to 
specify which window function shall be used via the "type" parameter as one of the types 
defined in the  rsWindowFunction::WindowType enum class. */

class rsWindowedFilterDesigner
{

public:


  template<class T>
  static void hilbert(T* h, int numTaps, rsWindowFunction::WindowType type);

  template<class T>
  static void hilbertSmoothed(T* h, int numTaps, rsWindowFunction::WindowType type, 
    bool evenNominalLength);




};



#endif