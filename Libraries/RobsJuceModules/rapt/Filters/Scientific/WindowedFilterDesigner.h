#ifndef RAPT_WINDOWEDFILTERDESIGNER_H_INCLUDED
#define RAPT_WINDOWEDFILTERDESIGNER_H_INCLUDED

class rsWindowedFilterDesigner
{

public:

  template<class T>
  static void hilbert(T* h, int numTaps, rsWindowFunction::WindowType type);

  template<class T>
  static void hilbertSmoothed(T* h, int numTaps, rsWindowFunction::WindowType type);

};



#endif