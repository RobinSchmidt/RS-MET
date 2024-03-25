#ifndef RAPT_WINDOWEDFILTERDESIGNER_H_INCLUDED
#define RAPT_WINDOWEDFILTERDESIGNER_H_INCLUDED

class rsWindowedFilterDesigner
{

public:

  template<class T>
  void hilbert(T* h, int numTaps, RAPT::rsWindowFunction::WindowType type);

  template<class T>
  void hilbertSmoothed(T* h, int numTaps, RAPT::rsWindowFunction::WindowType type);

};



#endif