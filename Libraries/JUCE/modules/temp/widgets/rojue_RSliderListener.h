#ifndef rojue_RSliderListener_h
#define rojue_RSliderListener_h

namespace rojue
{

  class RSlider;

  class RSliderListener
  {

  public:

    virtual void rSliderValueChanged(RSlider* rSlider) = 0;

  };

}

#endif
