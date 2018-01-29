#ifndef rosic_PolyModule_h
#define rosic_PolyModule_h

namespace rosic
{




class rsPolyModule
{

public:

  virtual rsFloat64x2 getSample(rsFloat64x2 in) = 0;

protected:


};


}

#endif