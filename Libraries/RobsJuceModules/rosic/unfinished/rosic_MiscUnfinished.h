#pragma once

namespace rosic
{


//=================================================================================================

template<class TSig, class TPar> // move to RAPT
class rsVectorMixer
{

public:

  enum class Mode
  {
    linear,
    sinCosApprox,
    sinCos
  };

  void setX(TPar newX) { x = newX; }

  void setY(TPar newY) { y = newY; }

  void getGains(TSig* topLeft, TSig* topRight, TSig* bottomLeft, TSig* bottomRight)
  {
    TPar xx = TPar(0.5) * (x + TPar(1));  // map -1..+1 to 0..1
    TPar yy = TPar(0.5) * (y + TPar(1));

    // todo: switch between modes - this is for linear:
    TPar right  = xx;
    TPar left   = TPar(1) - xx;
    TPar top    = yy;
    TPar bottom = TPar(1) - yy;

    *topLeft     = TSig(top    * left);
    *topRight    = TSig(top    * right);
    *bottomLeft  = TSig(bottom * left);
    *bottomRight = TSig(bottom * right);
  }

protected:

  Mode mode = Mode::linear;

  TPar x = TPar(0), y = TPar(0);

};
// todo: maybe have separate coordinates for left and right channel xL, yL, xR, yR
// may de-templatize it -> have xL, yL, xR, yR as double and the 4 gains are computed as 
// rsFloat64x2


class rsVectorMixerPoly : public rsPolyModule, public rsVectorMixer<rsFloat64x2, double>
{

public:




  /*
  void processFrame(const double* in, int numIns, double* out, int numOuts, int voice) override
  {
    RAPT::rsAssert(numOuts == 4);
    out[0] = out[1] = out[2] = out[3] = 1.0; // preliminary
  }
  */

protected:

};

//=================================================================================================

class rsModulatorArrayPoly : public rsPolyModule
{

public:

protected:

};

//=================================================================================================

class rsTriSawOscPoly : public rsPolyModule, public RAPT::rsTriSawOscillator<rsFloat64x2>
{

public:


  rsFloat64x2 getSample(const rsFloat64x2& in, int voice) override
  {
    return rsFloat64x2(0, 0);
  }


protected:

  std::vector<RAPT::rsTriSawVoice<rsFloat64x2>> voices;
  
  //RAPT::rsTriSawVoice<rsFloat64x2> voices; // the code for this is in the wrong file - move over

};

//=================================================================================================

class rsLadderFilterPoly : public rsPolyModule
{

public:

protected:

};

//=================================================================================================

class rsAttackDecayEnvelopePoly : public rsPolyModule
{

public:

protected:

};







}