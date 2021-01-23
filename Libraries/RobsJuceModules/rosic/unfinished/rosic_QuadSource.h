#ifndef rosic_QuadSource_h
#define rosic_QuadSource_h

namespace rosic
{

/** A class that combines four sources and lets their amplitude be adjusted by a vector mixer. 
However, for more flexibility in subsequent processing stages, the outputs of the 4 source are 
not mixed internally. The sources themselves are assumed to be stereo, so this module produces
8 output signals in total (as 4 stereo pairs). */

class rsQuadSourcePoly : public rsPolyModule
{

public:

  rsQuadSourcePoly();

  virtual ~rsQuadSourcePoly();

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Replaces the module in the given slot with the passed new module and takes ownership over
  the newModule object. */
  void setModule(int slot, rsPolyModule* newModule)
  {
    RAPT::rsAssert(slot >= 0 && slot <= numSources, "invalid slot index");
    RAPT::rsAssert(newModule != nullptr, "must be a valid object, or else access violation");
    delete sources[slot];
    sources[slot] = newModule;
  }


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry


  /*
  rsPolyModule* getModulePointer(int slot)
  {

  }
  */


  //-----------------------------------------------------------------------------------------------
  // \name Processing


  rsFloat64x2 getSample(const rsFloat64x2& in, int voice) override
  { 
    rsFloat64x2 gains[4];
    mixer->getGains(&gains[0], &gains[1], &gains[2], &gains[3]);

    rsFloat64x2 out(0, 0);
    for(int i = 0; i < numSources; i++)
      out += gains[i] * sources[i]->getSample(in, voice);
    return out;
  }


  /*
  void processFrameNoMix(const double* in, int numIns, double* out, int numOuts, int voice)
  {
    RAPT::rsAssert(numOuts == 8); // 4 stereo outs

    // Each source is supposed to produce a stereo pair of outputs:
    sources[0]->processFrame(in, numIns, &out[0], 2, voice);
    sources[1]->processFrame(in, numIns, &out[2], 2, voice);
    sources[2]->processFrame(in, numIns, &out[4], 2, voice);
    sources[3]->processFrame(in, numIns, &out[6], 2, voice);

    // Apply the gains of the vector mixer, but we don't mix them internally to let the subsequent 
    // DualFilter be more flexible with the routing:
    double a[4];
    mixer->processFrame(in, numIns, a, 4, voice);
    out[0] *= a[0];
    out[1] *= a[0];
    out[2] *= a[1];
    out[3] *= a[1];
    out[4] *= a[2];
    out[5] *= a[2];
    out[6] *= a[3];
    out[7] *= a[3];
    // maybe factor out the weighting by the gain factors
  }

  void processFrame(const double* in, int numIns, double* out, int numOuts, int voice) override
  {
    RAPT::rsAssert(numOuts == 2); // stereo out
    double tmp[8];
    processFrameNoMix(in, numIns, tmp, 8, voice);
    out[0] = tmp[0] + tmp[2] + tmp[4] + tmp[6];
    out[1] = tmp[1] + tmp[3] + tmp[5] + tmp[7];
  }
  */

  // todo: implement processBlock


protected:

  static const int numSources = 4;

  rsPolyModule*      sources[numSources];  // the sources
  rsVectorMixerPoly* mixer;                // maybe use direct member (no pointer)
  // maybe use smart pointers - which? std::shared_ptr? std::unique_ptr?
  // the sources and the mixer are owned here (is this a good idea)

};


}

#endif