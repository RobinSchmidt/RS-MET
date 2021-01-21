#ifndef rosic_NewSynth_h
#define rosic_NewSynth_h

namespace rosic
{

/** A semimodular syntheyizer with 4 sources that can be mixed via a vector mixer and interact 
through sync, frequency modulation, ring modulation, a dual filter and an arbitrary number of 
freely routable modulators. The sources and filters are actually slots into which various kinds
of source and filter modules can be plugged. 

todo: 
-maybe allow for the same modules to be inserted into source and filter slots - we need no 
 formal distinction - sources can have an input (they may ignore it, if they don't need it) and 
 filters can self-oscillate. 
-also, there should be no distinction between modulators and audio sources at least on the code
 side -> treat all kinds of modules uniformly
-maybe allow source1+2 to be routed into filter1 and source2+3 into filter2
*/

class rsNewSynth : public rsPolyModuleManager // rsPolyModuleObserver
{

public:




  INLINE rsFloat64x2 getSample(const rsFloat64x2& x)
  {
    return 0;
    //return filter.getSample(source.getSample(x));
  }

  INLINE void getSampleFrameStereo(double* left, double* right)
  {
    *left = *right = 0;
    //rsFloat64x2 tmp(*left, *right);
    //*left  = tmp[0];
    //*right = tmp[1];
  }

//protected:

  rsQuadSourcePoly source;
  rsDualFilterPoly filter;
  rsModulatorArrayPoly polyMod;

};


}

#endif