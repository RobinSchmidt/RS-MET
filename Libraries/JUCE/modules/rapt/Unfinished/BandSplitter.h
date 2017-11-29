#ifndef RAPT_BANDSPLITTER_H_INCLUDED
#define RAPT_BANDSPLITTER_H_INCLUDED

/** A filter pair to split an incoming signal into lowpass- and a highpass part. */

template<class TSig, class TPar>
class rsTwoBandSplitter
{

public:

protected:

};

//=================================================================================================

/** Uses an arbitrary number of two-band splitters to split a signal into an arbitrary number of 
bands. */

template<class TSig, class TPar>
class rsMultiBandSplitter
{

public:

protected:

  std::vector<rsTwoBandSplitter<TSig, TPar>*> splitters;

};



#endif