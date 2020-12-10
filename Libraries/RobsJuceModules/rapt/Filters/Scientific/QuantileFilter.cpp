#undef small

template<class T>
void rsQuantileFilterCore<T>::setLengthAndReadPosition(int newLength, int newPosition, bool hard)
{
  int C = getMaxLength(); // capacity
  rsAssert(newLength   <= C,           "Length cannot exceed capacity");
  rsAssert(newLength   >= 2,           "We require L >= 2");
  rsAssert(newPosition >= 1,           "We require p >= 1");
  rsAssert(newPosition <= newLength-1, "We require p <= L-1");
  if(hard) {
    L = newLength;
    p = newPosition;
    dblHp.setData(&small[0], p, C, &large[0], L-p, C);
    reset();  }
  else
    modulateLengthAndReadPosition(newLength, newPosition);
}

template<class T>
void rsQuantileFilterCore<T>::reset()
{
  for(int n = 0; n < L; n++) {
    dblHp.atIndex(n).value  = T(0);
    int k = dblHp.indexToKey(n);
    dblHp.atIndex(n).bufIdx = n;
    keyBuf[n] = k; }
  bufIdx = 0;
  if(sigBuf)
    sigBuf->reset();
}

/*
template<class T>
T rsQuantileFilterCore<T>::readOutputLongerBy1()
{
  rsAssert(this->sigBuf != nullptr, "To use this feature, the input buffer must be assigned.");
  T xL = (*this->sigBuf)[this->L];   // should be x[n-L], client code must assure this
  return readOutputWithOneMoreInput(xL);
}
// Seems like some member functions are not instantiated when we do the instantiation in the
// rs_testing module because they are never called from there. They are called, however, in the
// test project - but that's apparently too late - the rs_testing module is compiled before that.
// https://stackoverflow.com/questions/31157196/template-member-function-is-instantiated-only-if-called
*/

/*
template<class T>
T rsQuantileFilterCore<T>::readOutputWithOneMoreInput(T x)
{
  int p1;                                                 // read position
  T w1, xS, xL;                                           // weight, xLarge, xSmall
  T q = getQuantile();
  lengthAndQuantileToPositionAndWeight(L+1, q, &p1, &w1);
  T S0 = this->dblHp.getLargestSmallValue().value;        // do we need the "this->"?
  T L0 = this->dblHp.getSmallestLargeValue().value;
  if(p1 == p) {                                           // additional slot is in the large heap
    T S1 = get2ndLargestSmallOrX(x);
    if(     x > L0) { xS = S0; xL = L0; }
    else if(x > S0) { xS = S0; xL = x;  }
    else if(x > S1) { xS = x;  xL = S0; }
    else            { xS = S1; xL = S0; } }
  else {                                                  // additional slot is in the small heap
    rsAssert(p1 == p+1);                                  // sanity check
    T L1 = get2ndSmallestLargeOrX(x);
    if(     x < S0) { xS = S0; xL = L0; }
    else if(x < L0) { xS = x;  xL = L0; }
    else if(x < L1) { xS = L0; xL = x;  }
    else            { xS = L0; xL = L1; } }
  T y = (T(1)-w1)*xS + w1*xL;
  return y;
}
*/
// ToDo: maybe refactor such that a subclass can save time by avoiding calling 
// lengthAndQuantileToPositionAndWeight (it's expensive) when the p1,w1 did not change between 
// samples. They change only when the settings L,w (or L,q) change, so these values could be 
// computed once when the setting is changed and then cached in member variables of a subclass. 
// This can be supported by factoring out a function readOutputWithOneMoreInput(T x, int p1, T w1);
// do we need the "this->"?



/*
template<class T>
void rsQuantileFilterCore<T>::lengthAndQuantileToPositionAndWeight(int L, T q, int* p, T* w)
{
  T P = q * (L-1);            // non-integer read position
  *w  = P - floor(P);         // weight for value right to P
  *p  = floor(P) + 1;         // +1 because p is actually the value to the right of P
  if(*p > L-1) {              // quantile == 1 (maximum) needs special care
    *p = L-1; *w = T(1);  }
}
*/
// todo: optimize - don't call floor twice



template<class T>
void rsQuantileFilterCore<T>::modulateLengthAndReadPosition(int newLength, int newPosition)
{
  int numSmall = newPosition;
  int numLarge = newLength - numSmall;
  while(L > newLength) discardOldestSample();
  while(L < newLength) addOlderSample();
  while(dblHp.getNumSmallValues() < numSmall) moveFirstLargeToSmall();
  while(dblHp.getNumLargeValues() < numLarge) moveFirstSmallToLarge();
  p = dblHp.getNumSmallValues();
}

template<class T>
bool rsQuantileFilterCore<T>::moveFirstLargeToSmall()
{
  if(dblHp.getNumLargeValues() <= 1)
    return false;
  Node n = dblHp.getSmallestLargeValue();   // the node we want to move
  int i = n.bufIdx;                         // index in keyBuf of to-be-moved node
  int k = dblHp.getNumSmallValues();        // the new key for the moved node
  n  = dblHp.large.extractFirst();          // shuffles large heap and buf
  keyBuf[i] = k;                            // set one key value in keyBuf
  dblHp.small.insert(n);                    // shuffles small heap and buf
  return true;
}

template<class T>
bool rsQuantileFilterCore<T>::moveFirstSmallToLarge()
{
  if(dblHp.getNumSmallValues() <= 1)
    return false;
  Node n = dblHp.getLargestSmallValue();
  int i = n.bufIdx;
  int k = dblHp.getNumLargeValues() | firstBitOnly; // set the 1st bit to indicate L-key
  n  = dblHp.small.extractFirst();
  keyBuf[i] = k;
  dblHp.large.insert(n);
  return true;
}

template<class T>
void rsQuantileFilterCore<T>::discardOldestSample()
{
  rsAssert(L > 2); if(L <= 2) return;
  int k = keyBuf[bufIdx];
  dblHp.remove(k);
  bufIdx = wrap(bufIdx+1);              // we need to do a step forward
  L--;
}

template<class T>
void rsQuantileFilterCore<T>::addOlderSample()
{
  T x;
  if(sigBuf) x = sigBuf->fromNewest(L); // use actual old sample, if possible/available
  else       x = readOutput();          // otherwise make up something based on stored values
  bufIdx--;                             // we need to go a step backward
  if(bufIdx < 0)
    bufIdx = (int)keyBuf.size() - 1;    // backward wrap-around at 0
  Node n(x, bufIdx);
  int k = dblHp.getPreliminaryKey(n);   // corresponds to the end of one the heaps
  keyBuf[bufIdx] = k;                   // this may get changed during the actual insert
  k = dblHp.insert(n);                  // lets the new node float up, modifies keyBuf
  L++;
}

template<class T>
bool rsQuantileFilterCore<T>::isStateConsistent()
{
  bool r = true;

  // for converting the loop variable to a buffer index:
  auto convert = [=](int i)->int{ return (i + bufIdx) % (int)keyBuf.size(); };

  // Check that all nodes in the double-heap have the correct back-link to their index in the
  // delay-buffer. This back-link is needed to update the delay-buffer when nodes in the heap are
  // swapped during floatUp/Down:
  for(int i = 0; i < L; i++)
    r &= isBufferSlotConsistent(convert(i));

  // Check that each key occurs in the buf exactly once:
  std::vector<int> tmp(L);
  for(int i = 0; i < L; i++)
    tmp[i] = dblHp.keyToIndex(keyBuf[convert(i)]);
  // r &= isIndexPermutation(&tmp[0], (int)tmp.size()); 
  // commented bcs isIndexPermutation is still in the test-section, not in the library - can be
  // re-activated, when it's dragged over - but this isStateConsistent function was scaffolding 
  // anyway (i.e. only relevant during construction of the class)

  // todo:
  // -check that the sum of the heap sizes matches the buffer size
  // -check that each buffer index occurs in the heap exactly once

  return r;
}

template<class T>
bool rsQuantileFilterCore<T>::isNodeConsistent(const rsQuantileFilterCore<T>::Node& n)
{
  bool r = true;

  for(size_t i = 0; i < this->buf.size(); i++) {
    int  k  = keyBuf[i];      // key of node in i-th buffer slot, use keyBuf.data[i]
    Node n2 = dblHp.atKey(k); // retrieve the node
    if(n2 == n)
      r &= n.bufIdx == i; }

  return r;
}

template<class T>
bool rsQuantileFilterCore<T>::isBufferSlotConsistent(int i)
{
  int k = keyBuf[i];
  Node n = dblHp.atKey(k); // retrieve the node
  bool result = n.bufIdx == i;
  return result;
  // The keyBuf at index i contains a key to a node in the double-heap. That node contains a
  // buffer-index and the buffer index stored in the node should always match the index i. If it
  // doesn't, something is wrong, i.e. the state is inconsistent. This is the way, we map back and
  // forth between heap-nodes sorted by age (in the keyBuf, which is essentially a delayline) and
  // partially sorted by value (in the double-heap). We must be able to figure out, where a node
  // with a given delay is in the heap (to be able to replace the oldest node, when a new sample
  // comes in) and we must also be able to figure out, where a given heap-node's key is in the
  // delay-buffer in order to reshuffle the delay-buffer along, when the heap gets reshuffled due
  // to the replaced node floating up or down.
}

//=================================================================================================

/*
template<class T>
void rsQuantileFilter<T>::convertParameters(
  T length, T quantile, T sampleRate, int* L, int* p, T* w, T* q)
{
  rsAssert(quantile >= T(0) && quantile <= T(1), "Quantile needs to be between 0 and 1");
  *L  = (int) round(length * sampleRate);  // length of filter in samples (maybe use floor later?)
  *L  = rsMax(*L, 2);                      // ...needs to be at least 2

  //*q  = quantile * sampleRate * (*L - 1);  // old - why did this ever work?
  *q  = quantile * (*L - 1);               // readout position in sorted array


  *p  = (int) floor(*q);                   // integer part (floor)
  *w  = *q - *p;                           // fractional part
  *p += 1;                                 // algo wants the next one
  if(*p > *L - 1) {                        // quantile == 1 (maximum) needs special care
    *p = *L - 1; *w = T(1);  }


  *q *= 0.5;                               // found empirically - todo: verify theoretically!
  //*q *= 0.5 * sampleRate;     // needs test
}
// obsolete - can be deleted soon
*/
// It's confusing to use q here - the output *q is actually the delay
//
// part of the code is now implemented in
// rsQuantileFilterCore<T>::lengthAndQuantileToPositionAndWeight, so we should try to get rid of 
// the duplications here

/*

todo: 
-implement a getShortenedOutput similar to getElongatedOutput - this may simplify the 
 implementation of non-integer lengths and even allow for lengths < 2 (down to L=1)
 ...done?

Ideas:

Create a version that generates a pseudo-resonance by:
-running a rsMinMaxFilter in parallel that extracts the min and max of the segment of the same 
 length
-the total output is formed by wq*quantile + wMin*min + wMax*max
-wq and wMin+wMax are obtained by a user crossfade: 0: wq=1, wMin+wMax=0, 1: wq=0, wMin+wMax=1
 this crossfades between normal filter output and "resonance"
-the relative weights of wMin,wMax are derived by some other feature of the signal, for example:
 -the current number of stored samples in the min- and max buffers
 -the position of the current quantile sample in the doubleheap
 -the absolute values of min and max
 -which feature is used can be selected by the user via a "ResoMode" parameter
-the resulting resonance waveform will be some sort of pulse-wave which fits nicely to the general
 character of this type of filter
-we may need to implement an rsMinMaxFilter of non-integer length (by crossfading between two 
 buffered samples, min and 2nd-to-min, ditto for max)

*/