#undef small

template<class T>
void rsQuantileFilterCore<T>::setLengthAndReadPosition(int newLength, int newPosition, bool hard)
{
  int C = getMaxLength(); // capacity
  rsAssert(newLength   <= C,           "Length cannot exceed capacity");
  rsAssert(newLength   >= 2,           "We require L >= 2");
  rsAssert(newPosition >= 1,           "We require p <= 1");
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
}

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
  r &= isIndexPermutation(&tmp[0], (int)tmp.size());

  // todo: 
  // -check that the sum of the heap sizes matches the buffer size
  // -check that each buffer index occurs in the heap exactly once

  return r;
}

template<class T>
bool rsQuantileFilterCore<T>::isNodeConsistent(const rsQuantileFilterCore<T>::Node& n)
{
  bool r = true;

  for(size_t i = 0; i < buf.size(); i++) {
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

template<class T>
void rsQuantileFilter<T>::convertParameters(
  T length, T quantile, T sampleRate, int* L, int* p, T* w, T* q)
{
  rsAssert(quantile >= T(0) && quantile <= T(1), "Quantile needs to be between 0 and 1");
  *L  = (int) round(length * sampleRate);  // length of filter in samples
  *L  = rsMax(*L, 2);                      // ...needs to be at least 2
  *q  = quantile * sampleRate * (*L - 1);  // readout position in sorted array
  *p  = (int) floor(*q);                   // integer part (floor)
  *w  = *q - *p;                           // fractional part
  *p += 1;                                 // algo wants the next one
  if(*p > *L - 1) {                        // quantile == 1 (maximum) needs special care
    *p = *L - 1; *w = T(1);  }
  *q *= 0.5;                               // found empirically - todo: verify theoretically!
}