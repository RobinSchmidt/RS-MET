template <class T>
bool defaultLess(const T& left, const T& right)
{
  return left < right;
}

template <class T>
void rsMaxHeapify(T *buffer, int length, int i, int heapSize, 
  bool (*less)(const T& left, const T& right))
{
  int left  = 2*i+1;
  int right = 2*i+2;
  int largest;
  if(left <= heapSize-1 && less(buffer[i], buffer[left]))
    largest = left;
  else
    largest = i;

  if(right <= heapSize-1 && less(buffer[largest], buffer[right]))
    largest = right;

  if(largest != i) {
    rsSwap(buffer[i], buffer[largest]);
    rsMaxHeapify(buffer, length, largest, heapSize, less); 
  }
}
// todo: try to get rid of the recursive call to itself -> convert to iteration

template <class T>
void rsBuildMaxHeap(T *buffer, int length, 
  bool (*less)(const T& left, const T& right))
{
  for(int i = length/2-1; i >= 0; i--)
    rsMaxHeapify(buffer, length, i, length, less);
}

template <class T>
void rsHeapSort(T *buffer, int length, bool (*less)(const T& left, const T& right))
{
  rsBuildMaxHeap(buffer, length, less);
  int heapSize = length;
  for(int i = length-1; i >= 1; i--) {
    rsSwap(buffer[0], buffer[i]);
    heapSize--;
    rsMaxHeapify(buffer, length, 0, heapSize, less);
  }
}

/*
template <class T>
bool rsIsSortedAscending(T *buffer, int length)
{
  for(int i = 0; i < length-1; i++) {
    if(!(buffer[i] <= buffer[i+1]))
      return false;
  }
  return true;
}

template <class T>
bool rsIsSortedStrictlyAscending(T *buffer, int length)
{
  for(int i = 0; i < length-1; i++) {
    if(!(buffer[i] < buffer[i+1]))
      return false;
  }
  return true;
}
*/
// moved to rsArrayTools

template <class T>
std::vector<int> rsFindAllOccurencesOf(T* buffer, int bufferLength,
  T* pattern, int patternLength)
{
  std::vector<int> table(patternLength + 1, -1);
  for(int i = 1; i <= patternLength; i++)
  {
    int pos = table[i-1];
    while(pos != -1 && pattern[pos] != pattern[i-1])
      pos = table[pos];
    table[i] = pos + 1;
  }

  std::vector<int> matches;
  int bufferIndex  = 0;
  int patternIndex = 0;
  while(bufferIndex < bufferLength)
  {
    while(patternIndex != -1 && (patternIndex == patternLength
      || pattern[patternIndex] != buffer[bufferIndex]))
      patternIndex = table[patternIndex];
    patternIndex++;
    bufferIndex++;
    if(patternIndex == patternLength)
      matches.push_back(bufferIndex - patternLength);
  }

  return matches;
}

template <class T>
int rsFindFirstOccurrenceOf(const T *buffer, const int length,
  const T &elementToFind, const int searchStart)
{
  for(int i = searchStart; i < length; i++)
  {
    if(buffer[i] == elementToFind)
      return i;
  }
  return -1;
}

template <class T>
int rsFindFirstOccurrenceOf(const T *buffer, const int length, const T *patternToMatch,
  const int patternLength, const int searchStart)
{
  for(int shift = searchStart; shift <= length-patternLength; shift++)
  {
    if(rsArrayTools::equal(&buffer[shift], patternToMatch, patternLength))
      return shift;
  }
  return -1;
}

template<class T>
int rsFindHighestPeakIndex(T *buffer, int length)
{
  int maxIndex = 0;
  T maxValue   = -RS_INF(T);
  T value;
  for(int i = 1; i < length-1; i++)
  {
    value = buffer[i];
    if(value > buffer[i-1] && value >= buffer[i+1] && value >= maxValue)
    {
      maxValue = value;
      maxIndex = i;
    }
  }
  return maxIndex;
}

template <class T>
int rsFindLastOccurrenceOf(const T *buffer, const int length, const T &elementToFind,
  const int searchStart)
{
  int i = searchStart;
  if(i < 0)
    i = length-1;
  while(i >= 0)
  {
    if(buffer[i] == elementToFind)
      return i;
    i--;
  }
  return -1;
}

template<class T>
T rsFindNearestUpwardZeroCrossing(T* buffer, int length, T searchStart)
{
  // initialize the search pointers for above and below the search start
  int hi1 = (int)ceil(searchStart);
  int hi2 = hi1 + 1;
  int lo1 = (int)floor(searchStart);
  int lo2 = lo1 - 1;

  // catch some special conditions:
  if(lo2 < 0 || hi2 >= length)
    return searchStart;

  // check whether the search start itself is in between two values which surround a zero
  // crossing:
  if(buffer[lo1] <= 0.0 && buffer[hi1] > 0.0)
  {
    return lo1 + buffer[lo1] / (buffer[lo1]-buffer[hi1]);
  }

  // search to the right (ascending direction):
  while(hi2 < length)
  {
    if(buffer[hi1] <= 0.0 && buffer[hi2] > 0.0)
      break;
    hi1++;
    hi2++;
  }

  // search to the left (descending direction):
  while(lo2 >= 0)
  {
    if(buffer[lo2] <= 0.0 && buffer[lo1] > 0.0)
      break;
    lo1--;
    lo2--;
  }

  // \todo optimization: search to both directions simultaneously and stop when a zero crossing
  // is found in either one of the directions

  // if the search was without success, we return the starting point of the search:
  if(lo2 == 0 && hi2 == length-1)
    return searchStart;

  // select the closer from both and return the zero crossing of the connecting line:
  if((searchStart-lo2) <= -(searchStart-hi1))
    return lo2 + buffer[lo2] / (buffer[lo2]-buffer[lo1]);
  else
    return hi1 + buffer[hi1] / (buffer[hi1]-buffer[hi2]);
}