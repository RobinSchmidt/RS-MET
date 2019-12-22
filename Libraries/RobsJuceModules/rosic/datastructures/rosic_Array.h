#ifndef rosic_Array_h
#define rosic_Array_h

namespace rosic
{

/** This class implements a simple array which can dynamically grow and shrink as needed.  */

// rename to rsDynamicArray to avoid name-clash with RAPT::rsArray - and maybe move to rapt

template<class ElementType>
class rsArray
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. If initialSize is nonzero, memory will be reserved to hold just that number of
  elements and the array will filled with default elements as they are constructed by the default
  constructor of the element class. */
  rsArray(const int initialSize = 0)
  {
    elements     = NULL;
    numUsed      = 0;
    numAllocated = 0;
    //granularity  = 8;
    if(initialSize > 0)
    {
      elements     = new ElementType[initialSize];
      numAllocated = initialSize;
      numUsed      = initialSize;
    }
  }

  /** Copy-constructor - creates a deep copy of the other array. */
  rsArray(const rsArray<ElementType>& other)
  {
    elements     = new ElementType[other.numUsed];
    numUsed      = other.numUsed;
    numAllocated = other.numUsed;
    //granularity  = other.granularity;
    for(int i=0; i<numUsed; i++)
      elements[i] = other.elements[i];
  }

  /** Destructor. */
  ~rsArray()
  {
    if(elements != NULL)
      delete[] elements;
  }

  //---------------------------------------------------------------------------------------------
  // element access:


  /** Copies the content of the array into a std::vector (from the STL). */
  void toVectorSTL(std::vector<ElementType> &vectorToCopyInto)
  {
    vectorToCopyInto.clear();
    vectorToCopyInto.reserve(numUsed);
    for(int i=0; i<numUsed; i++)
      vectorToCopyInto.push_back(elements[i]);
  }

  /** Copies the content of a std::vector (from the STL). */
  void fromVectorSTL(std::vector<ElementType> &vectorToCopyFrom)
  {
    clear();
    ensureAllocatedSize(vectorToCopyFrom.size());
    for(int i=0; i<vectorToCopyFrom.size(); i++)
      elements[i] = vectorToCopyFrom[i];
    numUsed = vectorToCopyFrom.size();
  }

  /** Adds an element to this array. */
  void appendElement(const ElementType newElement)
  {
    ensureAllocatedSize(numUsed+1);
    elements[numUsed] = newElement;
    numUsed++;
  }

  /** Adds an element to this array when the element in question odes not already exist (via the == operator). */
  void appendIfNotAlreadyThere(const ElementType newElement)
  {
    if(!hasElement(newElement))
      appendElement(newElement);
  }

  /** Appends another array to this one. */
  void appendArray(const rsArray<ElementType> &arrayToAppend)
  {
    ensureAllocatedSize(numUsed + arrayToAppend.numUsed);
    for(int i=0; i<arrayToAppend.numUsed; i++)
      elements[numUsed+i] = arrayToAppend.elements[i];
    numUsed += arrayToAppend.numUsed;
  }

  /** Searches the array for the element (using the == operator) and returns the index or -1
  if the element is not found. */
  int findElement(const ElementType &elementToFind) const
  {
    for(int i=0; i<numUsed; i++)
    {
      if(elements[i] == elementToFind)
        return i;
    }
    return -1;
  }

  /** Returns true, if the array contains the element, false otherwise. */
  bool hasElement(const ElementType elementToFind) const
  {
    return findElement(elementToFind) != -1;
  }

  /** Inserts an element to a given position. All elements after this position will be moved one
  position to the right to make up space. */
  void insertElement(const int index, const ElementType newElement)
  {
    if(index >= 0 && index < numUsed)
    {
      ensureAllocatedSize(numUsed+1);
      for(int i=numUsed-1; i>=index; i--)
        elements[i+1] = elements[i];
      elements[index] = newElement;
      numUsed++;
    }
    else if(index >= numUsed)
      appendElement(newElement);
    else
      DEBUG_BREAK;  // negative array index
  }

  /** Removes an element from a given position. All elements after this position will be moved
  one position to the left to fill the gap. */
  void removeElementByIndex(const int index)
  {
    if(index >= 0 && index < numUsed)
    {
      for(int i=index; i<numUsed-1; i++)
        elements[i] = elements[i+1];
      numUsed--;
    }
    else
    {
      DEBUG_BREAK; // invalid array index
    }
  }

  /** Removes an element from the array, if present. */
  void removeElementByValue(const ElementType elementToRemove)
  {
    int index = findElement(elementToRemove);
    if(index != -1)
      removeElementByIndex(index);
  }

  /** Clears the array. */
  void clear()
  {
    numUsed = 0;
    // optionally free  memory here
  }

  /** Replaces an element with a new one - if there is no element with the given index, nothing
  will be done. */
  void replaceElement(const int index, const ElementType newElement)
  {
    if(index >= 0 && index < numUsed)
      elements[index] = newElement;
    else
      DEBUG_BREAK; // invalid array index
  }

  /** Returns the element at a given index - this functions does bounds-checking and will return
  the neutral element of the ElementType class (as created by its standard-constructor). */
  ElementType getElement(const int index) const
  {
    if(index >= 0 && index < numUsed)
      return elements[index];
    else
    {
      DEBUG_BREAK;         // invalid array index
      return ElementType();
    }
  }

  //---------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the number of elements in the array. */
  INLINE int getNumElements() const { return numUsed; }

  /** Returns true if this array is empty, false otherwise. */
  INLINE bool isEmpty() const { return getNumElements() == 0; }

  //---------------------------------------------------------------------------------------------
  // array operators:

  /** Allows fast indexed acces to array elements - warning no bounds checking is done. */
  inline ElementType& operator[] (const int index)
  {
    return elements[index];
  }

  inline const ElementType& operator[] (const int index) const
  {
    return elements[index];
  }

  // const_reference operator[]( size_type pos ) const;

  /** Assignment operator - creates a deep copy of this array and returns it. */
  rsArray& operator= (const rsArray& other)
  {
    // catch self-assignment:
    if(this == &other)
      return *this;

    // delete old element-array:
    if(elements != NULL)
      delete[] elements;

    // create new array and copy the data:
    elements     = new ElementType[other.numUsed];
    numUsed      = other.numUsed;
    numAllocated = other.numUsed;
    for(int i=0; i<numUsed; i++)
      elements[i] = other.elements[i];

    return *this;
  }

  //---------------------------------------------------------------------------------------------
  // memory management:

  /** Ensures that the array can hold at least 'minNumElements' elements. */
  void ensureAllocatedSize(const int minNumElements)
  {
    if(numAllocated < minNumElements)
    {
      int sizeToAllocate = RAPT::rsNextPowerOfTwo(minNumElements);

      // we must enlarge the allocated memory:
      ElementType* tmpElements = new ElementType[sizeToAllocate];

      // copy the elements form the old array to the new larger array:
      for(int i=0; i<numUsed; i++)
        tmpElements[i] = elements[i];

      // delete the old array which has become too small:
      if(elements != NULL)
      {
        delete[] elements;
        elements = NULL;
      }

      // re-direct the elements pointer to the new allocated memory:
      elements = tmpElements;

      // update the allocated size:
      numAllocated = sizeToAllocate;
    }
  }

  //-----------------------------------------------------------------------------------------------
  // compatibility functions to let it be a replacement for std::vector

  size_t size() const { return (size_t)numUsed; }
  bool empty()  const { return isEmpty(); }
  void push_back(ElementType x) { appendElement(x); }

protected:

  ElementType* elements;
  int numUsed;

  int numAllocated;
    // if we always shrink the number of allocated slots when possible, we can get rid of that member and instead use the invariant:
    // numAllocated = nextPowerOfTwo(numUsed) - that maybe wrapped into a function getNumAllocated
    // but this requires to really enforce this invariant (there are places (copy-constructor, etc.) where we allocate just enough memory
    // to hold the required number of elements

};







} // end namespace rosic

#endif //  rosic_Array_h
