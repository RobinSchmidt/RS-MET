#ifndef RS_ARRAY_H
#define RS_ARRAY_H

namespace RSLib
{

  /**

  This class implements a dynamic array which can grow and shrink as needed. Whenever an element
  is added that does not fit into the currently allocated memory, a new memory area is allocated,
  the size of which is given by the next power of two of the number of required elements. Whenever
  the number of used elements falls below half the size of the allocated memory, the allocated
  memory is shrunken (also to the next power of two of required elements). Access to array-elements
  is not thread-safe If you need to access the array concurrently, you must take care for the
  synchronization yourself. The maximum length is 2^31-1 (maximum value of a 32-bit signed
  integer) - we use a signed integer for indexing in order to use -1 as special return value for
  the functions that return an index, for example, when searching for the index of a certain value,
  to indicate, that no such value was found in the array.

  */

  template<class ElementType>
  class rsArrayTools
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. If sizeToAllocate is nonzero, memory will be reserved to hold the given
    number of elements but the array will still be considered empty. */
    rsArrayTools(const int sizeToAllocate = 0)
    {
      initMembers();
      if( sizeToAllocate > 0 )
      {
        elements     = new ElementType[sizeToAllocate];
        numAllocated = sizeToAllocate;
      }
    }

    /** Constructor. Creates an rsArrayTools from the given c-array of elements. */
    rsArrayTools(const ElementType *initialElements, const int numberOfInitialElements)
    {
      initMembers();
      if( numberOfInitialElements > 0 )
      {
        elements     = new ElementType[numberOfInitialElements];
        numAllocated = numberOfInitialElements;
        numUsed      = numberOfInitialElements;
        rsCopyBuffer(initialElements, elements, numUsed);
      }
    }

    /** Copy-constructor - creates a deep copy of the other array. */
    rsArrayTools(const rsArrayTools<ElementType>& other)
    {
      initMembers();
      copyDataFrom(other);
    }

    /** Conversion constructor. Converts a std::vector of the given type into an rsArrayTools of the
    same type. */
    rsArrayTools(const std::vector<ElementType> stdVector)
    {
      initMembers();
      ensureAllocatedSize((int)stdVector.size());
      for(unsigned int i = 0; i < stdVector.size(); i++)
        appendElement( stdVector[i] );
    }

    /** Destructor. */
    ~rsArrayTools()
    {
      delete[] elements;
    }


    /** \name Element Insertion */

    /** Appends an element to the end of this array. */
    void appendElement(const ElementType newElement)
    {
      ensureAllocatedSize(numUsed+1);
      elements[numUsed] = newElement;
      numUsed++;
    }

    /** Appends an element to this array when the element is not already there (based on '=='
    operator of ElementType). */
    void appendIfNotAlreadyThere(const ElementType newElement)
    {
      if( !hasElement(newElement) )
        appendElement(newElement);
    }

    /** Inserts an element to a given position. All elements after this position will be moved one
    position to the right to make up
     space. */
    void insertElement(const int index, const ElementType newElement)
    {
      if( index >= 0 && index < numUsed )
      {
        ensureAllocatedSize(numUsed+1);
        for(int i = numUsed-1; i >= index; i--)
          elements[i+1] = elements[i];
        elements[index] = newElement;
        numUsed++;
      }
      else if( index >= numUsed )
        appendElement(newElement);
      else
        RS_DEBUG_BREAK;  // negative array index
    }

    /** Appends another array to the end of this one. */
    void appendArray(const rsArrayTools<ElementType>& other)
    {
      ensureAllocatedSize(numUsed+other.numUsed);
      for(int i = 0; i < other.numUsed; i++)
        elements[numUsed+i] = other.elements[i];
      numUsed += other.numUsed;
    }

    // todo: write insertArray


    /** \name Element Removal */

    /** Removes an element from a given position. All elements after this position will be moved
    one position to the left to fill the gap. */
    void removeElementByIndex(const int index)
    {
      rsAssert( index >= 0 && index < numUsed );
      for(int i = index; i < numUsed-1; i++)
        elements[i] = elements[i+1];
      numUsed--;
      shrinkMemoryOccupationIfAppropriate();
    }

    /** Removes an element from the array, if present. */
    void removeElementByValue(const ElementType elementToRemove)
    {
      int index = findElement(elementToRemove);
      if( index != -1 )
        removeElementByIndex(index);
    }

    /** Removes the range between (and including) startIndex and endIndex.  */
    void removeRange(const int startIndex, const int endIndex)
    {
      rsAssert( startIndex >= 0 && startIndex < numUsed );
      rsAssert( endIndex   >= 0 && endIndex   < numUsed );
      rsAssert( endIndex   >= startIndex );
      int numElementsToRemove = endIndex-startIndex+1;
      for(int i = startIndex; i < numUsed-numElementsToRemove; i++)
        elements[i] = elements[i+numElementsToRemove];
      numUsed -= numElementsToRemove;
      shrinkMemoryOccupationIfAppropriate();
    }

    /** Removes all elements from the array that match one of the elements in the parameter array
    (via the '==' operator of ElementType).
    \todo maybe return the number of elements removed == numUsed (old) - numUsed (new)
    */
    void removeMatchingElements(const rsArrayTools& elementsToMatch)
    {
      numUsed = rsCopyIfNotMatching(elements, elements, numUsed, elementsToMatch.elements,
        elementsToMatch.numUsed);
      shrinkMemoryOccupationIfAppropriate();
    }

    /** Removes all elements from the array that do not match one of the elements in the parameter
    array (via the '==' operator of ElementType). */
    void keepOnlyMatchingElements(const rsArrayTools elementsToMatch)
    {
      numUsed = rsCopyIfMatching(elements, elements, numUsed, elementsToMatch.elements,
        elementsToMatch.numUsed);
      shrinkMemoryOccupationIfAppropriate();

    }

    /** Clears the array. */
    void clear()
    {
      numUsed = 0;
      minimizeMemoryOccupation();
    }

    /** Clears the array without de-allocating the reserved memory. */
    void clearWithoutDeAllocation()
    {
      numUsed = 0;
    }

    /** Repeats this array the given number of times. You may also pass zero in which case the
    result will be an empty array. */
    void repeat(const int numberOfTimes)
    {
      rsArrayTools tmp = *this;
      clear();
      ensureAllocatedSize(numberOfTimes*numUsed);
      for(int i = 0; i < numberOfTimes; i++)
        appendArray(tmp);
      // \TODO: move to String
    }

    /** Replaces an element with a new one - if there is no element with the given index, nothing
    will be done. */
    void replaceElement(const int index, const ElementType newElement)
    {
      rsAssert( index >= 0 && index < numUsed ); // invalid array index
      elements[index] = newElement;
    }

    /** Returns a copy of the element at the given index. */
    ElementType at(const int index) const
    {
      rsAssert( index >= 0 && index < numUsed ); // invalid array index
      return elements[index];
    }


    /** \name Manipulations */

    void reverse()
    {
      rsReverse(elements, numUsed);
    }


    /** \name Inquiry */

    /** Searches the array for the element (using the == operator) and returns the index or -1, if
    the element is not found. */
    int findElement(const ElementType elementToFind)
    {
      for(int i = 0; i < numUsed; i++)
      {
        if( elements[i] == elementToFind )
          return i;
      }
      return -1;
      // \todo: turn this into a template function, rename into findFirstOccurrenceOf, extend
      // signature to inclue a search start (optional, default zero)
    }

    /** Searches for the last occurrence of an element in this array and returns its index or -1,
    if the element is not found. The optional searchStart parameter can be used to specify, where
    the search should start. If a value < 0 is passed (for example, -1), the search will start at
    the very end of the array, which is the default. */
    int rsFindLastOccurrenceOf(const ElementType &elementToFind, const int searchStart = -1) const
    {
      return RSLib::rsFindLastOccurrenceOf(elements, numUsed, elementToFind, searchStart);
    }

    /** Returns true, if the array contains the element, false otherwise. */
    bool hasElement(const ElementType elementToFind)
    {
      return findElement(elementToFind) != -1;
    }

    /** Returns the number of elements in the array. */
    int getNumElements() const
    {
      return numUsed;
    }

    /** Returns the number of allocated (but not necessarily used) elements in the array - client
    code should actually not really be concerned about this, this function merely facilitates
    testing. */
    int getNumAllocatedElements() const
    {
      return numAllocated;
    }

    /** Returns true when the array is empty, false otherwise. */
    bool isEmpty() const
    {
      return numUsed == 0;
    }

    /** Searches for the first occurrence of a sub-array and returns the index where it starts,
    -1 if the searched sub-array is not found in this array. The optional searchStart parameter can
    be used to skip an initial section in the search.
    \todo: switch to a more efficient algorithm (Knuth/Morris/Pratt) */
    int findFirstOccurrenceOf(const rsArrayTools<ElementType>& subArrayToFind,
                              const int searchStart = 0) const
    {
      return RSLib::rsFindFirstOccurrenceOf(elements, numUsed, subArrayToFind.elements,
                                             subArrayToFind.numUsed, searchStart);
    }
    // \TODO: move to String

    /** Returns a copy of the element at the given index. */
    ElementType getElement(const int index) const { return elements[index]; }

    /** Returns the sub-array between (and including) startIndex and endIndex.  */
    rsArrayTools<ElementType> getSubArray(const int startIndex, const int endIndex) const
    {
      rsAssert( startIndex >= 0 && startIndex < numUsed );
      rsAssert( endIndex   >= 0 && endIndex   < numUsed );
      rsAssert( endIndex   >= startIndex );
      rsArrayTools<ElementType> result;
      int resultLength = endIndex-startIndex+1;
      result.ensureAllocatedSize(resultLength);
      for(int i = 0; i < resultLength; i++)
        result.appendElement(elements[startIndex+i]);
      return result;
    }

    /** Returns a pointer to the raw data. Use this with care, since it breaks the
    data-encapsulation of the class */
    ElementType* getRawData()
    {
      return elements;
    }

    /** Copies the array elements into the passed buffer. The caller is responsible for ensuring,
    that the buffer is large enough. */
    void copyDataTo(ElementType *buffer)
    {
      rsCopyBuffer(elements, buffer, numUsed);
    }

    /** \name Operators */

    /** Allows fast indexed acces to array elements - warning no bounds checking is done. */
    inline ElementType& operator[] (const int index)
    {
      rsAssert( index >= 0 && index < numUsed ); // invalid array index
      return elements[index];
    }

    /** Assignment operator - creates a deep copy of this array and returns it. */
    rsArrayTools<ElementType>& operator= (const rsArrayTools<ElementType>& other)
    {
      if( this != &other )   // nothing to do in case of self-assignment
        copyDataFrom(other);
      return *this;
    }

    /** Compares two arrays for equality. Two arrays are considered equal, if they have the same
    length and match element by element. The reserved memory size is irrelevant for this
    comparison.  */
    bool operator==(const rsArrayTools& other) const
    {
      if( numUsed != other.numUsed )
        return false;
      else
        return rsAreBuffersEqual(elements, other.elements, numUsed);
    }

    /** Compares two arrays for inequality. */
    bool operator!=(const rsArrayTools& other) const
    {
      return !(*this == other);
    }


    /** \name Memory Management */

    /** Ensures that the array can hold at least 'minNumElements' elements. */
    void ensureAllocatedSize(const int minNumElements)
    {
      if( numAllocated < minNumElements )
      {
        int newSizeToAllocate = rsNextPowerOfTwo(minNumElements);
        reAllocateMemoryAndMoveData(newSizeToAllocate);
      }
    }

    /** Reduces the amount of occupied memory to the number of actually used elements. */
    void minimizeMemoryOccupation()
    {
      int newSizeToAllocate = numUsed;
      reAllocateMemoryAndMoveData(newSizeToAllocate);
    }


  protected:

    /** Allocates a new memory area for the data, copies the data to the new location and releases
    the old memory area. */
    void reAllocateMemoryAndMoveData(int numElementsToAllocate)
    {
      ElementType* tmpElements = new ElementType[numElementsToAllocate];
      for(int i = 0; i < numUsed; i++)
        tmpElements[i] = elements[i];
      delete[] elements;
      elements     = tmpElements;
      numAllocated = numElementsToAllocate;
    }

    /** Shrinks the allocated memory size to the next power of two of the actual number of used
    elements. */
    void shrinkMemoryOccupationIfAppropriate()
    {
      if(numUsed <= numAllocated / 2)
      {
        int newSizeToAllocate = rsNextPowerOfTwo(numUsed);
        reAllocateMemoryAndMoveData(newSizeToAllocate);
      }
    }

    /** Copies the data from another array into this one. */
    void copyDataFrom(const rsArrayTools& sourceArray)
    {
      ensureAllocatedSize(sourceArray.numUsed);
      numUsed = sourceArray.numUsed;
      for(int i = 0; i < numUsed; i++)
        elements[i] = sourceArray.elements[i];
      shrinkMemoryOccupationIfAppropriate();
    }

    /** Initializes the data members. */
    void initMembers()
    {
      elements     = NULL;
      numUsed      = 0;
      numAllocated = 0;
    }


    /** \name Data */

    ElementType* elements;  // the C-array containing the elements
    int numUsed;            // number of used elements
    int numAllocated;       // number of allocated elements

  };

}

#endif
