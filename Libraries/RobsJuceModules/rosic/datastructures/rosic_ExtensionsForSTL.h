#ifndef rosic_ExtensionsForSTL_h
#define rosic_ExtensionsForSTL_h

//#include <vector>
//#include "../basics/rosic_FunctionTemplates.h"

// todo: merge with RAPT StandardContainerFunctions

namespace rosic
{

  /**

  This file contains functions that extend the functionality of container classes of the C++ Standard Template Library (STL). Whenever you 
  think that an STL container ought to have some member-function but doesn't, look here. For example, you might want to inquire whether or
  not a std::vector contains a certain element and realize in disappointment that a function std::vector::containsElement() (or similar), 
  doesn't exist. No problem - simply use containsElement(const std::vector<ElementType> &vectorToLookIn, const ElementType &elementToFind) 
  from this function-library. Likewise, when you want to implement new "pseudo-member" functions for STL-containers yourself, here is the 
  right place to put them. Some of the functions here may simply restate std's member functions with a different syntax and act as 
  replacements (like, for example, appendElement can replace push_back) - we do this in the name of code-cleanliness, because it is ugly to
  mix calls to the std-member functions with calls to the pseudo-member functions defined here. Moreover, the functions here are hoped to 
  be more self-explanatory and easier to use than the std-functions.

  Why this kludge?
  In an ideal world, we would prefer to define our own containers and directly define all member functions we need in the way we like. But 
  on the other hand, taking advantage of the STL makes some things easier, so we have to live with this ugly antipattern. Another 
  possibility would be to define a template-class that inherits from std::vector or yet another would be to define a template class that 
  has a std::vector as member and then delegate all the requests. This solution here was chosen because it makes it easier to debug 
  programs in MS Visual Studio - the debugger there can directly show the contents of STL containers, so - for ease of debugging - it's 
  beneficial to directly use STL containers inside the code whenever possible. Burying the underlying STL container in a class (either by 
  inheritance or by composition) would add another level to the debugger's tree-view for the variables and thus, make it clunkier to 
  inspect the container contents.

  */

    
    
  /** Appends an element at the end of the vector. */
  template<class ElementType>
  void appendElement(std::vector<ElementType> &vectorToAppendElementTo, const ElementType &newElement);
    
  /** Appends an element at the end of the vector but only when the element in question does not already exist (via the == operator). */
  template<class ElementType>
  void appendIfNotAlreadyThere(std::vector<ElementType> &vectorToAppendElementTo, const ElementType &newElement);

  /** Appends the second vector to the first one. */
  template<class ElementType>
  void appendVector(std::vector<ElementType> &vectorToAppendTo, const std::vector<ElementType> &appendixVector);
    
  /** Returns true, if the vector contains the element, false otherwise. */  
  template<class ElementType>
  bool containsElement(const std::vector<ElementType> &vectorToLookIn, const ElementType &elementToFind);

  /** Searches the vector for the element (using the == operator) and returns the index or -1 if the element is not found. The function 
  assumes no particular order of the elements - it searches linearly through the whole vector and thus has complexity O(N). */
  template<class ElementType>
  int findElement(const std::vector<ElementType> &vectorToLookIn, const ElementType &elementToFind);

  /** Returns the number of elements in the passed vector. */
  template<class ElementType>
  int numberOfElements(const std::vector<ElementType> &vector);

  /** Removes an element from a vector at a given position. All elements after this position will be moved one position up to fill the 
  gap. */  
  template<class ElementType>
  void removeElementByIndex(std::vector<ElementType> &vectorToRemoveElementFrom, const int index);

  /** Removes an element from the vector, if present. If the element is present multiple times, it will remove only the first occurence. */
  template<class ElementType>  
  void removeElementByValue(std::vector<ElementType> &vectorToRemoveElementFrom, const ElementType &elementToRemove);

  /** Removes all occurences of the given element from the vector. */
  template<class ElementType>  
  void removeAllOccurencesOfElement(std::vector<ElementType> &vectorToRemoveElementsFrom, const ElementType &elementToRemove);

  /** Removes all occurences of elements in the first vector that are also present in the second vector. Or - to put it the other way 
  around - it retains only those elements in the first vector that are not at the same time elements of the second vector. */
  template<class ElementType>  
  void removeAllOccurencesOfElements(std::vector<ElementType> &vectorToRemoveElementsFrom, 
                                     const std::vector<ElementType> &elementsToRemove);

  //=======================================================================================================================================
  // implementation:


  template<class ElementType>
  void appendElement(std::vector<ElementType> &vectorToAppendElementTo, const ElementType &newElement)
  {
    vectorToAppendElementTo.push_back(newElement);
  }

  template<class ElementType>
  void appendIfNotAlreadyThere(std::vector<ElementType> &vectorToAppendElementTo, const ElementType &newElement)
  {
    if( !containsElement(vectorToAppendElementTo, newElement) )
      appendElement(vectorToAppendElementTo, newElement);
  }

  template<class ElementType>
  void appendVector(std::vector<ElementType> &vectorToAppendTo, const std::vector<ElementType> &appendixVector)
  {
    vectorToAppendTo.reserve(vectorToAppendTo.size() + appendixVector.size());
    vectorToAppendTo.insert(vectorToAppendTo.end(), appendixVector.begin(), appendixVector.end());
  }

  template<class ElementType>
  bool containsElement(const std::vector<ElementType> &vectorToLookIn, const ElementType &elementToFind)
  {
    return findElement(vectorToLookIn, elementToFind) != -1;
  }

  template<class ElementType>
  int findElement(const std::vector<ElementType> &vectorToLookIn, const ElementType &elementToFind)
  {
    for(int i = 0; i < (int) vectorToLookIn.size(); i++)
    {
      if( vectorToLookIn[i] == elementToFind )
        return i;
    }
    return -1;
  }

  template<class ElementType>
  int numberOfElements(const std::vector<ElementType> &vector)
  {
    return vector.size();
  }

  template<class ElementType>
  void removeElementByIndex(std::vector<ElementType> &vectorToRemoveElementFrom, const int index)
  {
    vectorToRemoveElementFrom.erase(vectorToRemoveElementFrom.begin() + index);
  }

  template<class ElementType>  
  void removeElementByValue(std::vector<ElementType> &vectorToRemoveElementFrom, const ElementType &elementToRemove)
  {
    int index = findElement(vectorToRemoveElementFrom, elementToRemove);
    if( index != -1 )
      removeElementByIndex(vectorToRemoveElementFrom, index);
  }

  template<class ElementType>  
  void removeAllOccurencesOfElement(std::vector<ElementType> &vectorToRemoveElementsFrom, const ElementType &elementToRemove)
  {
    int i = (int) vectorToRemoveElementsFrom.size() - 1;
    while( i >= 0 )
    {
      if( vectorToRemoveElementsFrom[i] == elementToRemove )
        removeElementByIndex(vectorToRemoveElementsFrom, i);
      i--;
    }
  }

  template<class ElementType>  
  void removeAllOccurencesOfElements(std::vector<ElementType> &vectorToRemoveElementsFrom, 
                                     const std::vector<ElementType> &elementsToRemove)
  {
    for(unsigned int i = 0; i < elementsToRemove.size(); i++ )
      removeAllOccurencesOfElement(vectorToRemoveElementsFrom, elementsToRemove.at(i));
  }


} // end namespace rosic

#endif 
