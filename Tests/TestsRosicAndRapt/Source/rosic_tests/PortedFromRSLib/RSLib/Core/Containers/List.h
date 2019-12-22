#ifndef RS_LIST_H
#define RS_LIST_H

namespace RSLib
{

  /**

  This class implements an item/element of a doubly linked list of data-elements of arbitrary type. 
  The item contains the data itself (the type of which is determined by the template parameter) and 
  a pointer to the next and previous item. When the 'next' pointer is NULL, the item is the last 
  one in the list, when the previous pointer is NULL, the item is the first in the list.

  */

  template<class ElementType>
  class rsListItem
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. Creates a list item without predecessor an successor (next and previous 
    pointers are NULL). The value that should be associated with this item can be passed in here. 
    If nothing is passed, a default object is created via the default constructor of class 
    'ElementType' - so the class must provide a default constructor. */
    rsListItem(ElementType dataAtThisItem = ElementType())
    {
      predecessor = NULL;
      successor   = NULL;
      itemData    = dataAtThisItem;
    }


    /** \name Element Insertion */

    /** Inserts the passed new item after this one such that it becomes the new successor of this 
    one. Any existing successor will then become the successor of the item that is being 
    inserted. */
    void insertItemAfter(rsListItem<ElementType>* newItem)
    {
      newItem->predecessor = this;
      newItem->successor   = this->successor;
      if( successor != NULL )
        successor->predecessor = newItem;
      this->successor = newItem;
    }

    /** Inserts the passed new item before this one such that it becomes the new predecessor of 
    this one. Any existing predecessor will then become the predecessor of the item that is being 
    inserted. */
    void insertItemBefore(rsListItem<ElementType>* newItem)
    {
      newItem->successor   = this;
      newItem->predecessor = this->predecessor;
      if( predecessor != NULL )
        predecessor->successor = newItem;
      this->predecessor = newItem;
    }


    /** \name Element Access */

    /** Returns a pointer to the predecessor of this item. If this item does not have a predecessor 
    (i.e. it is the first item), it will return a NULL pointer. */
    rsListItem<ElementType>* getPredecessor() const 
    { 
      return predecessor; 
    }

    /** Returns a pointer to the successor of this item. If this item does not have a successor 
    (i.e. it is the last item), it will return a NULL pointer. */
    rsListItem<ElementType>* getSuccessor() const 
    { 
      return successor; 
    }

    /** Returns the data at this item by value, so this is recommended only for small data-objects 
    at the items. */
    ElementType getItemDataByValue() const 
    { 
      return itemData; 
    }

    /** Returns a reference (i.e. pointer) to the data at this item. You may use the pointer to 
    manipulate tha data but don't delete it as it is owned here. */
    ElementType* getItemDataByReference() 
    { 
      return &itemData; 
    }


    /** \name Inquiry */

    /** Returns true when this item is the first item in the list, false otherwise. */
    bool isFirst() const 
    { 
      return predecessor == NULL; 
    }

    /** Returns true when this item is the last item in the list, false otherwise. */
    bool isLast() const 
    { 
      return successor == NULL; 
    }

  //protected:
     
    /** \name Data */

    ElementType             itemData;      // data/value at this node
    rsListItem<ElementType> *predecessor;  // pointer to previous item (NULL in case of first item)
    rsListItem<ElementType> *successor;    // pointer to next item (NULL in case of last item)

  };


  //===============================================================================================

  /**

  This class implements a doubly linked list of data-elements of arbitrary type.

  \todo: this class is still incomplete and untested
  \todo: write subclasses rsStack and rsQueue inheriting privately and defining 
         pushElement/popElement, enqueueElement/dequeueElement
  \todo: maybe have classes rsPointerList, rsPointerArray, rsPointerStack, etc that optionally 
         delete the pointed-to objects on destruction - avoids the template-associated (binary)
         code-bloat because the container-logic needs to exist only once in the dll. Moreover, 
         there's will be much less implicit copying (via copy constructors) going on when setting 
         and retrieving elements

  */

  template<class ElementType>
  class rsList
  {

  public:


    /** \name Construction/Destruction */

    /** Constructor. Creates a list item without predecessor an successor (next and previous 
    pointers are NULL). The value that should be associated with this item can be passed in here. 
    I nothing is passed, a default object is created via the default constructor of class 
    'ElementType' - so the class must provide a default constructor. */
    rsList()
    {
      firstItem = NULL;
      lastItem  = NULL;
    }

    /** Destructor. */
    ~rsList()
    {
      deleteAllItems();
    }


    /** \name Element Insertion */

    /** Appends the passed item at the end of the list. */
    void appendItem(ElementType newItem) { appendItem(new rsListItem<ElementType>(newItem)); }

    /** Appends the passed item at the end of the list. The List object takes over responsibility for 
    eventual destruction.
    \todo: maybe move into protected area and leave only prependItem(ElementType newItem) public.
    */
    void appendItem(rsListItem<ElementType>* newItem)
    {
      if( isEmpty() )
      {
        firstItem = lastItem = newItem;
        newItem->predecessor = NULL;
        newItem->successor   = NULL;
      }
      else
      {
        lastItem->successor  = newItem;
        newItem->predecessor = lastItem;
        newItem->successor   = NULL;
        lastItem             = newItem;
      }
    }

    /** Prepends the passed item at the beginning of the list. */
    void prependItem(ElementType newItem) { prependItem(new rsListItem<ElementType>(newItem)); }

    /** Prepends the passed item at the beginning of the list. The List object takes over 
    responsibility for eventual destruction. 
    \todo: maybe move into protected area and leave only prependItem(ElementType newItem) public.
    */
    void prependItem(rsListItem<ElementType>* newItem)
    {
      if( isEmpty() )
      {
        firstItem = lastItem = newItem;
        newItem->predecessor = NULL;
        newItem->successor   = NULL;
        // code duplication - see appendItem ...maybe write function insertInitialItem or something
      }
      else
      {
        firstItem->predecessor = newItem;
        newItem->successor     = firstItem;
        newItem->predecessor   = NULL;
        firstItem              = newItem;
      }
    }

    /** Prepends the passed item at the beginning of the list. */
    void insertSorted(ElementType newItem, bool insertAfterEqualItems = true)
    { insertSorted(new rsListItem<ElementType>(newItem), insertAfterEqualItems); }

    /** Inserts the new item into the list such that a presumed existing ascending order is 
    preserved. If the second argument is true, the new element will be inserted after equal items 
    (if any) and before them otherwise. */
    void insertSorted(rsListItem<ElementType>* newItem, bool insertAfterEqualItems = true)
    {
      if( insertAfterEqualItems == true )
      {
        rsListItem<ElementType>* newPredecessor = getLastSortedPredecessorFor(newItem);
        if( newPredecessor == NULL )
          prependItem(newItem);
        else
        {
          newPredecessor->insertItemAfter(newItem);
          if( newItem->successor == NULL )
            lastItem = newItem;
        }
      }
      else
      {
        rsListItem<ElementType>* newSuccessor = getFirstSortedSuccessorFor(newItem);
        if( newSuccessor == NULL )
          appendItem(newItem);
        else
        {
          newSuccessor->insertItemBefore(newItem);
          if( newItem->predecessor == NULL )
            firstItem = newItem;
        }
      }
    }


    /** \name Element Deletion */

    // \todo: appendList, prependList, mergeLists(bool sortedMerge), removeFirstOccurenceOf, 
    // removeLastOccurenceOf, removeAllOccurencesOf, removeFirst, removeLast, getFirst, getLast

    /** Deletes all the items in the list. */
    void deleteAllItems()
    {
      if( isEmpty() )
        return;
      else
      {
        rsListItem<ElementType>* currentItem = firstItem;
        rsListItem<ElementType>* nextItem    = NULL;
        while( currentItem->getSuccessor() != NULL )
        {
          nextItem = currentItem->getSuccessor();
          delete currentItem;
          currentItem = nextItem;
        }
        delete currentItem; // last item has successor == NULL but must be deleted anyway
      }
    }


    /** \name Inquiry */

    /** Returns a pointer to the first item in this list. If the list is empty, it will return a 
    NULL pointer. */
    rsListItem<ElementType>* getFirstItem() const 
    { 
      return firstItem; 
    }

    /** Returns a pointer to the last item in this list. If the list is empty, it will return a 
    NULL pointer. */
    rsListItem<ElementType>* getLastItem() const 
    { 
      return lastItem; 
    }

    /** Returns true when this list is empty, false otherwise. */
    bool isEmpty() const
    {
      if( firstItem == NULL && lastItem == NULL )
        return true;
      else
        return false;
      // return firstItem == lastItem == NULL; ...does not work - why?
    }

    /** Returns the total number of items in this list.
    \todo: maybe store the number of items as member variable.  */
    int getNumberOfItems()
    {
      if( isEmpty() )
        return 0;
      else
      {
        int result = 1;
        rsListItem<ElementType>* currentItem = firstItem;
        while( currentItem->getSuccessor() != NULL )
        {
          result++;
          currentItem = currentItem->getSuccessor();
        }
        return result;
      }
    }

    /** Returns a pointer to the last item in this list after which the new item can be inserted 
    without destroying an ascending order (as defined by the '>' operator of class ElementType). 
    If the list is empty or even the first item in the list is greater then the new item, it will 
    return a NULL pointer. For this to make sense, the list should be already sorted. */
    rsListItem<ElementType>* getLastSortedPredecessorFor(rsListItem<ElementType>* newItem)
    {
      if( isEmpty() || firstItem->itemData > newItem->itemData )
        return NULL;
      else
      {
        rsListItem<ElementType>* currentItem = lastItem;
        while( currentItem->getPredecessor() != NULL && currentItem->itemData > newItem->itemData )
          currentItem = currentItem->getPredecessor();
        return currentItem;
      }
    }

    /** Returns a pointer to the first item in this list before which the new item can be inserted 
    without destroying an ascending order (as defined by the '<' operator of class ElementType). If 
    the list is empty or even the last item in the list is less then the new item, it will return a 
    NULL pointer. For this to make sense, the list should be already sorted. */
    rsListItem<ElementType>* getFirstSortedSuccessorFor(rsListItem<ElementType>* newItem)
    {
      if( isEmpty() || lastItem->itemData < newItem->itemData )
        return NULL;
      else
      {
        rsListItem<ElementType>* currentItem = firstItem;
        while( currentItem->getSuccessor() != NULL && currentItem->itemData < newItem->itemData )
          currentItem = currentItem->getSuccessor();
        return currentItem;
      }
    }

    /** Converts this list into an Array and returns it. */
    rsArrayTools<ElementType> getAsArray()
    {
      if( isEmpty() )
        return rsArrayTools<ElementType>(); // default constructor of Array constructs an empty array
      else
      {
        rsArrayTools<ElementType> listAsArray;
        listAsArray.ensureAllocatedSize(getNumberOfItems());
        rsListItem<ElementType>* currentItem = firstItem;
        listAsArray.appendElement(currentItem->getItemDataByValue());
        while( currentItem->getSuccessor() != NULL )
        {
          currentItem = currentItem->getSuccessor();
          listAsArray.appendElement(currentItem->getItemDataByValue());
        }
        return listAsArray;
      }
    }

        
    /** \name Operators */

    /** Compares two lists for equality. Two lists are considered equal, if they have the same 
    number of items and the data at each of the items matches. */
    bool operator==(const rsList& other) const
    {
      if( getNumberOfItems() != other.getNumberOfItems() )
        return false;
      if( isEmpty() )
        return true; // at this point, we know that both lists are empty

      // compare the lists item-by-item (we already know that they are non-empty and the lengths match):
      rsListItem<ElementType>* currentItemInThisList  = firstItem;
      rsListItem<ElementType>* currentItemInOtherList = other.firstItem;
      while( currentItemInThisList->getSuccessor() != NULL )
      {
        if( *currentItemInThisList != *currentItemInOtherList )
          return false;
        currentItemInThisList  = currentItemInThisList->getSuccessor();
        currentItemInOtherList = currentItemInOtherList->getSuccessor();
      }

      return true;
    }

    /** Compares two lists for inequality by negation of the equality operator. */
    bool operator!=(const rsList& other) const
    {
      return !(*this == other);
    }

  protected:
    
    /** \name Data */

    rsListItem<ElementType> *firstItem;   // pointer to first item (NULL in case of an empty list)
    rsListItem<ElementType> *lastItem;    // pointer to last item (NULL in case of an empty list)

  };

}

#endif
