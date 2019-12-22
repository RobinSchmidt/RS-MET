#ifndef RS_TREE_H
#define RS_TREE_H

namespace RSLib
{

  /**

  This class implements a tree structure of data-elements of arbitrary type. Each node can have any 
  number of child-nodes (aka sub-nodes or sub-trees) and has either one or zero parent nodes (where 
  in the latter case, the node represents the root of the tree). Nodes that do not have any 
  sub-nodes are considered as leaf nodes.

  \todo: maybe split the responsibilities between a TreeNode and Tree class in a similar way as 
  with ListItem and List

  */

  template<class ElementType>
  class rsTree
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. Creates a node that doesn't have any sub-nodes nor a parent node - thus the 
    'trivial' tree is root and leaf at the same time. The value that should be associated with this 
    node can be passed in here. If nothing is passed, a default object is created via the default 
    constructor of class 'ElementType' - so the class must provide a default constructor. */
    rsTree(ElementType dataAtThisNode = ElementType())
    {
      parentNode = NULL;
      nodeData   = dataAtThisNode;
    }

    /** Destructor. */
    ~rsTree()
    {
      deleteAllSubTrees();
    }


    /** \name Node Insertion/Deletion */

    /** Adds a subtree as child node to this Tree - 'this' object takes over responsibility for 
    deleting it eventually. */
    void hangInSubTree(rsTree<ElementType>* subTree)
    {
      childNodes.appendIfNotAlreadyThere(subTree);
      subTree->parentNode = this;
    }

    void deleteAllSubTrees()
    {
      for(int i = 0; i < childNodes.getNumElements(); i++)
        delete childNodes[i];
      childNodes.clear();
    }

    // \todo void deleteChildNode


    /** \name Node Access */

    /** Returns a pointer to one of the direct child nodes (subtrees) of this node. Note that this 
    node still owns the returned child-node, so don't delete the pointer. If this node does not 
    have a child-node with given index, it will return a NULL pointer. */
    rsTree<ElementType>* getSubTree(int index)
    {
      if( index >= 0 && index < childNodes.getNumElements() )
        return childNodes[index];
      else
        return NULL;
    }

    /** Returns a pointer to the parent node of this node. If this node does not have a parent-node 
    (i.e. it is the root node), it will return a NULL pointer. */
    rsTree<ElementType>* getParentNode() const 
    { 
      return parentNode; 
    }

    /** Returns the data at this node by value, so this is recommended only for small data-objects 
    at the nodes. */
    ElementType getNodeDataByValue() const 
    { 
      return nodeData; 
    }

    /** Returns a reference (i.e. pointer) to the data at this node. You may use the pointer to 
    manipulate the data but don't delete it as it is owned here. */
    ElementType* getNodeDataByReference() 
    { 
      return &nodeData; 
    }


    /** \name Inquiry */

    /** Returns true when this tree/node is the root node, false otherwise. */
    bool isRoot() const 
    { 
      return parentNode == NULL; 
    }

    /** Returns true when this tree/node is a leaf node, false otherwise. */
    bool isLeaf() const 
    { 
      return childNodes.isEmpty(); 
    }

    /** Returns the total number of nodes in this tree (including the node on which this function 
    is called itself). */
    int getTotalNumberOfNodes()
    {
      int result = 1;
      for(int i = 0; i < childNodes.getNumElements(); i++)
        result += childNodes[i]->getTotalNumberOfNodes();
      return result;
    }

    /** Returns the number of child nodes which are directly below this node, so this function is 
    non-recursive in the sense that it doesn't count the children's children ...and so forth. */
    int getNumberOfDirectChildNodes() { return childNodes.getNumElements(); }

        
    /** \name Operators */

    /** Compares two trees for equality. Two trees are considered equal, if their nodeData fields 
    at the root notes are equal to each other, the trees have the same number of subtrees and all 
    the subtrees satisfy the equality condition recursively. */
    bool operator==(const rsTree& other) const
    {
      if( nodeData != other.nodeData )
        return false;
      if( childNodes.getNumElements() != other.childNodes.getNumElements() )
        return false;
      for(int i = 0; i < childNodes.getNumElements(); i++)
      {
        if( !(*childNodes.getElement(i) == *other.childNodes.getElement(i)) )
          return false;
      }
      return true;
    }

    /** Compares two trees for inequality by negation of the equality operator. */
    bool operator!=(const rsTree& other) const
    {
      return !(*this == other);
    }

  protected:

    /** \name Data */

    ElementType                     nodeData;    // the data/value at this node
    rsTree<ElementType>             *parentNode; // pointer to the parent node (NULL in case of root node)
    rsArrayTools< rsTree<ElementType>* > childNodes;  // pointers to the child-nodes (empty in case of leaf node)

    // \todo: maybe factor out a subclass PointerTree where the data is just a pointer to void 
    // (maybe with an additional byteSize field) - the PointerTree class does then not need to be 
    // a template class and can be implemented in the .cpp file. a derived template class could 
    // then just use typecasts in the node-data getters

  };

}

#endif
