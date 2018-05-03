#ifndef RS_KEYVALUESTRINGTREE_H
#define RS_KEYVALUESTRINGTREE_H

namespace RSLib
{


  /**

  This class augments the rsKeyValueStringArray class by an additional name.

  \todo: rename to rsKeyValueTreeNode
  \todo maybe absorb this class in KeyValueStringArray

  */

  class RSLib_API rsKeyValueStringNodeData : public rsKeyValueStringArray
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. You must pass a name for the node here. */
    rsKeyValueStringNodeData(const rsString& nodeName = rsString::empty) 
    { 
      name = nodeName; 
    }


    /** \name Setup */

    /** Changes the name of this node. */
    void setName(const rsString& newName) { name = newName; }


    /** \name Inquiry */

    /** Returns the name of this node. */
    rsString getName() const { return name; }

    /** Check whether the node has the given name. */
    bool hasName(const rsString& nameToCheck) const { return nameToCheck == name; }


    /** \name Operators */

    /** Compares two KeyValueStringNodeData objects for equality. */
    bool operator==(const rsKeyValueStringNodeData& other) const
    {
      return keyValuePairs == other.keyValuePairs && name == other.name;
    }

    /** Compares two KeyValueStringNodeData objects for inequality. */
    bool operator!=(const rsKeyValueStringNodeData& other) const
    {
      return !(*this == other);
    }

  protected:

    /** \name Data */

    rsString name;

  };

  //===============================================================================================

  /**

  A class for representing a tree structure where each node has a name and an array of key/value 
  pairs - thus, this class can be seen as an abstract form of the structure of XML documents - the 
  name of the node corrresponds to the XML tag-name, the array of key/value pairs corresponds to 
  the XML attributes and the subtrees correspond to the XML child elements.

  \todo: factor out a class XmlDocumentCreator/Parser ...or something...the parser possibly with 
  template methods that are implemented here to set up things from attributes/child-elements

  */

  class RSLib_API rsKeyValueStringTree : public rsTree<rsKeyValueStringNodeData>
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. You must pass a name for the node here. A KeyValueStringTree is always 
    constructed as root-node, but it may become a child-node later by adding it to another node. 
    This implies that you must construct KeyValueStringTree objects on the heap (via new) if they
    are going to become child-nodes of some other nodes because the parent node will then take over 
    ownership of the object and eventually delete it. */
    rsKeyValueStringTree(const rsString& rootNodeName) { nodeData.setName(rootNodeName); }


    /** \name Setup */

    // these accessor methods just delegate the access through to the 'value' member inherited from 
    // rsTree<rsKeyValueStringNodeData> - the 'value' member is thus the actual key/value array

    /** Adds a new string value with the given key/value to the array or updates an existing value 
    with the given key. */
    void setStringValue(const rsString& key, const rsString& value) 
    { 
      nodeData.setStringValue(key, value); 
    }

    /** Adds a new real number value with the given key/value to the array or updates an existing 
    value with the given key. */
    void setRealNumberValue(const rsString& key, const double value) 
    { 
      nodeData.setRealNumberValue(key, value); 
    }

    /** Adds a new integer number value with the given key/value to the array or updates an 
    existing value with the given key. */
    void setIntegerValue(const rsString& key, const int value) 
    { 
      nodeData.setIntegerValue(key, value); 
    }

    /** Adds a new boolean (true/false) value with the given key/value to the array or updates an 
    existing value with the given key. */
    void setBooleanValue(const rsString& key, const bool value) 
    { 
      nodeData.setBooleanValue(key, value); 
    }

    /** Clears this KeyValueStringTree which means deleting all key/value pairs and all 
    subtrees. */
    void clear() 
    { 
      deleteAllKeyValueStringPairs(); deleteAllSubTrees(); 
    }

    /** Deletes all key/value pairs of this node (but leaves alone the subtrees, if any). */
    void deleteAllKeyValueStringPairs() 
    { 
      nodeData.clear(); 
    }


    /** \name Inquiry */

    /** Returns the string value at the given key (if key present, otherwise returns the 
    'defaultReturnValue'). */
    rsString getStringValue(const rsString& key, 
                            const rsString& defaultReturnValue = rsString::empty)
    { 
      return nodeData.getStringValue(key, defaultReturnValue); 
    }

    /** Returns the real number value at the given key (if key present, otherwise returns the 
    'defaultReturnValue'). */
    double getRealNumberValue(const rsString& key, const double defaultReturnValue = 0.0)
    { 
      return nodeData.getRealNumberValue(key, defaultReturnValue); 
    }

    /** Returns the integer number value at the given key (if key present, otherwise returns the 
    'defaultReturnValue'). */
    int getIntegerValue(const rsString& key, const int defaultReturnValue = 0)
    { 
      return nodeData.getIntegerValue(key, defaultReturnValue); 
    }

    /** Returns the boolean (true/false) value at the given key (if key present, otherwise returns 
    the 'defaultReturnValue'). */
    bool getBooleanValue(const rsString& key, const bool defaultReturnValue = false)
    { 
      return nodeData.getBooleanValue(key, defaultReturnValue); 
    }


    /** \name XML Document Creation/Parsing */

    /** Returns a string that represents this tree in form of an XML document including the 
    xml-prolog. */
    rsString getAsXmlDocument();

    /** Sets up this KeyValueStringTree from the passed xml-document. */
    void setFromXmlDocument(const rsString& xmlDocument);

  protected:

    /** \name Internal Functions */

    /** Returns the document element for the XML representation of this tree. This forms the actual 
    body of th xml document, that is everything except the xml-prolog. */
    rsString getXmlDocumentElement(int indentionLevel = 0);

    /** Return the prolog for XML documents. */
    rsString getXmlProlog();

    /** Applies some pre-processing steps to the an xml-document to facilitate subsequent 
    parsing. */
    void preProcessXmlDocumentForParsing(rsString& stringToPreProcess);

    /** Sets this node up from the the given xml document (which is assumed to be pre-processed by 
    preProcessXmlDocumentForParsing). */
    void setupFromXmlString(const rsString& documentToParse);

    /** Removes the xml-prolog form the passed string. */
    void removeXmlProlog(rsString& stringToRemovePrologFrom);

    /** Returns the name of the root xml-element in a pre-processed xml-document. */
    rsString getRootXmlElementName(const rsString& documentToParse);

    /** Scans a pre-processed xml-document for an xml-element with the given element-name and 
    returns the number of attributes which this
    element has. */
    int getNumAttributes(const rsString& documentToParse, const rsString& elementName);

    /** Extracts the attributes of the root xml-element in a pre-processed xml-document and stores 
    them here in our nodeData member. The return value is the position inside the document of the 
    opening angle bracket of the first child element or the opening angle bracket of xml 
    end-tag. */
    int extractAttributes(const rsString& documentToParse);

    /** Given the index of an opening angle bracket '<' of a start-tag inside some xml document, 
    this function returns the index of the
    final character '>' in the corresponding end-tag. */
    int findMatchingEndTag(const rsString& document, const int startTagIndex);

    /** Given the index of an opening angle bracket '<' inside some xml document, this function 
    returns the name of the tag. */
    rsString getTagName(const rsString& document, const int openingBracketIndex);

    /** Creates a child-node and sets it up from the passed xml-string. */
    void createChildNodeFromXmlString(const rsString& xmlString);

  };

}

#endif
