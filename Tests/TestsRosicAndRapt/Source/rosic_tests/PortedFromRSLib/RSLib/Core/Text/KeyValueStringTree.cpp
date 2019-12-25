using namespace RSLib;

rsString rsKeyValueStringTree::getAsXmlDocument()
{
  rsString xmlString = getXmlProlog();
  xmlString += getXmlDocumentElement();
  xmlString += rsString("\n");
  return xmlString;
}

rsString rsKeyValueStringTree::getXmlDocumentElement(int indentionLevel)
{
  rsString theString;
  rsString indentionSpace = rsString(" ");
  indentionSpace.repeat(indentionLevel);
  int i;

  // start-tag with attributes:
  theString += indentionSpace + rsString("<") + nodeData.getName();
  for(i=0; i<nodeData.getNumKeyValueStringPairs(); i++)
  {
    rsKeyValueStringPair pair = nodeData.getKeyValueStringPair(i);
    theString += rsString("\n") + indentionSpace + rsString(" ") + pair.getKey()
               + rsString("=\"") + pair.getStringValue() + rsString("\"");
  }
  theString += rsString(">");

  // child elements:
  for(i=0; i<getNumberOfDirectChildNodes(); i++)
  {
    rsTree<rsKeyValueStringNodeData> *subTree = getSubTree(i);
    rsKeyValueStringTree *castedSubTree = static_cast<rsKeyValueStringTree*> (subTree);
    theString += rsString("\n") + castedSubTree->getXmlDocumentElement(indentionLevel+1);
  }

  // end tag:
  theString += rsString("\n") + indentionSpace + rsString("</") + nodeData.getName() 
             + rsString(">");

  return theString;
}

void rsKeyValueStringTree::setFromXmlDocument(const rsString& xmlDocument)
{
  clear();
  rsString tmpString = xmlDocument;
  preProcessXmlDocumentForParsing(tmpString);
  nodeData.setName(getRootXmlElementName(tmpString));
  setupFromXmlString(tmpString);
}

rsString rsKeyValueStringTree::getXmlProlog()
{
  return rsString("<?xml version=\"1.0\"?>\n");
}

void rsKeyValueStringTree::preProcessXmlDocumentForParsing(rsString& stringToPreProcess)
{
  removeXmlProlog(stringToPreProcess);
  stringToPreProcess.replaceCharacter('\n', ' ');    // replaces newline-characters with spaces
  stringToPreProcess.removeRepeatedCharacters(' ');  // remove repeated spaces
  stringToPreProcess.removeLeadingCharacters(' ');
  stringToPreProcess.removeTrailingCharacters(' ');
}

void rsKeyValueStringTree::setupFromXmlString(const rsString& documentToParse)
{
  int position = extractAttributes(documentToParse);
  while( documentToParse.at(position+1) != '/' )
  {
    // the bracket at 'position' opens a new child element - extract substring for this child 
    // element, create child node and let the child node setup itself from the substring:
    int      subNodeEnd     = findMatchingEndTag(documentToParse, position);
    rsString childSubString = documentToParse.getSubString(position, subNodeEnd);
    createChildNodeFromXmlString(childSubString);
    position = subNodeEnd+2;
  }
}

void rsKeyValueStringTree::removeXmlProlog(rsString& stringToRemovePrologFrom)
{
  rsRange<int> range = stringToRemovePrologFrom.findRangeEnclosedBy(rsString("<"), 
                         rsString(">"), true);
  stringToRemovePrologFrom.removeRange(range.getMin(), range.getMax());
}

rsString rsKeyValueStringTree::getRootXmlElementName(const rsString& documentToParse)
{
  // \todo: handle invalid documents (no space, no closing bracket....)...or better: validate 
  // document beforehand

  int start = documentToParse.findFirstOccurrenceOf(rsString("<"))+1;

  // scan for whitespace ' ' or closing angle bracket '>':
  int end = start+1;
  while( documentToParse.at(end) != '>' && documentToParse.at(end) != ' ' )
    end++;
  end--; // one too many

  return documentToParse.getSubString(start, end);
}

int rsKeyValueStringTree::getNumAttributes(const rsString& documentToParse, 
                                           const rsString& elementName)
{
  // the number of attributes is given by the number of whitespaces before the next '<' 
  // character minus 1:
  int count = 0;
  int start = documentToParse.findFirstOccurrenceOf(elementName) + elementName.getLength();
  int i     = start;
  while( documentToParse.at(i) != '<' )
  {
    if( documentToParse.at(i) == ' ' )
      count++;
    i++;
  }
  return count-1;
}

int rsKeyValueStringTree::extractAttributes(const rsString& documentToParse)
{
  rsString elementName = getRootXmlElementName(documentToParse);
  int numToExtract     = getNumAttributes(documentToParse, elementName);
  int start            = documentToParse.findFirstOccurrenceOf(elementName) 
                         + elementName.getLength();
  //int i              = start;
  rsString attributeName;
  rsString attributeValue;
  for(int i=1; i<=numToExtract; i++)
  {
    attributeName = documentToParse.getSubStringEnclosedBy(rsString(" "), rsString("="), 
                                                           false, start);
    start += attributeName.getLength()+2;
    attributeValue = documentToParse.getSubStringEnclosedBy(rsString("\""), rsString("\""), 
                                                            false, start);
    start += attributeValue.getLength()+2;
    setStringValue(attributeName, attributeValue);
  }

  // some debug checks:
  //char bracket = documentToParse.at(start+2);
  rsAssert( documentToParse.at(start+2) == '<' ); // should be an opening angle bracket

  return start+2;
}

int rsKeyValueStringTree::findMatchingEndTag(const rsString& document, const int startTagIndex)
{
  rsString tagName  = getTagName(document, startTagIndex);
  rsString startTag = rsString("<") + tagName;
  rsString endTag   = rsString("</") + tagName + rsString(">");
  int searchStart = startTagIndex + startTag.getLength();
  int endTagIndex = document.findFirstOccurrenceOf(endTag, searchStart);

  // endTagIndex so found must be interpreted as a first candidate - if there are alike opening 
  // start-tags before this end-tag we must continue our search for correspondingly later end-tags:
  rsArrayTools<int> otherOccurrencesOfStartTag = document.findAllOccurrencesOf(startTag, searchStart, 
                                                                          endTagIndex);
  int nestingLevel = otherOccurrencesOfStartTag.getNumElements();
  for(int i=0; i<nestingLevel; i++)
    endTagIndex = document.findFirstOccurrenceOf(endTag, endTagIndex+endTag.getLength());

  endTagIndex += endTag.getLength()-1; // return position of last character in end-tag
  return endTagIndex;
}

rsString rsKeyValueStringTree::getTagName(const rsString& document, const int openingBracketIndex)
{
  int tagEnd = openingBracketIndex+1;
  while( document.at(tagEnd) != '>' && document.at(tagEnd) != ' ' )
    tagEnd++;
  return document.getSubString(openingBracketIndex+1, tagEnd-1);
}

void rsKeyValueStringTree::createChildNodeFromXmlString(const rsString& xmlString)
{
  rsString childNodeName = getTagName(xmlString, 0);
  rsKeyValueStringTree *childNode = new rsKeyValueStringTree(childNodeName);
  childNode->setupFromXmlString(xmlString);
  hangInSubTree(childNode);
}
