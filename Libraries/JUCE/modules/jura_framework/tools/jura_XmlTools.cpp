XmlElement* stringToXml(const String& xmlStr)
{
  XmlDocument xmlDoc(xmlStr);
  return xmlDoc.getDocumentElement();
}

XmlElement* getXmlFromFile(const File &fileToLoadFrom)
{
  if( fileToLoadFrom.existsAsFile() )
  {
    XmlDocument myDocument(fileToLoadFrom);
    XmlElement *xml = myDocument.getDocumentElement();
    return xml;
  }
  else
    return NULL;
}

XmlElement* getXmlFromFile(const String &fileNameToLoadFrom)
{
  return getXmlFromFile(File(fileNameToLoadFrom));
}

bool saveXmlToFile(const XmlElement &xmlToSave, const File &fileToSaveTo, bool askForOverwrite)
{
  // turn the XmlElement thing into a text document:
  String myXmlDoc = xmlToSave.createDocument(String::empty);

  if( fileToSaveTo.existsAsFile() )
  {
    bool overwrite; 
    if( askForOverwrite == false )
      overwrite = true;
    else
    {
      overwrite = AlertWindow::showOkCancelBox(AlertWindow::WarningIcon, String("Overwrite?"), 
      String("File already exists - overwrite?") );
    }
    if( overwrite == false )
      return false;
    fileToSaveTo.deleteFile();
  }

  fileToSaveTo.create();
  fileToSaveTo.appendText(myXmlDoc);

  return true;
}

int findChildElementsWithTagName(Array<XmlElement*> &results, const XmlElement& parentElement,
  const String& tagNameToLookFor)
{
  int numResults = 0;
  //forEachXmlChildElementWithTagName(parentElement, child, tagNameToLookFor)
  //{
  //  results.add(child);
  //  numResults++;
  //}
  XmlElement* child = NULL;
  for(int c=0; c<parentElement.getNumChildElements(); c++)
  {
    child = parentElement.getChildElement(c);
    if( child->hasTagName(tagNameToLookFor) )
    {
      results.add(child);
      numResults++;
    }
  }
  return numResults;
}

XmlElement* getChildElementByNameAndIndexAmongNameSakes(const XmlElement& xml, 
  const juce::String& name, int index)
{
  if( index < 0 )
    return NULL;
  else if( index == 0 )
    return xml.getChildByName(name); 
  else
  {
    int i = 0;
    XmlElement* child = xml.getChildByName(name);
    while(child != NULL && i < index)
    {
      child = child->getNextElement();
      i++;
    }  
    return child;
  }
}
