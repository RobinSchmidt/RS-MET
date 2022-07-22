XmlElement* stringToXml(const String& xmlStr)
{
  XmlDocument xmlDoc(xmlStr);

  //return xmlDoc.getDocumentElement();  // old
  return new XmlElement(*xmlDoc.getDocumentElement().get()); 
  // new: returns a deep copy - preliminary, todo: return the std::unique_ptr now returned from 
  // getDocumentElement directly
}

XmlElement* getXmlFromFile(const File &fileToLoadFrom)
{
  if( fileToLoadFrom.existsAsFile() )
  {
    XmlDocument myDocument(fileToLoadFrom);
    //XmlElement *xml = myDocument.getDocumentElement(); // old
    XmlElement *xml = new XmlElement(*myDocument.getDocumentElement().get()); // new, preliminary
    return xml;
  }
  else
    return nullptr;
}

XmlElement* getXmlFromFile(const String &fileNameToLoadFrom)
{
  return getXmlFromFile(File(fileNameToLoadFrom));
}

bool saveXmlToFile(const XmlElement &xmlToSave, const File &fileToSaveTo, bool askForOverwrite)
{
  // turn the XmlElement thing into a text document:
  //String myXmlDoc = xmlToSave.createDocument(String());
  String myXmlDoc = xmlToSave.toString();

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
  // Maybe instead of appendText, we should use replaceWithText?
  // https://docs.juce.com/master/classFile.html#ace1786a60bfc469bfa75b3638923e163

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

void addAttributesFromMap(XmlElement& xml, std::map<std::string, std::string>& map)
{
  std::map<std::string, std::string>::iterator it = map.begin();
  while(it != map.end()) {
    std::string key   = it->first;
    std::string value = it->second;
    xml.setAttribute(juce::String(key), juce::String(value));
    it++;
  }
  // see here:
  // https://thispointer.com/how-to-iterate-over-a-map-in-c/
}

std::map<std::string, std::string> getAttributesAsMap(const XmlElement& xml)
{
  std::map<std::string, std::string> map;
  for(int i = 0; i < xml.getNumAttributes(); i++) {
    std::string name  = xml.getAttributeName(i).toStdString();
    std::string value = xml.getAttributeValue(i).toStdString();
    map.emplace(name, value);
  }
  return map;
}

/*
bool isValidXmlAttributeName(const juce::String& s)
{
  int N = s.length();
  if(N == 0)
    return false;

  // this works - let's keep it simple!
  typedef rosic::rsString S;
  if(!(S::isLetter(s[0]) || s[0] == '_'))            // allow letters or underscores for 1st char
    return false;
  for(int i = 1; i < N; i++)
    if(!(S::isLetterOrDigit(s[i]) || s[i] == '_'))   // allow letters, digits or underscores
      return false;

  return true;
}
// superfluous XmlElement::isValidXmlName should work for attributes also
// see here
// https://docstore.mik.ua/orelly/xml/xmlnut/ch02_04.htm
*/