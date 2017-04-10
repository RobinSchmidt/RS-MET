#ifndef jura_XmlTools_h
#define jura_XmlTools_h

/** Returns an XmlElement from a file - returns NULL if the file doesn't exist or is not a
valid xml file. The caller is responsible for deleting the XmlElement. */
JUCE_API XmlElement* getXmlFromFile(const juce::File& fileToLoadFrom);

/** Returns an XmlElement from a file - returns NULL if the file doesn't exist or is not a
valid xml file. The caller is responsible for deleting the XmlElement. */
JUCE_API XmlElement* getXmlFromFile(const juce::String& fileNameToLoadFrom);

/** Saves the passed XmlElement to a file. If you set the final argument to false, it will 
overwrite an already existing file without asking, so be careful with that. */
JUCE_API bool saveXmlToFile(const XmlElement& xmlToSave, const juce::File& fileToSaveTo, 
  bool askForOverwrite = true);

/** Scans the passed 'parentElement' for child elements with tag name 'tagNameToLookFor' and
adds the pointers to these elements in the array 'results'. The return value will tell the number
of found/added child elements. */
JUCE_API int findChildElementsWithTagName(juce::Array<XmlElement*> &results, 
  const XmlElement& parentElement, const juce::String& tagNameToLookFor);

/** When the XmlElement passed in"xml" has multiple child-elements with the same name, this 
function returns the "index"-th one. The passed index should be zero for the first child-element 
with given name and then count upward. The function will return a NULL pointer when the number of 
child-elements with given name is smaller than the given index plus one. */
JUCE_API XmlElement* getChildElementByNameAndIndexAmongNameSakes(const XmlElement& xml, 
  const juce::String& name, int index);

#endif