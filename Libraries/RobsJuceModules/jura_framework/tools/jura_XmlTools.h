#ifndef jura_XmlTools_h
#define jura_XmlTools_h

// todo: instead of returning raw pointers that the caller has to delete, return 
// std::unique_ptr<XmlElement> because that's what juce::XmlDocument::getDocumentElement returns 
// nowadays: https://docs.juce.com/master/classXmlDocument.html
// we could perhaps do something dirty by retrieving the managed object from the std::unique_ptr
// and return the raw-pointer from our function anyway, defeating the purpose of std::unique_ptr 
// but this should then perhaps be considered deprecated legacy-compatibility code:
// https://en.cppreference.com/w/cpp/memory/unique_ptr/release

/** Converts a string that represents an xml document to an actual XmlElement object. The caller is
responsible for deleting the object eventually. Convenience function to save the user from going 
through the XmlDocument class. */
JUCE_API XmlElement* stringToXml(const String& xmlStr);

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

/** When the XmlElement passed as "xml" has multiple child-elements with the same name, this 
function returns the "index"-th one. The passed index should be zero for the first child-element 
with given name and then count upward. The function will return a NULL pointer when the number of 
child-elements with given name is smaller than the given index plus one. */
JUCE_API XmlElement* getChildElementByNameAndIndexAmongNameSakes(const XmlElement& xml, 
  const juce::String& name, int index);

/** Takes a std::map of std::strings (for key and value) and for each key/value pair in the map,
adds a cooresponding attribute to the given XmlElement. */
JUCE_API void addAttributesFromMap(XmlElement& xml, std::map<std::string, std::string>& map);

/** Returns a std::map of std::string containing all attributes in the given XmlElement as 
key/value pairs. */
JUCE_API std::map<std::string, std::string> getAttributesAsMap(const XmlElement& xml);


//JUCE_API bool isValidXmlAttributeName(const juce::String& s);

//XmlElement::isValidXmlName

#endif