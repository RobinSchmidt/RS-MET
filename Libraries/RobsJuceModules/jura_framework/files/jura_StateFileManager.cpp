//#include "rojue_StateFileManager.h"
//using namespace rojue;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

StateFileManager::StateFileManager() 
{
  wildcardPatterns = String("*.xml");
  defaultExtension = String(".xml");
}

StateFileManager::~StateFileManager()
{

}

//-------------------------------------------------------------------------------------------------
// FileManager overrides:

bool StateFileManager::loadFile(const File& fileToLoad)
{ 
  bool result = loadStateFromXmlFile(fileToLoad);
  if( result == true ) 
    markFileAsClean(true);    
  notifyListeners();
  return result;
}

bool StateFileManager::saveToFile(const File& fileToSaveTo)
{
  bool result = saveStateToXmlFile(fileToSaveTo);
  if( result == true ) 
    markFileAsClean(true);    
  notifyListeners();
  return result;
}

//-------------------------------------------------------------------------------------------------
// others:

bool StateFileManager::loadStateFromXmlFile(const File& fileToLoadFrom)
{
  if( fileToLoadFrom.existsAsFile() )
  {
    XmlDocument myDocument(fileToLoadFrom);
    //XmlElement* xmlState = myDocument.getDocumentElement(); // old
    XmlElement* xmlState = new XmlElement(*myDocument.getDocumentElement().get()); // new, preliminary
    if( xmlState != nullptr )
    {
      setStateFromXml(*xmlState, fileToLoadFrom.getFileNameWithoutExtension(), true);
      updateFileList();
      delete xmlState;
      setStateName(fileToLoadFrom.getFileNameWithoutExtension(), true);
      //markStateAsClean();
      return true;
    }
    else
      return false;
  }
  else
    return false;
}

bool StateFileManager::saveStateToXmlFile(const File& fileToSaveTo)
{
  // maybe we can use rojue::saveXmlToFile here?

  XmlElement* xmlState = getStateAsXml(fileToSaveTo.getFileNameWithoutExtension(), true);
  if( xmlState == NULL )
    return false;
  String myXmlDoc = xmlState->toString();
  if( fileToSaveTo.existsAsFile() )
  {
    File tmpFile = fileToSaveTo.getNonexistentSibling();
    tmpFile.create();
    tmpFile.appendText(myXmlDoc);
    fileToSaveTo.deleteFile();
    tmpFile.moveFileTo(fileToSaveTo);
  }
  else
  {
    fileToSaveTo.create();
    fileToSaveTo.appendText(myXmlDoc);
  }

  jassert( fileToSaveTo.existsAsFile() ); // file could not be created?

  updateFileList();
  delete xmlState;
  setStateName(fileToSaveTo.getFileNameWithoutExtension(), true);
  //markStateAsClean();
  return true;
}


