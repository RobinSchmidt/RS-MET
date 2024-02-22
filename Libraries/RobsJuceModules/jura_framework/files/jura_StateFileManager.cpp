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
  if(result == true)
  {
    markFileAsClean(true);
    FileManager::setActiveFile(fileToLoad);   // new
  }
  notifyListeners();        // Why is this not inside the "if(result == true)" ?
  return result;
  // I think, we need to update the activeFileIndex. I think we need to call 
  // setActiveFileIfInList or setActiveFile
}

bool StateFileManager::saveToFile(const File& fileToSaveTo)
{
  bool result = saveStateToXmlFile(fileToSaveTo);
  if(result == true)
  {
    markFileAsClean(true);
    FileManager::setActiveFile(fileToSaveTo);   // new
  }
  notifyListeners();       // Why is this not inside the "if(result == true)" ?
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
      // I think, updating the file list is only necessarry when the new fileToLoadFrom resides in
      // a different directory than the filed that was previously loaded. Maybe put the update into
      // a conditional and then test, test, test


      delete xmlState;
      setStateName(fileToLoadFrom.getFileNameWithoutExtension(), true);

      //markStateAsClean();  
      // Why is this commented out? I think, it's because loadStateFromXmlFile is called only from
      // loadFile which calls markStateAsClean - but the method is public so client code could 
      // potentially call it...hmm...maybe we can make it private? ...OK...done

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
  if( xmlState == nullptr )
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
  //markStateAsClean();  // Why is this commented out? See comment in load...
  return true;
}


