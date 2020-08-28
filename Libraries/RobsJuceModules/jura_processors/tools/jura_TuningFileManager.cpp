//-------------------------------------------------------------------------------------------------
// construction/destruction:

TuningFileManager::TuningFileManager()
{
  wildcardPatterns = juce::String("*.tun;*.xml");
  //wildcardPatterns = juce::String("*.tun");
  defaultExtension = juce::String(".tun");

  // set the current directory for the tuning files:
  //setActiveDirectory(getApplicationDirectory() + juce::String("/TuningTables") );
  setActiveDirectory(getSupportDirectory() + "/TuningTables");
  theTable = nullptr;
}

TuningFileManager::~TuningFileManager()
{

}

//-------------------------------------------------------------------------------------------------
// others:

XmlElement* tuningTableStateToXml(TuningTable* tuningTable, XmlElement* xmlElementToStartFrom)
{
  // the XmlElement which stores all the relevant state-information:
  XmlElement* xmlState;
  if( xmlElementToStartFrom == NULL )
    xmlState = new XmlElement(juce::String("TuningTable")); 
  else
    xmlState = xmlElementToStartFrom;

  xmlState->setAttribute("Name",         juce::String(tuningTable->getName()) );
  xmlState->setAttribute("MasterTuneA4", tuningTable->getMasterTuneA4() );
  for(int i=0; i<128; i++)
    xmlState->setAttribute(juce::String(i), tuningTable->getFrequency(i) );

  return xmlState;
}

bool tuningTableStateFromXml(TuningTable* tuningTable, const XmlElement &xmlState)
{
  bool success = true;

  tuningTable->resetToDefaults();

  juce::String name = xmlState.getStringAttribute("Name", "12-TET");
  char* nameC = toZeroTerminatedString(name);
  tuningTable->setName(nameC);
  delete[] nameC;
  tuningTable->setMasterTuneA4(xmlState.getDoubleAttribute("MasterTuneA4", 440.0) );
  if( name == juce::String("12-TET") )
  {
    return true; 
    // leave it in default state and return (avoids roundoff by loading the default values from 
    // the XML)
  }

  /*
  for(int i=0; i<128; i++)
  {
  tuningTable->assignFrequency(i,      
  xmlState.getDoubleAttribute(juce::String(i), pitchToFreq((double) i)));
  }
  */

  double basicInterval    = xmlState.getDoubleAttribute("BasicInterval",            2.0);
  int    basicKeyDistance = xmlState.getIntAttribute(   "BasicIntervalKeyDistance", 12); 
  double numerator        = xmlState.getDoubleAttribute("Numerator",                1.0);
  double denominator      = xmlState.getDoubleAttribute("Denominator",              1.0);   
  double factor           = numerator/denominator;  

  for(int i=0; i<128; i++)
  {
    if( xmlState.hasAttribute("K" + juce::String(i)) )
      tuningTable->assignFrequency(i, factor * xmlState.getDoubleAttribute("K" + juce::String(i), 20.0));
    else
    {
      // if there is no frequency stored for this note number, we try to extrapolate by octave 
      // equivalence from some other assigned note number to this one:
      bool equivalentFound = false;;
      int  j = 1;
      int  d = j*basicKeyDistance;

      // check note-numbers below the current one i:
      while( (i-d) >= 0 )
      {
        if( xmlState.hasAttribute("K" + juce::String(i-d)) )
        {
          tuningTable->assignFrequency(i, 
            pow(basicInterval, j)*factor*xmlState.getDoubleAttribute("K" + juce::String(i-d), 20.0));
          equivalentFound = true;
          break;
        }
        j++;
        d += basicKeyDistance;
      }

      if( equivalentFound )
        continue;

      // check note-numbers above the current one i:
      j = 1;
      d = j*basicKeyDistance;
      while( (i+d) < 128 )
      {
        if( xmlState.hasAttribute("K" + juce::String(i+d)) )
        {
          tuningTable->assignFrequency(i, 
            pow(basicInterval, -j)*factor*xmlState.getDoubleAttribute("K" + juce::String(i+d), 20.0));
          equivalentFound = true;
          break;
        }
        j++;
        d += basicKeyDistance;
      }

      if( equivalentFound )
        continue;

      // nothing found to extrapolate from - use 20 Hz as fallback frequency:
      tuningTable->assignFrequency(i, 20.0);
    }
  }

  return success;
}

void TuningFileManager::assignTuningTable(rosic::TuningTable *newTable)
{
  theTable = newTable;
}

bool TuningFileManager::loadTuningFromFile(const juce::File &fileToLoadFrom)
{
  if( theTable == NULL )
    return false;

  if( fileToLoadFrom.hasFileExtension(juce::String("xml")) )
  {
    auto xmlTuning = getXmlFromFile(fileToLoadFrom);
    if( xmlTuning == NULL  )
      return false;
    tuningTableStateFromXml(theTable, *xmlTuning);
  }
  else if( fileToLoadFrom.hasFileExtension(juce::String("tun")) )
  {
    juce::String fullPath = fileToLoadFrom.getFullPathName();
    char*  pathC    = toZeroTerminatedString(fullPath);
    theTable->loadFromTunFile(pathC);
    delete[] pathC;
  }

  juce::String fileName = fileToLoadFrom.getFileName();
  char*  nameC    = toZeroTerminatedString(fileName);
  theTable->setName(nameC);
  delete[] nameC;

  return true;
}

//-------------------------------------------------------------------------------------------------
// FileManager overrides:

bool TuningFileManager::loadFile(const juce::File& fileToLoad)
{ 
  bool success = loadTuningFromFile(fileToLoad);
  notifyListeners();
  return success;
}

bool TuningFileManager::saveToFile(const juce::File& fileToSaveTo)
{
  jassertfalse;
  return false; // not yet functional
}













