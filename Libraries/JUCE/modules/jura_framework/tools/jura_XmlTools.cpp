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

//-------------------------------------------------------------------------------------------------

XmlElement* automatableModuleStateToXml(const AutomatableModule* device, 
  XmlElement* xmlElementToStartFrom)
{
  // the XmlElement which stores all the relevant state-information:
  XmlElement* xmlState;
  if( xmlElementToStartFrom == NULL )
    xmlState = new XmlElement(juce::String("AutomatableModuleState")); 
  else
    xmlState = xmlElementToStartFrom;

  // create an XmlElement which will be added as child element when there are any relevant 
  // controller-mappings to be stored
  XmlElement* xmlMapping = new XmlElement( juce::String("ControllerMapping") );

  int numParameters = device->getNumParameters();
  Parameter            *p;
  AutomatableParameter *ap;
  for(int i=0; i<numParameters; i++)
  {
    p  = device->getParameterByIndex(i);
    ap = dynamic_cast<AutomatableParameter*> (p);

    // check the pointer against NULL to be on the safe side for cases in which another thread 
    // changes the number of parameters while we are looping through them:
    if( ap != NULL )
    {
      // store only non-default mappings:
      if( ap->isInDefaultState() == false )
      {
        XmlElement* xmlParameterSetup = new XmlElement(juce::String(ap->getName()));
        xmlParameterSetup->setAttribute(juce::String("MidiCC"), ap->getAssignedMidiController());
        xmlParameterSetup->setAttribute(juce::String("Min"),    ap->getLowerAutomationLimit()  );
        xmlParameterSetup->setAttribute(juce::String("Max"),    ap->getUpperAutomationLimit()  );
        xmlMapping->addChildElement(xmlParameterSetup);
      }
    }
  }

  // when there is actually something stored in the xmlMapping element, add it as child element to
  // the xmlState, otherwise delete it (in the former case, the parent will take care for the 
  // deletion):
  if( xmlMapping->getNumChildElements() != 0 )
    xmlState->addChildElement(xmlMapping);
  else
    delete xmlMapping;

  // having stored the controller mappings for the automatable parameters, we now proceed to store 
  // their current values:
  for(int i=0; i<numParameters; i++)
  {
    p = device->getParameterByIndex(i);

    // check the pointer against NULL to be on the safe side for cases in which another thread 
    // changes the number of parameters while we are looping through them:
    if( p != NULL )
    {
      if( p->shouldBeSavedAndRecalled() && !p->isCurrentValueDefaultValue() )
      {
        // retrieve the name and store the value of the parameter
        if( p->isStringParameter() )
          xmlState->setAttribute( juce::String(p->getName()), juce::String(p->getStringValue()) );
        else
          xmlState->setAttribute( juce::String(p->getName()), juce::String(p->getValue()) );
      }
    }
  }

  return xmlState;
}

bool automatableModuleStateFromXml(AutomatableModule* device, const XmlElement &xmlState)
{
  bool success = true;

  // try to restore the values of the automatable parameters:
  juce::String name;
  int numParameters = device->getNumParameters();
  Parameter* p;


  double test;
  for(int i=0; i<numParameters; i++)
  {
    p = device->getParameterByIndex(i);

    // check the pointer against NULL to be on the safe side for cases in which another thread 
    // changes the number of parameters while we are looping through them:
    if( p != NULL )
    {
      // retrieve the name of the parameter:
      name = juce::String(p->getName());

      // look, if there is an attribute stored with this name, if so, set up the device with the
      // value, if not fall back to the default value:
      if( p->shouldBeSavedAndRecalled() )
      {
        if( p->isStringParameter() )
          p->setStringValue(xmlState.getStringAttribute(
            name, p->getDefaultStringValue()), true, true);
        else
        {
          // test:
          test =  xmlState.getDoubleAttribute(name, -10000);
          //if( test == -10000 )
          //  DEBUG_BREAK;

          p->setValue(xmlState.getDoubleAttribute(name, p->getDefaultValue()), true, true);
        }
      }
      //and the value
      //xmlState->setAttribute( juce::String(p->getName()), juce::String(p->getValue()) );
    }
  }

  // reset the controller mappings to the defaults:
  device->revertToDefaultMapping();

  // look for the controller mapping element - it is supposed to exist as child element in the
  // xmlState element:
  XmlElement* xmlMapping = xmlState.getChildByName(juce::String("ControllerMapping") );
  if( xmlMapping == NULL )
    return success;

  // it was not NULL so there are some controller mappings stored in this xmlState - we now 
  // retrieve them:
  numParameters = xmlMapping->getNumChildElements(); // we actually do not use this?
  int    midiCC = -1;
  char*  nameC;
  double min, max;
  //AutomatableParameter* p;

  // the control mapping settings for the individual parameters are stored as child elements in 
  // the xmlMapping element, so we loop through all those child elements now:
  forEachXmlChildElement(*xmlMapping, xmlParameterSetup)
  {
    name  = xmlParameterSetup->getTagName();
    nameC = toZeroTerminatedString(name);
    p     = device->getParameterByName(nameC);

    // check if a parameter in the device has this name (pointer is non-null in this case) and if 
    // so, restore its settings from the xmlParameterSetup:
    AutomatableParameter *ap = dynamic_cast<AutomatableParameter*> (p);
    if( ap != NULL )
    {
      midiCC = xmlParameterSetup->getIntAttribute(   juce::String("MidiCC"), -1);
      min    = xmlParameterSetup->getDoubleAttribute(juce::String("Min"),    ap->getMinValue() );
      max    = xmlParameterSetup->getDoubleAttribute(juce::String("Max"),    ap->getMaxValue() );    
      ap->assignMidiController(midiCC);
      ap->setLowerAutomationLimit(min);
      ap->setUpperAutomationLimit(max);
    }

    // we must take care to delete the c-string which is created by toZeroTerminatedString():
    delete[] nameC;
  }

  return success;
}
