rsParameterSetupBase::rsParameterSetupBase(rsAutomatableWidget* widgetToSetup,
  MetaParameterManager* metaManagerToUse)
  : rsDeletionRequester(widgetToSetup), widget(widgetToSetup), metaManager(metaManagerToUse)
{
  addWidget( closeButton = new RClickButtonNotifyOnMouseUp(RButton::CLOSE) );
  closeButton->setClickingTogglesState(false);
  closeButton->addRButtonListener(this);

  ColourSchemeComponent* csc = dynamic_cast<ColourSchemeComponent*>
    (widget->getWrappedWidget()->getParentComponent());
  if(csc){
    editorColourScheme = csc->getEditorColourScheme();
    widgetColourScheme = csc->getWidgetColourScheme();
    plotColourScheme   = csc->getPlotColourScheme();
  }
}

//=================================================================================================

rsMetaMapEditor::rsMetaMapEditor(MetaParameterManager* metaManagerToUse)
{
  metaManager = metaManagerToUse;
  param = nullptr;

  //setAlwaysOnTop(true); // hmm - this doesn't seem to help - it still gets obscured when moving the
                        // parameter slider
}

rsMetaMapEditor::~rsMetaMapEditor()
{
  if(param)
    param->deRegisterParameterObserver(this);
}

void rsMetaMapEditor::setParameterToControl(MetaControlledParameter* p) 
{ 
  if(param) 
    param->deRegisterParameterObserver(this);
  param = p; 
  if(param)
    param->registerParameterObserver(this);
}

MetaParameter* rsMetaMapEditor::getAttachedMetaParameter()
{
  if(param)
    return param->getAttachedMeta();
  return nullptr;
}

void rsMetaMapEditor::paint(Graphics& g)
{
  rsNodeBasedFunctionEditor::paint(g);

  // new - shows line always:
  if(param != nullptr)
  {
    double nv = param->getNormalizedValue();
    float pixelX = (float) xyMapper.mapX(nv);
    g.drawLine(pixelX, 0.f, pixelX, (float) getHeight(), 1.f);
  }

  // old - shows line onle when a meta is attached:
  /*
  MetaParameter* mp = getAttachedMetaParameter();
  if(mp != nullptr)
  {
    double mv = mp->getMetaValue();
    float pixelX = (float) xyMapper.mapX(mv);
    g.drawLine(pixelX, 0.f, pixelX, (float) getHeight(), 1.f);
  }
  */
}

void rsMetaMapEditor::parameterChanged(Parameter* p)
{
  if(p == param)
    repaintOnMessageThread(); // to update position of vertical line
  else
  {
    updateParameter();
    rsNodeBasedFunctionEditor::parameterChanged(p);
  }
}

int rsMetaMapEditor::addNode(double pixelX, double pixelY)
{
  int result = rsNodeBasedFunctionEditor::addNode(pixelX, pixelY);
  updateParameter();
  return result;
}

bool rsMetaMapEditor::removeNode(int index)
{
  bool result = rsNodeBasedFunctionEditor::removeNode(index);
  updateParameter();
  return result;
}

int rsMetaMapEditor::moveNodeTo(int index, int pixelX, int pixelY)
{
  int result = rsNodeBasedFunctionEditor::moveNodeTo(index, pixelX, pixelY);
  updateParameter();
  return result;
}

rsNodeBasedFunctionEditor::NodeParameterSet* rsMetaMapEditor::addNodeParameters(
  rsDraggableNode* node)
{
  rsNodeBasedFunctionEditor::NodeParameterSet* params = 
    rsNodeBasedFunctionEditor::addNodeParameters(node);
  setupNodeParameterObservation(node, params);
  return params;
}

void rsMetaMapEditor::setupNodeParameterObservation(
  rsDraggableNode* node, NodeParameterSet* params)
{
  params->x          ->registerParameterObserver(this);
  params->y          ->registerParameterObserver(this);
  params->shapeType  ->registerParameterObserver(this);
  params->shapeAmount->registerParameterObserver(this);
  // register, so we can call updateParameter() on changes
}

void rsMetaMapEditor::updateParameter()
{
  if(param)
  {
    //double metaValue = param->getNormalizedValue();
    //param->setFromMetaValue(metaValue, true, true);
    // this doesn't work because in MetaControlledParameter::setNormalizedValue there's a check
    // if(normalizedValue == newNormalizedValue) which lets the function return early

    param->metaMapChanged();
  }
  // we need to call this also when the node position is changed via the sliders and/or the node
  // cuver is changed
}

//=================================================================================================

rsAutomationSetup::rsAutomationSetup(rsAutomatableWidget* widgetToAutomate,
  MetaParameterManager* metaManagerToUse)
  : rsParameterSetupBase(widgetToAutomate, metaManagerToUse)
{
  closeButton->setDescription("Closes the automation setup window");
  createWidgets();
  setAlwaysOnTop(true); // should be stay on top, if parameter slider is moved
  metaMapEditor->registerObserver(this);
  setSize(300, 250);
}

rsAutomationSetup::~rsAutomationSetup() 
{
  metaMapEditor->deRegisterObserver(this);
  //delete boxMetaAttach; // why? shouldn't that happen automatically? ..ok, no memleak warning
}

void rsAutomationSetup::rButtonClicked(RButton *b)
{
  if(b == closeButton)
    requestDeletion();  // maybe move to baseclass
}

void rsAutomationSetup::rComboBoxChanged(RComboBox* cb)
{
  if(cb == boxMetaAttach)
  {
    // ...retrieve selected item and attach meta accordingly
  }
  else if(cb == boxShapeType)
  {
    // set shape type fo selecetd node
  }
}

void rsAutomationSetup::resized()
{
  int d   = sliderDistance;
  int x   = d;
  int y   = d;
  int w   = getWidth();
  int h   = getHeight();
  int sh  = sliderHeight;
  int inc = sh+d;

  closeButton->setBounds(w-16, 0, 16, 16);
  //automationLabel->setBounds(x, y, w-8-16, sh); y += inc; 

  y = closeButton->getBottom() + d;

  // top row:
  w = getWidth()/2 - 2*d;
  boxMetaAttach->setBounds(  x, y, w, sh);
  x = getWidth()/2 + d;
  sliderSmoothing->setBounds(x, y, w, sh);

  // bottom row:
  x = d;
  y = getHeight() - 2*sh - d + 2;
  w = getWidth()/2 - 2*d;
  sliderNodeX->setBounds(x, y,      w, sh);
  sliderNodeY->setBounds(x, y+sh-2, w, sh);
  x = getWidth()/2 + d;
  boxShapeType    ->setBounds(x, y,      w, sh);
  sliderShapeParam->setBounds(x, y+sh-2, w, sh);

  // plot:
  y = sliderSmoothing->getBottom() + d;
  w = getWidth();
  h = sliderNodeX->getY()-y-d;
  metaMapEditor->setBounds(0, y, w, h);
}

void rsAutomationSetup::setVisible(bool shouldBeVisible)
{
  if(shouldBeVisible)
    metaMapEditor->updateDraggableNodesArray(); 
    // to keep them in sync with the corresponding nodes in RAPT::rsNodeBasedFunction
  rsParameterSetupBase::setVisible(shouldBeVisible);
}

void rsAutomationSetup::nodeWasAdded(rsNodeEditor* editor, int nodeIndex)
{
  // not needed - maybe make it optional to override this
}

void rsAutomationSetup::nodeWillBeRemoved(rsNodeEditor* editor, int nodeIndex)
{
  // dito
}

void rsAutomationSetup::nodeWasSelected(rsNodeEditor* editor, int nodeIndex)
{
  assignNodeParameterWidgets(nodeIndex);
}

void rsAutomationSetup::assignNodeParameterWidgets(int i)
{
  if(i != -1) {
    rsDraggableNode* node = metaMapEditor->getNode(i);
    sliderNodeX     ->assignParameter(node->getParameterX());
    sliderNodeY     ->assignParameter(node->getParameterY());
    boxShapeType    ->assignParameter(node->getNodeParameter(0));
    sliderShapeParam->assignParameter(node->getNodeParameter(1));
    // we get an access violation when selecting an existing node, but none when inserting a new node
  }
  else {
    sliderNodeX     ->assignParameter(nullptr);
    sliderNodeY     ->assignParameter(nullptr);
    boxShapeType    ->assignParameter(nullptr);
    sliderShapeParam->assignParameter(nullptr);
  }
  updateWidgetVisibility();
}

void rsAutomationSetup::createWidgets()
{
  addWidget(metaMapEditor = new rsMetaMapEditor(metaManager));
  metaMapEditor->setDescription("Mapping between meta-parameter and parameter");
  metaMapEditor->setValueRange(0, 1, 0, 1);
  metaMapEditor->setClipCoordinatesToRange(true);

  MetaControlledParameter* mcp = widget->getMetaControlledParameter();
  if(mcp != nullptr)
  {
    metaMapEditor->setParameterToControl(mcp);
    metaMapEditor->setFunctionToEdit(mcp->getMetaMapper());
    metaMapEditor->setMutexToUse(mcp->getUsedMutex());
  }

  addWidget(boxMetaAttach = new RNamedComboBox("", "Meta:"));
  boxMetaAttach->setDescription("Select meta parameter to attach");
  // todo: fill the box

  addWidget(sliderSmoothing = new RSlider("Smoothing"));
  sliderSmoothing->setDescription("Smoothing time in milliseconds");

  addWidget(sliderNodeX = new RSlider("X"));
  sliderNodeX->setDescription("X-coordinate of selected node");

  addWidget(sliderNodeY = new RSlider("Y"));
  sliderNodeY->setDescription("Y-coordinate of selected node");

  addWidget(boxShapeType = new RNamedComboBox("", "Shape:"));
  boxShapeType->setDescription("Shape of segment approaching selected node");
  // fill the box: left/right/nearest neighbour, linear, cubic, etc.

  addWidget(sliderShapeParam = new RSlider("Shape Parameter"));
  sliderShapeParam->setDescription("Parameter of curve shape");

  updateWidgetVisibility();
}

void rsAutomationSetup::updateWidgetVisibility()
{
  sliderNodeX->setVisible(sliderNodeX->hasAssignedParameter());
  sliderNodeY->setVisible(sliderNodeY->hasAssignedParameter());
  sliderShapeParam->setVisible(sliderShapeParam->hasAssignedParameter());
  boxShapeType->setVisible(boxShapeType->hasAssignedParameter());
}

//=================================================================================================

rsModulationSetup::rsModulationSetup(rsAutomatableWidget* widgetToModulate, 
  MetaParameterManager* metaManagerToUse)
  : rsParameterSetupBase(widgetToModulate, metaManagerToUse)
{
  closeButton->setDescription("Closes the modulation setup window");

  addWidget( modulationsLabel = 
    new RTextField(widgetToModulate->getParameterName() + " Modulations") );
  modulationsLabel->setNoBackgroundAndOutline(true);
  modulationsLabel->setDescription("Modulation setup");

  addWidget( addButton = new RButton("Add") );
  addButton->setDescription(juce::String("Adds a new modulation connection"));
  addButton->setClickingTogglesState(false);
  addButton->addRButtonListener(this);

  addWidget( removeButton = new RButton("Remove") );
  removeButton->setDescription(juce::String("Removes an existing modulation connection"));
  removeButton->setClickingTogglesState(false);
  removeButton->addRButtonListener(this);

  addWidget( clipMinField = new RLabeledTextEntryField("Min:") );
  clipMinField->setEntryFieldText(String(getClipMin()));
  clipMinField->setDescription(juce::String("Range minimum for modulated value"));
  clipMinField->setLabelWidth(36);
  clipMinField->getTextEntryField()->registerTextEntryFieldObserver(this);

  addWidget( clipMaxField = new RLabeledTextEntryField("Max:") );
  clipMaxField->setEntryFieldText(String(getClipMax()));
  clipMaxField->setDescription(juce::String("Range maximum for modulated value"));
  clipMaxField->setLabelWidth(36);
  clipMaxField->getTextEntryField()->registerTextEntryFieldObserver(this);

  updateConnectionWidgetsArray();

  setAlwaysOnTop(true); // should not be hidden by main gui window

  //modManager = nullptr;
  ModulatableParameter* mp = widgetToModulate->getModulatableParameter();
  if(mp != nullptr)
    mp->registerModulationTargetObserver(this);
}

rsModulationSetup::~rsModulationSetup()
{
  ModulatableParameter* mp = widget->getModulatableParameter();
  if(mp != nullptr)
    mp->deRegisterModulationTargetObserver(this);

  delete connectableSourcesPopUp;
  delete removableSourcesPopUp;
}

//void rsModulationSetup::paint(Graphics& g)
//{
//  ColourSchemeComponent::paint(g);
//}

void rsModulationSetup::resized()
{
  int d   = sliderDistance;
  int x   = d;
  int y   = d;
  int w   = getWidth();
  int h   = getHeight();
  int sh  = sliderHeight;
  int inc = sh+d;

  closeButton->setBounds(w-16, 0, 16, 16);
  modulationsLabel->setBounds(x, y, w-8-16, sh); y += inc; 

  for(int i = 0; i < size(connectionWidgets); i++) {
    connectionWidgets[i]->setBounds(x, y, w-8, sh); y += inc; }

  y = h - sh - d;
  addButton->setBounds(x, y, 40, 16);
  x = addButton->getRight() + d;
  removeButton->setBounds(x, y, 60, 16);

  y -= inc;
  x  = 0;
  clipMinField->setBounds(x+d, y, w/2-2*d, 16);
  x = w/2;
  clipMaxField->setBounds(x+d, y, w/2-2*d, 16);
}

void rsModulationSetup::rButtonClicked(RButton *button)
{
  if(button == closeButton)
    requestDeletion();
  else if(button == addButton)
    showConnectableSourcesPopUp();
  else if(button == removeButton)
    showRemovableSourcesPopUp();
  for(int i = 0; i < size(connectionWidgets); i++){
    if(button == connectionWidgets[i]->removeButton)
      removeConnection(i);
  }
}

void rsModulationSetup::rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged)
{
  int id;
  if(menuThatHasChanged == connectableSourcesPopUp)
  {
    id = connectableSourcesPopUp->getSelectedIdentifier();
    if(id > 0)
      addConnection(id-1);
  }
  else if(menuThatHasChanged == removableSourcesPopUp)
  {
    id = removableSourcesPopUp->getSelectedIdentifier();
    if(id > 0)
      removeConnection(id-1);
  }
}

void rsModulationSetup::textChanged(RTextEntryField *rTextEntryFieldThatHasChanged)
{
  //jassertfalse;
  // this doesn't work yet: String::getDoubleValue doesn't parse -inf correctly - it returns inf

  if(rTextEntryFieldThatHasChanged == clipMinField->getTextEntryField())
    setClipMin(toDouble(clipMinField->getTextEntryField()->getText()));
  else if(rTextEntryFieldThatHasChanged == clipMaxField->getTextEntryField())
    setClipMax(toDouble(clipMaxField->getTextEntryField()->getText()));
}

void rsModulationSetup::modulationsChanged()
{
  updateConnectionWidgetsArray();
}

void rsModulationSetup::addConnection(int index)
{
  ModulatableParameter* mp = widget->getModulatableParameter();
  if(mp != nullptr)
  {
    std::vector<ModulationSource*> sources = mp->getDisconnectedSources();
    mp->addModulationSource(sources[index]);
    //addWidgetsForConnection(mp->getConnectionTo(sources[index]));
    // not necessary anymore - we do now a full update in the modulationsChanged callback
  }
}

void rsModulationSetup::removeConnection(int index)
{
  ModulatableParameter* mp = widget->getModulatableParameter();
  if(mp != nullptr)
  {
    //removeWidgetsForConnection(index); // not necessary anymore, see comment in addConnection
    connectionWidgets[index]->depthSlider->assignParameter(nullptr);
    std::vector<ModulationSource*> sources = mp->getConnectedSources();
    mp->removeModulationSource(sources[index]);
    updateSize();
  }
}

void rsModulationSetup::updateConnectionWidgetsArray()
{
  clearConnectionWidgets();
  ModulatableParameter* mp = widget->getModulatableParameter();
  if(mp != nullptr)
  {
    std::vector<ModulationConnection*> connections = mp->getConnections();
    for(int i = 0; i < size(connections); i++)
    {
      if(!hasSlider(connections[i]->getDepthParameter()))
        addWidgetsForConnection(connections[i]);
    }
  }
  updateSize();
}

void rsModulationSetup::showConnectableSourcesPopUp()
{
  // create popup, if necessary:
  if(connectableSourcesPopUp == nullptr)
  {
    connectableSourcesPopUp = new RPopUpMenu(this); // maybe attach to the addButton instead of this?
    connectableSourcesPopUp->registerPopUpMenuObserver(this);
    connectableSourcesPopUp->setDismissOnFocusLoss(true);
  }

  // populate it:
  connectableSourcesPopUp->clear();
  ModulatableParameter* mp = widget->getModulatableParameter();
  if(mp != nullptr)
  {
    std::vector<ModulationSource*> sources = mp->getDisconnectedSources();
    for(int i = 0; i < size(sources); i++)
    {
      juce::String name = sources[i]->getModulationSourceDisplayName();
      if(sources[i]->hasConnectedTargets())
        name += "*";
      connectableSourcesPopUp->addItem(i+1, name);     // +1 bcs 0 is not allowed for the id
    }
  }

  // show it:
  int w = connectableSourcesPopUp->getRequiredWidth(true);
  int h = connectableSourcesPopUp->getRequiredHeight(true);
  connectableSourcesPopUp->show(true, RPopUpComponent::BELOW, w, h); // showModally = true
}

void rsModulationSetup::showRemovableSourcesPopUp()
{
  // lots of code duplication from showConnectableSourcesPopUp - can this be refactored?

  if(removableSourcesPopUp == nullptr)
  {
    removableSourcesPopUp = new RPopUpMenu(this);
    removableSourcesPopUp->registerPopUpMenuObserver(this);
    removableSourcesPopUp->setDismissOnFocusLoss(true);
  }

  removableSourcesPopUp->clear();
  ModulatableParameter* mp = widget->getModulatableParameter();
  if(mp != nullptr)
  {
    std::vector<ModulationSource*> sources = mp->getConnectedSources();
    for(int i = 0; i < size(sources); i++)
    {
      juce::String name = sources[i]->getModulationSourceDisplayName();
      removableSourcesPopUp->addItem(i+1, name);
    }
  }

  int w = removableSourcesPopUp->getRequiredWidth(true);
  int h = removableSourcesPopUp->getRequiredHeight(true);
  removableSourcesPopUp->show(true, RPopUpComponent::BELOW, w, h);
}

bool rsModulationSetup::hasSlider(MetaControlledParameter* p)
{
  for(int i = 0; i < size(connectionWidgets); i++)
  {
    Parameter* ps = connectionWidgets[i]->depthSlider->getAssignedParameter();
    if(ps == p)
      return true;
  }
  return false;
}

void rsModulationSetup::addWidgetsForConnection(ModulationConnection* c)
{
  rsModulationConnectionWidget* w = new rsModulationConnectionWidget(c, this);
  connectionWidgets.push_back(w);
  w->depthSlider->assignParameter(c->getDepthParameter());
  w->depthSlider->setSliderName(c->getSource()->getModulationSourceDisplayName());
  w->removeButton->addRButtonListener(this);
  addWidget(w);
  updateSize();
}

/*
void rsModulationSetup::removeWidgetsForConnection(int i)
{
  connectionWidgets[i]->removeButton->removeRButtonListener(this);
  deleteObject(connectionWidgets[i]);              // mark for later deletion
  removeWidget(connectionWidgets[i], true, false); // false, to not delete it immediately
  remove(connectionWidgets, i);
}
*/

void rsModulationSetup::clearConnectionWidgets()
{
  for(int i = 0; i < size(connectionWidgets); i++)
  {  
    // maybe factor out these 3 lines into a function, they appear also in 
    // removeWidgetsForConnection:
    connectionWidgets[i]->removeButton->removeRButtonListener(this);
    deleteObject(connectionWidgets[i]);              // mark for later deletion
    removeWidget(connectionWidgets[i], true, false); // false, to not delete it immediately
  }
  connectionWidgets.clear();
}

void rsModulationSetup::updateSize()
{
  int width  = 250;  // maybe we should use the widget's width...but maybe not
  int height = 100;  // preliminary

  height  = (sliderHeight+sliderDistance) * size(connectionWidgets);
  height += 68;

  setSize(width, height); 
  resized(); // needed during development - might be redundant when finished
}

void rsModulationSetup::setClipMin(double newMin) 
{ 
  widget->getModulatableParameter()->setModulationRangeMin(newMin); 
}

void rsModulationSetup::setClipMax(double newMax) 
{ 
  widget->getModulatableParameter()->setModulationRangeMax(newMax); 
}

double rsModulationSetup::getClipMin()
{ 
  return widget->getModulatableParameter()->getModulationRangeMin(); 
}

double rsModulationSetup::getClipMax()
{ 
  return widget->getModulatableParameter()->getModulationRangeMax(); 
}

//=================================================================================================

rsAutomatableWidget::rsAutomatableWidget(RWidget *widgetToWrap)
{
  wrappedWidget = widgetToWrap;
}

rsAutomatableWidget::~rsAutomatableWidget()
{
  delete rightClickPopUp;
  delete modSetup;
  delete metaSetup;
}

bool rsAutomatableWidget::isPopUpOpen()
{
  return popUpIsOpen;

  // The code below doesn't work because apparently, when the user clicks on a widget while the
  // popup is open, the popup gets closed first and only after that, the mouseDown callback of the
  // widget is received, so isPopUpOpen would always return false in the mouseDown callback. The
  // desired behavior is that one right-click on the widget opens the popup and a second click
  // closes it. This behavior now requires that we maintain a popUpIsOpen flag here and that flag
  // should be set in the mouseDown method of the widget. Without that, a second right-click would
  // make the menu disappear for a fraction of a second and immediately reappear.

  //if(rightClickPopUp == nullptr)
  //  return false;
  //else
  //  return rightClickPopUp->isOpen();
}

void rsAutomatableWidget::rPopUpMenuChanged(RPopUpMenu* menu)
{
  if(menu != rightClickPopUp)
    return;

  int selectedIdentifier = rightClickPopUp->getSelectedIdentifier();

  if(selectedIdentifier == AUTOMATION_SETUP) { showAutomationSetup(); return; }
  if(selectedIdentifier == MODULATION_SETUP) { showModulationSetup(); return; }

  MetaControlledParameter* mcp = getMetaControlledParameter();
  if(mcp != nullptr)
  {
    switch(selectedIdentifier)
    {
    case META_ATTACH: mcp->attachToMetaParameter(
      (int)wrappedWidget->openModalNumberEntryField(mcp->getMetaParameterIndex()), false); break;
    case META_ATTACH_FLAT: mcp->attachToMetaParameter(
      (int)wrappedWidget->openModalNumberEntryField(mcp->getMetaParameterIndex()), true ); break;
    case META_DETACH: mcp->detachFromMetaParameter();    break;
    }
  }

  AutomatableParameter* ap = getAutomatableParameter();
  if(ap != nullptr)
  {
    if(selectedIdentifier != MIDI_LEARN)
      ap->switchIntoMidiLearnMode(false); // turn off, if currently on and something else is selected
    switch(selectedIdentifier)
    {
    case MIDI_LEARN:  ap->switchIntoMidiLearnMode();             break;
    case MIDI_ASSIGN:
    {
      int result = (int)wrappedWidget->openModalNumberEntryField(ap->getAssignedMidiController());
      result = (int)RAPT::rsClip(result, 0, 127);
      ap->assignMidiController(result);
    } break;
    case MIDI_MIN:    ap->setLowerAutomationLimit(ap->getValue());   break;
    case MIDI_MAX:    ap->setUpperAutomationLimit(ap->getValue());   break;
    case MIDI_REVERT: ap->revertToDefaults(false, false, false);     break;
    }
  }
}

void rsAutomatableWidget::rPopUpMenuDismissed(RPopUpMenu* menuThatwasDismissed)
{
  popUpIsOpen = false;
}

void rsAutomatableWidget::updatePopUpMenu()
{
  if(rightClickPopUp == nullptr)
  {
    // popup used the 1st time - we need to create it:
    rightClickPopUp = new RPopUpMenu(wrappedWidget);
    rightClickPopUp->registerPopUpMenuObserver(this);
    rightClickPopUp->setDismissOnFocusLoss(true);
    wrappedWidget->addChildWidget(rightClickPopUp, false, false);
  }
  rightClickPopUp->clear();
  addPopUpMenuItems();
}

void rsAutomatableWidget::addPopUpMenuItems()
{
  addPopUpMetaItems();
  addPopUpMidiItems();
  addPopUpModulationItems();
}

void rsAutomatableWidget::addPopUpMidiItems()
{
  AutomatableParameter* ap = getAutomatableParameter();
  if(ap != nullptr)
  {
    // prepare some strings for the popup menu:
    int cc = ap->getAssignedMidiController();
    String ccString;
    if(cc > -1)
      ccString = "(currently CC" + String(cc) + ")";
    else
      ccString = "(currently none)";

    int defaultCc = ap->getDefaultMidiController();
    String defaultString;
    if(defaultCc > -1)
      defaultString = "CC" + String(defaultCc);
    else
      defaultString = "none";
    String minString = wrappedWidget->stringConversionFunction(ap->getLowerAutomationLimit());
    String maxString = wrappedWidget->stringConversionFunction(ap->getUpperAutomationLimit());

    rightClickPopUp->addItem(MIDI_LEARN,  "MIDI learn " + ccString);
    rightClickPopUp->addItem(MIDI_ASSIGN, "MIDI assign");
    rightClickPopUp->addItem(MIDI_MIN,    "use value as lower limit (currently " + minString + String(")"));
    rightClickPopUp->addItem(MIDI_MAX,    "use value as upper limit (currently " + maxString + String(")"));
    rightClickPopUp->addItem(MIDI_REVERT, "revert MIDI mapping to defaults");
  }
}

void rsAutomatableWidget::addPopUpMetaItems()
{
  MetaControlledParameter* mcp = getMetaControlledParameter();
  if(mcp != nullptr)
  {
    int mi = mcp->getMetaParameterIndex();
    String miString;
    if(mi > -1)
      miString = "(currently " + mcp->getMetaParameterName() + ")";
    else
      miString = "(currently none)";

    rightClickPopUp->addItem(META_ATTACH,      "Meta attach " + miString);
    rightClickPopUp->addItem(META_DETACH,      "Meta detach");
    rightClickPopUp->addItem(META_ATTACH_FLAT, "Meta attach (flat)");
    rightClickPopUp->addItem(AUTOMATION_SETUP, "Automation setup");
  }
}

void rsAutomatableWidget::addPopUpModulationItems()
{
  ModulatableParameter* mp = getModulatableParameter();
  if(mp != nullptr)
    rightClickPopUp->addItem(MODULATION_SETUP, "Modulation setup");
}

void rsAutomatableWidget::openRightClickPopupMenu()
{
  updatePopUpMenu();
  int w = jmax(wrappedWidget->getWidth(), rightClickPopUp->getRequiredWidth(true));
  int h = jmin(200,                       rightClickPopUp->getRequiredHeight(true));
  //rightClickPopUp->show(false, RPopUpComponent::BELOW, w, h); // showModally = false
  rightClickPopUp->show(true, RPopUpComponent::BELOW, w, h); // showModally = true
  // If we don't show it modally (1st parameter = true), it will be immediately dismissed
  // after opening (so it appears as if it doesn't open at all). We could avoid it by calling
  // setDismissOnFocusLoss(false) in our constructor, but then it will stay open all the time
  // until we choose some option.

  popUpIsOpen = true;
}

void rsAutomatableWidget::closePopUp()
{
  if(rightClickPopUp != nullptr)
    rightClickPopUp->dismiss();
  popUpIsOpen = false;
}

void rsAutomatableWidget::showAutomationSetup()
{
  int wh = wrappedWidget->getHeight();       // widget height
  int x  = wrappedWidget->getScreenX();
  int y  = wrappedWidget->getScreenY() + wh; // preliminary

  if(metaSetup == nullptr)
    metaSetup = new rsAutomationSetup(this, getMetaParameterManager());
  metaSetup->setTopLeftPosition(x, y);
  metaSetup->addToDesktop(ComponentPeer::windowHasDropShadow | ComponentPeer::windowIsTemporary);
  metaSetup->setVisible(true);
  metaSetup->toFront(true);
}
// todo: get rid of code duplication
void rsAutomatableWidget::showModulationSetup()
{
  int wh = wrappedWidget->getHeight();       // widget height
  int x  = wrappedWidget->getScreenX();
  int y  = wrappedWidget->getScreenY() + wh; // preliminary

  if(modSetup == nullptr)
    modSetup = new rsModulationSetup(this, getMetaParameterManager());
  modSetup->setTopLeftPosition(x, y);
  modSetup->addToDesktop(ComponentPeer::windowHasDropShadow | ComponentPeer::windowIsTemporary);
  modSetup->setVisible(true);
  modSetup->toFront(true);
}

void rsAutomatableWidget::deleteObject(rsDeletionRequester* objectToDelete)
{
  if(objectToDelete == modSetup)
    modSetup->setVisible(false);
  else if(objectToDelete == metaSetup)
    metaSetup->setVisible(false);
  else
    jassertfalse;
}

AutomatableParameter* rsAutomatableWidget::getAutomatableParameter()
{
  return dynamic_cast<AutomatableParameter*> (wrappedWidget->assignedParameter);
}

MetaControlledParameter* rsAutomatableWidget::getMetaControlledParameter()
{
  return dynamic_cast<MetaControlledParameter*> (wrappedWidget->assignedParameter);
}

ModulatableParameter* rsAutomatableWidget::getModulatableParameter()
{
  return dynamic_cast<ModulatableParameter*> (wrappedWidget->assignedParameter);
}

String rsAutomatableWidget::getParameterName()
{
  return wrappedWidget->assignedParameter->getName();
}

MetaParameterManager* rsAutomatableWidget::getMetaParameterManager()
{
  MetaControlledParameter* mcp = 
    dynamic_cast<MetaControlledParameter*> (wrappedWidget->assignedParameter);
  if(mcp)
    return mcp->getMetaParameterManager();
  else
    return nullptr;
}

//=================================================================================================

rsAutomatableSlider::rsAutomatableSlider()
  : rsAutomatableWidget(this)
{

}

void rsAutomatableSlider::mouseDown(const MouseEvent& e)
{
  if(e.mods.isRightButtonDown())
  {
    openRightClickPopupMenu();
    //if(!isPopUpOpen())
    //  openRightClickPopupMenu();
    //else
    //  closePopUp();
  }
  else
  {
    RSlider::mouseDown(e);

    // doesn't work:
    //if(!isPopUpOpen())
    //  RSlider::mouseDown(e);
    //else
    //  closePopUp();
  }
}

void rsAutomatableSlider::rPopUpMenuChanged(RPopUpMenu* menuThatHasChanged)
{
  if( menuThatHasChanged != rightClickPopUp )
    return;

  RTreeViewNode *selectedItem = rightClickPopUp->getSelectedItem();
  if( selectedItem == NULL )
    return;
  int selectedIdentifier = selectedItem->getNodeIdentifier();

  switch( selectedIdentifier )
  {
  case ENTER_VALUE:   setValue(openModalNumberEntryField(getValue()),        true); break;
  case DEFAULT_VALUE: setValue(selectedItem->getNodeText().getDoubleValue(), true); break; //?
  default: rsAutomatableWidget::rPopUpMenuChanged(menuThatHasChanged);
  }
}

void rsAutomatableSlider::addPopUpMenuItems()
{
  addPopUpEnterValueItem();
  addPopUpDefaultValueItems();
  rsAutomatableWidget::addPopUpMenuItems();
}

void rsAutomatableSlider::addPopUpEnterValueItem()
{
  rightClickPopUp->addItem(ENTER_VALUE, "Enter Value");
}

void rsAutomatableSlider::addPopUpDefaultValueItems()
{
  if( defaultValues.size() > 0 )
  {
    RTreeViewNode *defaultValuesNode = new RTreeViewNode("Default Values");
    for(int i = 0; i < size(defaultValues); i++)
      defaultValuesNode->addChildNode(new RTreeViewNode(String(defaultValues[i]), DEFAULT_VALUE));
    rightClickPopUp->addTreeNodeItem(defaultValuesNode);
  }
}

//=================================================================================================

rsAutomatableComboBox::rsAutomatableComboBox()
  : rsAutomatableWidget(this)
{

}

void rsAutomatableComboBox::mouseDown(const MouseEvent& e)
{
  if( e.mods.isRightButtonDown() )
  {
    if(!isPopUpOpen())
      openRightClickPopupMenu();
    else
      rsAutomatableWidget::closePopUp();
  }
  else
    RComboBox::mouseDown(e);
}

void rsAutomatableComboBox::parameterChanged(Parameter* p)
{
  RWidget::parameterChanged(p);
  // not sure, why that's needed - isn't it supposed to be called anyway, i.e. if we don't override
  // parameterChanged, the RWidget baseclass method would be called? but somehow, it doesn't seem
  // to work
}

//=================================================================================================

rsAutomatableButton::rsAutomatableButton(const juce::String& buttonText)
  : RButton(buttonText), rsAutomatableWidget(this)
{

}

void rsAutomatableButton::mouseDown(const MouseEvent& e)
{
  if( e.mods.isRightButtonDown() )
  {
    if(!isPopUpOpen())
      openRightClickPopupMenu();
    else
      closePopUp();
  }
  else
    RButton::mouseDown(e);
}

void rsAutomatableButton::parameterChanged(Parameter* p)
{
  RWidget::parameterChanged(p);
}

//=================================================================================================

void rsModulationDepthSlider::rPopUpMenuChanged(RPopUpMenu* menu)
{
  int id = rightClickPopUp->getSelectedIdentifier();
  switch( id )
  {
  case MOD_DEPTH_MIN: setModDepthMin(openModalNumberEntryField(getModDepthMin())); break;
  case MOD_DEPTH_MAX: setModDepthMax(openModalNumberEntryField(getModDepthMax())); break;
  default: rsAutomatableSlider::rPopUpMenuChanged(menu);
  }
  juce::String text = menu->getSelectedText();
  typedef ModulationConnection::modModes MM;
  if(text == "Mode: Absolute")       setModMode(MM::ABSOLUTE);
  if(text == "Mode: Relative")       setModMode(MM::RELATIVE);
  if(text == "Mode: Exponential")    setModMode(MM::EXPONENTIAL);
  if(text == "Mode: Multiplicative") setModMode(MM::MULTIPLICATIVE);
}

void rsModulationDepthSlider::addPopUpMenuItems()
{
  rsAutomatableSlider::addPopUpMenuItems();
  addPopUpMinMaxAndModeItems();
}

void rsModulationDepthSlider::addPopUpMinMaxAndModeItems()
{
  rightClickPopUp->addItem(MOD_DEPTH_MIN, "Mod depth min");
  rightClickPopUp->addItem(MOD_DEPTH_MAX, "Mod depth max");
  typedef ModulationConnection::modModes MM;
  int m = getModMode();
  rightClickPopUp->addItem(MOD_MODE_ABSOLUTE,       "Mode: Absolute",       true, m == MM::ABSOLUTE);
  rightClickPopUp->addItem(MOD_MODE_RELATIVE,       "Mode: Relative",       true, m == MM::RELATIVE);
  rightClickPopUp->addItem(MOD_MODE_EXPONENTIAL,    "Mode: Exponential",    true, m == MM::EXPONENTIAL);
  rightClickPopUp->addItem(MOD_MODE_MULTIPLICATIVE, "Mode: Multiplicative", true, m == MM::MULTIPLICATIVE);
}

//=================================================================================================

rsModulationConnectionWidget::rsModulationConnectionWidget(ModulationConnection* connection,
  rsGarbageCollector* deletor)
  : rsDeletionRequester(deletor)
{
  addChildWidget(depthSlider  = new rsModulationDepthSlider(connection));
  addChildWidget(removeButton = new RClickButtonNotifyOnMouseUp(RButton::CLOSE));
}

void rsModulationConnectionWidget::resized()
{
  int w = getWidth();
  int h = getHeight();
  int buttonWidth  = h;
  int buttonMargin = 2;
  depthSlider->setBounds(0, 0, w-buttonWidth-buttonMargin, h);
  removeButton->setBounds(depthSlider->getRight()+buttonMargin, 0, buttonWidth, h);
}

//=================================================================================================

rsModulatableSlider::rsModulatableSlider() 
{

}

rsModulatableSlider::~rsModulatableSlider()
{
  ModulatableParameter* mp = dynamic_cast<ModulatableParameter*> (assignedParameter);
  if(mp)
    mp->deRegisterModulationTargetObserver(this);
}

void rsModulatableSlider::modulationsChanged()
{
  repaint();
}

void rsModulatableSlider::assignParameter(Parameter* p)
{
  rsAutomatableSlider::assignParameter(p);
  ModulatableParameter* mp = dynamic_cast<ModulatableParameter*> (p);
  if(mp)
    mp->registerModulationTargetObserver(this);
}

void rsModulatableSlider::paint(Graphics& g)
{
  rsAutomatableSlider::paint(g);
  if(hasModulation())
    g.fillAll(Colour::fromFloatRGBA(1.f, 0.f, 0.f, 0.125f)); // preliminary
}

bool rsModulatableSlider::hasModulation()
{
  ModulatableParameter* mp = dynamic_cast<ModulatableParameter*> (assignedParameter);
  if(mp && mp->hasConnectedSources())
    return true;
  return false;
}

//=================================================================================================

rsModulatableSliderAnimated::rsModulatableSliderAnimated() : rsModulatableSlider() 
{
  hasModConnections = hasModulation();

  // register at the repaint manager
}

rsModulatableSliderAnimated::~rsModulatableSliderAnimated()
{
  // de-register from the repaint manager
}

void rsModulatableSliderAnimated::modulationsChanged()
{
  hasModConnections = hasModulation(); // update our flag

  // hmmm...this functions seems to get called only when setting up modulations on the gui - not
  // when modulations are recalled from a preset load

  // or - even better than setting a flag would be to de/register with the RepaintManager!!

  // maybe we have trigger a repaint one last time, after all connections have been removed to get
  // rid of the indicator?
}

void rsModulatableSliderAnimated::paint(Graphics& g)
{
  rsModulatableSlider::paint(g);

  ModulatableParameter* mp = getModulatableParameter();
  if(mp != nullptr && hasModConnections)
  {
    double modValue = 0.5;  // preliminary
    modValue = mp->getModulatedValue();
    modValue = mp->valueToProportion(modValue); // unmap the 2nd (non-user defined) map

    float x = (float) (modValue * getWidth());
    float w = 4.f;  // mod-value indicator width
    g.setColour(Colours::red.withAlpha(0.5f));  // preliminary
    g.fillRect(x-w*0.5f, 0.f, w, (float) getHeight());
  }
}
