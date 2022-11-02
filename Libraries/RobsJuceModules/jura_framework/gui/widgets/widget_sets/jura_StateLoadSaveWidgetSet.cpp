
//-------------------------------------------------------------------------------------------------
// construction/destruction:

StateLoadSaveWidgetSet::StateLoadSaveWidgetSet(const String& newStateLoadSaveWidgetSetName) 
{
  addWidget( stateLabel = new RTextField( String("State:")) );
  stateLabel->setDescription(String("Name of current file (if any)"));
  stateLabel->setNoBackgroundAndOutline(true);

  addWidget(stateFileNameLabel = new RTextField(String()) );
  stateFileNameLabel->setDescription(String("Name of current file (if any)"));
  stateFileNameLabel->setNoBackgroundAndOutline(false);

  addWidget( stateLoadButton = new RClickButton(String("Load")) );
  stateLoadButton->addRButtonListener(this);
  stateLoadButton->setDescription(String("Load setting from file"));
  stateLoadButton->setClickingTogglesState(false);
  stateLoadButton->setToggleState(false, false);

  addWidget( stateSaveButton = new RClickButton(String("Save")) );
  stateSaveButton->addRButtonListener(this);
  stateSaveButton->setDescription(String("Save current setting to file"));
  stateSaveButton->setClickingTogglesState(false);
  stateSaveButton->setToggleState(false, false);

  //addWidget( stateMinusButton = new RButton(RButton::MINUS) );
  addWidget( stateMinusButton = new RClickButton(RButton::ARROW_LEFT) );
  stateMinusButton->addRButtonListener(this);
  stateMinusButton->setDescription(String("Skip to previous file in current directory"));
  stateMinusButton->setClickingTogglesState(false);
  stateMinusButton->setToggleState(false, false);

  //addWidget( statePlusButton = new RButton(RButton::PLUS) );
  addWidget( statePlusButton = new RClickButton(RButton::ARROW_RIGHT) );
  statePlusButton->addRButtonListener(this);
  statePlusButton->setDescription(String("Skip to next file in current directory"));
  statePlusButton->setClickingTogglesState(false);
  statePlusButton->setToggleState(false, false);

  layout = ONE_LINE;
}

StateLoadSaveWidgetSet::~StateLoadSaveWidgetSet()
{
  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// setup:

void StateLoadSaveWidgetSet::updateStateNameField()
{
  if( watchedStateManager == NULL )
    return;
  else
    stateFileNameLabel->setText(watchedStateManager->getStateNameWithStarIfDirty());
}

void StateLoadSaveWidgetSet::setDescriptionField(RTextField* newDescriptionField)
{
  stateLabel->setDescriptionField(newDescriptionField);
  stateFileNameLabel->setDescriptionField(newDescriptionField);
  stateSaveButton->setDescriptionField(newDescriptionField);
  stateLoadButton->setDescriptionField(newDescriptionField);
  statePlusButton->setDescriptionField(newDescriptionField);
  stateMinusButton->setDescriptionField(newDescriptionField);
}

void StateLoadSaveWidgetSet::setLayout(int newLayout)
{
  layout = newLayout;
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void StateLoadSaveWidgetSet::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( watchedStateManager == NULL )
    return;

  StateFileManager* underlyingStateFileManager = 
    dynamic_cast<StateFileManager*> (watchedStateManager);
  if( underlyingStateFileManager == NULL )
    return;

  if( buttonThatWasClicked == stateLoadButton )
  {
    underlyingStateFileManager->openLoadingDialog();
    sendChangeMessage();
  }
  else if( buttonThatWasClicked == stateSaveButton )
    underlyingStateFileManager->openSavingDialog();
  else if( buttonThatWasClicked == statePlusButton )
  {
    underlyingStateFileManager->loadNextFile();
    sendChangeMessage();
  }
  else if( buttonThatWasClicked == stateMinusButton )
  {
    underlyingStateFileManager->loadPreviousFile();
    sendChangeMessage();
  }
}

void StateLoadSaveWidgetSet::stateDirtyFlagChanged(StateManager *stateManager)
{
  updateStateNameField();
}

/*
void StateLoadSaveWidgetSet::paint(Graphics &g)
{
  g.fillAll(Colours::red); // only for test
}


void StateLoadSaveWidgetSet::paintOverChildren(Graphics &g)
{

}
*/

void StateLoadSaveWidgetSet::resized()
{
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();
  int x1, x2;

  if( layout == LABEL_AND_BUTTONS_ABOVE ) 
  {
    int h2 = h/2;
    stateLabel->setBounds(      0, y,  w, h2);
    x  = w-h2;
    statePlusButton->setBounds (x, y, h2, h2);
    x -= (h2-2);
    stateMinusButton->setBounds(x, y, h2, h2);
    x -= (40-2);
    stateLoadButton->setBounds (x, y, 40, h2);
    x -= (40-2);
    stateSaveButton->setBounds (x, y, 40, h2);
    x  = 0;
    y += (h2-2);
    stateFileNameLabel->setBounds(0, y, w, h2);
  }
  else // ONE_LINE
  {
    x2  = getWidth()-4;
    w   = 16;
    x1  = x2-w;
    statePlusButton->setBounds(x1, 0, w, 16);
    x1 -= w; 
    stateMinusButton->setBounds(x1+2, 0, w, 16);
    w   = 40;
    x1  = stateMinusButton->getX()-(w-2);
    stateLoadButton->setBounds(x1, 0, w, 16);
    x1 -= (w-2);
    stateSaveButton->setBounds(x1, 0, w, 16);
    x2 = stateSaveButton->getX()+2;
    x1 = 4;
    w  = x2-x1;
    stateFileNameLabel->setBounds(x1, 0, w, 16);
  }
}
