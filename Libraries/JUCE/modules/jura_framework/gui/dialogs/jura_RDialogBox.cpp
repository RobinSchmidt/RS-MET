RDialogBox::RDialogBox()
{
  setHeadlineText(String("Dialog"));
  setHeadlineStyle(Editor::MAIN_HEADLINE);

  addWidget( okButton = new RClickButtonNotifyOnMouseUp(String("OK")) );
  okButton->setClickingTogglesState(false);
  okButton->addRButtonListener(this);

  addWidget( cancelButton = new RClickButtonNotifyOnMouseUp(String("Cancel")) );
  cancelButton->setClickingTogglesState(false);
  cancelButton->addRButtonListener(this);

  setSize(200, 100);
}

RDialogBox::~RDialogBox()
{
  deleteAllChildren();
}
 
//-------------------------------------------------------------------------------------------------
// setup:

void RDialogBox::addListener(RDialogBoxListener* const listenerToAdd)
{
  listeners.getLock().enter();
  listeners.addIfNotAlreadyThere(listenerToAdd);
  listeners.getLock().exit();
}

void RDialogBox::removeListener(RDialogBoxListener* const listenerToRemove)
{
  listeners.getLock().enter();
  listeners.removeFirstMatchingValue(listenerToRemove);
  listeners.getLock().exit();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void RDialogBox::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( buttonThatWasClicked == cancelButton )
    sendCancelClickedNotification();
  else if( buttonThatWasClicked == okButton )
    sendOKClickedNotification();
}

void RDialogBox::resized()
{
  int y = getHeight() - 28;
  int x = getWidth() / 2;
  okButton->setBounds(    x-60-4, y, 60, 20);
  cancelButton->setBounds(x+4,    y, 60, 20);
}

//-------------------------------------------------------------------------------------------------
// others:

void RDialogBox::sendChangeNotification()
{
  listeners.getLock().enter();
  for(int i=0; i<listeners.size(); i++)
    listeners[i]->rDialogBoxChanged(this);
  listeners.getLock().exit();
}

void RDialogBox::sendOKClickedNotification()
{
  listeners.getLock().enter();
  for(int i=0; i<listeners.size(); i++)
    listeners[i]->rDialogBoxOKClicked(this);
  listeners.getLock().exit();
}

void RDialogBox::sendCancelClickedNotification()
{
  listeners.getLock().enter();
  for(int i=0; i<listeners.size(); i++)
    listeners[i]->rDialogBoxCancelClicked(this);
  listeners.getLock().exit();
}