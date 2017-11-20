#include "rojue_FileSelectionBox.h"
using namespace rojue;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

FileSelectionBox::FileSelectionBox(const String& componentName, FileManager *fileManagerToUse)
{
  jassert( fileManagerToUse != NULL )
  fileManager = fileManagerToUse;
  fileManager->addFileManagerListener(this);

  addWidget( fileLabel = new RTextField( T("File:")) );
  fileLabel->setDescription(T("Currently active file"));
  fileLabel->setNoBackgroundAndOutline(true);

  addWidget( fileNameBox = new RTextField(String::empty) );
  fileNameBox->setDescription( fileLabel->getDescription() );
  fileNameBox->setNoBackgroundAndOutline(false);

  addWidget( loadButton = new RButton(String("Load")) );
  loadButton->addRButtonListener(this);
  loadButton->setDescription(String(T("Load new file")));
  loadButton->setClickingTogglesState(false);
  loadButton->setToggleState(false, false);

  addWidget( saveButton = new RButton(String("Save")) );
  saveButton->addRButtonListener(this);
  saveButton->setDescription(String(T("Save to file (save as)")));
  saveButton->setClickingTogglesState(false);
  saveButton->setToggleState(false, false);

  addWidget( minusButton = new RButton(RButton::MINUS) );
  minusButton->addRButtonListener(this);
  minusButton->setDescription(String(T("Skip to previous file in current directory")));
  minusButton->setClickingTogglesState(false);
  minusButton->setToggleState(false, false);

  addWidget( plusButton = new RButton(RButton::PLUS) );
  plusButton->addRButtonListener(this);
  plusButton->setDescription(String(T("Skip to next file in current directory")));
  plusButton->setClickingTogglesState(false);
  plusButton->setToggleState(false, false);

  labelPosition   = LABEL_ABOVE;
  buttonsPosition = BUTTONS_ABOVE;
  boxHeight       = 16;
}

FileSelectionBox::~FileSelectionBox()
{
  fileManager->removeFileManagerListener(this);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// callbacks:

void FileSelectionBox::rButtonClicked(RButton *buttonThatWasClicked)
{
  if( buttonThatWasClicked == loadButton )
    fileManager->openLoadingDialog();
  else if( buttonThatWasClicked == saveButton )
    fileManager->openSavingDialog();
  else if( buttonThatWasClicked == plusButton )
    fileManager->loadNextFile();
  else if( buttonThatWasClicked == minusButton )
    fileManager->loadPreviousFile();
  updateBoxContent();
}

void FileSelectionBox::activeFileChanged(FileManager *fileManagerThatHasChanged)
{
  updateBoxContent();
  sendChangeMessage();
}

void FileSelectionBox::resized()
{
  int boxX, boxY, boxWidth, buttonsY;
  int labelWidth = rojue::boldFont10px.getTextPixelWidth(fileLabel->getText(), 
    rojue::boldFont10px.getDefaultKerning()) + 4;
  int buttonsWidth = 2*boxHeight + 40 - 4;
  if( saveButton->isVisible() )
    buttonsWidth += (40-2);

  if( labelPosition = LABEL_ABOVE && buttonsPosition == BUTTONS_BELOW )
  {
    boxX     = 0;
    boxY     = boxHeight-2;
    boxWidth = getWidth();
    buttonsY = boxY+boxHeight-2;
  }
  else if( labelPosition = LABEL_ABOVE && buttonsPosition == BUTTONS_ABOVE)
  {
    boxX     = 0;
    boxY     = boxHeight-2;
    boxWidth = getWidth();
    buttonsY = 0;    
  }
  else if( labelPosition = LABEL_ABOVE && buttonsPosition == BUTTONS_RIGHT)
  {
    boxX     = 0;
    boxY     = boxHeight-2;
    boxWidth = getWidth()-buttonsWidth-2;
    buttonsY = boxY;  
  }
  else
  {
    jassertfalse;  // for the selected combination of positions is no implementation available
    boxX = boxY = boxWidth = buttonsY = 0; 
  }
  
  fileLabel->setBounds(0, 0, labelWidth, boxHeight);
  fileNameBox->setBounds(boxX, boxY, boxWidth, boxHeight);

  int x = getWidth()-boxHeight;
  plusButton->setBounds( x, buttonsY, boxHeight, boxHeight);
  x -= (boxHeight-2);
  minusButton->setBounds(x, buttonsY, boxHeight, boxHeight);
  x -= (40-2);
  loadButton->setBounds( x, buttonsY, 40,        boxHeight);
  x -= (40-2);
  saveButton->setBounds( x, buttonsY, 40,        boxHeight);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// others:

void FileSelectionBox::updateBoxContent()
{
  fileNameBox->setText(fileManager->getActiveFile().getFileName());
}
