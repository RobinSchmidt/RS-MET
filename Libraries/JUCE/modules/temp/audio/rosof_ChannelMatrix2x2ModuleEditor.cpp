#include "rosof_ChannelMatrix2x2ModuleEditor.h"
using namespace rosof;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

ChannelMatrix2x2ModuleEditor::ChannelMatrix2x2ModuleEditor(CriticalSection *newPlugInLock, 
                                                           ChannelMatrix2x2AudioModule* newChannelMatrix2x2AudioModule) 
                                                           : AudioModuleEditor(newPlugInLock, newChannelMatrix2x2AudioModule)
{
  // set the plugIn-headline:
  setHeadlineText( juce::String(T("Channel Matrix 2x2")) );

  // assign the pointer to the rosic::ChannelMatrix2x2 object to be used as aduio engine:
  jassert(newChannelMatrix2x2AudioModule != NULL ); // you must pass a valid module here
  channelMatrix2x2AudioModule = newChannelMatrix2x2AudioModule;

  addWidget( labelLeftToLeft = new RTextField(juce::String(T("gLL="))) );
  labelLeftToLeft->setDescription(T("Amount by which left input goes to left output."));
  labelLeftToLeft->setDescriptionField(infoField);

  addWidget( editLabelLeftToLeft = new RTextEntryField(juce::String(T("1.0"))) );
  editLabelLeftToLeft->setDescription(labelLeftToLeft->getDescription());
  editLabelLeftToLeft->setDescriptionField(infoField);
  editLabelLeftToLeft->registerTextEntryFieldObserver(this);

  addWidget( labelRightToLeft = new RTextField(juce::String(T("gRL="))) );
  labelRightToLeft->setDescription(T("Amount by which right input goes to left output."));
  labelRightToLeft->setDescriptionField(infoField);

  addWidget( editLabelRightToLeft = new RTextEntryField(juce::String(T("1.0"))) );
  editLabelRightToLeft->setDescription(labelRightToLeft->getDescription());
  editLabelRightToLeft->setDescriptionField(infoField);
  editLabelRightToLeft->registerTextEntryFieldObserver(this);

  addWidget( labelLeftToRight = new RTextField(juce::String(T("gLR="))) );
  labelLeftToRight->setDescription(T("Amount by which left input goes to right output."));
  labelLeftToRight->setDescriptionField(infoField);

  addWidget( editLabelLeftToRight = new RTextEntryField(juce::String(T("1.0"))) );
  editLabelLeftToRight->setDescription(labelLeftToRight->getDescription());
  editLabelLeftToRight->setDescriptionField(infoField);
  editLabelLeftToRight->registerTextEntryFieldObserver(this);

  addWidget( labelRightToRight = new RTextField(juce::String(T("gRR="))) );
  labelRightToRight->setDescription(T("Amount by which right input goes to right output."));
  labelRightToRight->setDescriptionField(infoField);

  addWidget( editLabelRightToRight = new RTextEntryField(juce::String(T("1.0"))) );
  editLabelRightToRight->setDescription(labelRightToRight->getDescription());
  editLabelRightToRight->setDescriptionField(infoField);
  editLabelRightToRight->registerTextEntryFieldObserver(this);

  addWidget(leftEquationLabel = new RTextField(juce::String(T("yL = gLL*xL + gRL*xR"))) );
  leftEquationLabel->setJustification(Justification::centred);

  addWidget(rightEquationLabel = new RTextField(juce::String(T("yR = gLR*xL + gRR*xR"))) );
  rightEquationLabel->setJustification(Justification::centred);

  // set up the widgets:
  updateWidgetsAccordingToState();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void ChannelMatrix2x2ModuleEditor::textChanged(RTextEntryField *rTextEntryFieldThatHasChanged)
{
  if( channelMatrix2x2AudioModule == NULL )
    return;
  if( channelMatrix2x2AudioModule->wrappedChannelMatrix2x2 == NULL )
    return;

  if( rTextEntryFieldThatHasChanged == editLabelLeftToLeft )
  {
    channelMatrix2x2AudioModule->wrappedChannelMatrix2x2->setLeftToLeftGain(editLabelLeftToLeft->getText().getDoubleValue());
    editLabelLeftToLeft->setText(juce::String(channelMatrix2x2AudioModule->wrappedChannelMatrix2x2->getLeftToLeftGain()));
  }
  else if( rTextEntryFieldThatHasChanged == editLabelRightToLeft )
  {
    channelMatrix2x2AudioModule->wrappedChannelMatrix2x2->setRightToLeftGain(editLabelRightToLeft->getText().getDoubleValue());
    editLabelRightToLeft->setText(juce::String(channelMatrix2x2AudioModule->wrappedChannelMatrix2x2->getRightToLeftGain()));
  }
  else if( rTextEntryFieldThatHasChanged == editLabelLeftToRight )
  {
    channelMatrix2x2AudioModule->wrappedChannelMatrix2x2->setLeftToRightGain(editLabelLeftToRight->getText().getDoubleValue());
    editLabelLeftToRight->setText(juce::String(channelMatrix2x2AudioModule->wrappedChannelMatrix2x2->getLeftToRightGain()));
  }
  else if( rTextEntryFieldThatHasChanged == editLabelRightToRight )
  {
    channelMatrix2x2AudioModule->wrappedChannelMatrix2x2->setRightToRightGain(editLabelRightToRight->getText().getDoubleValue());
    editLabelRightToRight->setText(juce::String(channelMatrix2x2AudioModule->wrappedChannelMatrix2x2->getRightToRightGain()));
  }

  channelMatrix2x2AudioModule->markStateAsDirty();
}

void ChannelMatrix2x2ModuleEditor::updateWidgetsAccordingToState()
{
  if( channelMatrix2x2AudioModule == NULL )
    return;
  if( channelMatrix2x2AudioModule->wrappedChannelMatrix2x2 == NULL )
    return;

  // update the global widgets and automatable sliders:
  AudioModuleEditor::updateWidgetsAccordingToState();

  // update the non-automatable widgets:
  editLabelLeftToLeft->setText(  juce::String(channelMatrix2x2AudioModule->wrappedChannelMatrix2x2->getLeftToLeftGain()));
  editLabelRightToLeft->setText( juce::String(channelMatrix2x2AudioModule->wrappedChannelMatrix2x2->getRightToLeftGain()));
  editLabelLeftToRight->setText( juce::String(channelMatrix2x2AudioModule->wrappedChannelMatrix2x2->getLeftToRightGain()));
  editLabelRightToRight->setText(juce::String(channelMatrix2x2AudioModule->wrappedChannelMatrix2x2->getRightToRightGain()));
}

void ChannelMatrix2x2ModuleEditor::resized()
{
  linkPosition          = RIGHT_TO_HEADLINE;
  presetSectionPosition = BELOW_HEADLINE;
  AudioModuleEditor::resized();
  int x = 0;
  int y = getPresetSectionBottom()+8;
  int w = getWidth();
  int h = getHeight();

  labelLeftToLeft->setBounds(4, y, 48, 20);
  editLabelLeftToLeft->setBounds(labelLeftToLeft->getRight(), y, w/2-labelLeftToLeft->getWidth()-8, 20);
  labelRightToLeft->setBounds(w/2+4, y, 48, 20);
  editLabelRightToLeft->setBounds(labelRightToLeft->getRight(), y, w/2-labelRightToLeft->getWidth()-8, 20);

  y += 24;

  labelLeftToRight->setBounds(4, y, 48, 20);
  editLabelLeftToRight->setBounds(labelLeftToRight->getRight(), y, w/2-labelLeftToRight->getWidth()-8, 20);
  labelRightToRight->setBounds(w/2+4, y, 48, 20);
  editLabelRightToRight->setBounds(labelRightToRight->getRight(), y, w/2-labelRightToRight->getWidth()-8, 20);

  y += 28;

  leftEquationLabel->setBounds(4, y, w/2-8, 20);
  rightEquationLabel->setBounds(w/2+4, y, w/2-8, 20);
}

