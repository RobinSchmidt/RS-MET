//-------------------------------------------------------------------------------------------------
// construction/destruction:

RMessageBox::RMessageBox()
{
  headlineTextField = new RTextField("Message Box");
  headlineTextField->setJustification(Justification::centred);
  headlineTextField->setNoBackgroundAndOutline(true);
  addAndMakeVisible(headlineTextField);

  bodyTextField = new RTextEditor("MessageBoxBodyTextEditor");
  //bodyTextField->setNoBackgroundAndOutline(true);
  bodyTextField->setMultiLine(true, true);
  bodyTextField->setReadOnly(true);
  bodyTextField->setPopupMenuEnabled(false);
  bodyTextField->setInterceptsMouseClicks(false, false);
  //Colour tb = Colours::transparentBlack;
  //bodyTextField->setColours(Colours::white, tb, tb, Colours::black, tb, tb);
  //bodyTextField->setColourScheme(colourScheme);
  bodyTextField->setText("This is a box in which some message can be displayed");
  addAndMakeVisible(bodyTextField);

  okButton = new RButton("OK");
  //okButton->setClickingTogglesState(false);
  okButton->addRButtonListener(this);
  addWidget(okButton);
}

RMessageBox::~RMessageBox()
{
  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// setup:

void RMessageBox::setHeadlineText(const String& newText)
{
  headlineTextField->setText(newText);
}

void RMessageBox::setBodyText(const String& newText)
{
  bodyTextField->setText(newText, false);
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void RMessageBox::rButtonClicked(RButton* button)
{
  setVisible(false);
}

void RMessageBox::resized()
{
  headlineTextField->setBounds(4, 8, getWidth()-8, 16);
  bodyTextField->setBounds(16, 32, getWidth()-32, getHeight()-64);
  okButton->setBounds(getWidth()/2-20, getHeight()-28, 40, 20);
}
