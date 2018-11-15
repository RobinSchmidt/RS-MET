#ifndef jura_ImageSavingDialog_h
#define jura_ImageSavingDialog_h

/** This class is to be used as content component for a DialogWindow object for saving images - it 
is mainly intended to export the contents of a Component to image files. */

class ImageSavingDialog : public Component, public Button::Listener, public ComboBox::Listener,
	public Label::Listener
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  ImageSavingDialog(rsPlot *owner = NULL, int defaultWidth = 320,
    int defaultHeight = 320, const String &defaultFormat = String("png"),
    const File& defaultTargetFile = File::nonexistent);

  /** Destructor. */
  virtual ~ImageSavingDialog();

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the desired pixel-width which as chosen by the user in the dialog. */
  virtual int getSelectedPixelWidth();

  /** Returns the desired pixel-height which as chosen by the user in the dialog. */
  virtual int getSelectedPixelHeight();

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  /** Implements the purely virtual rButtonClicked()-method of the ButtonListener base-class. */
  virtual void buttonClicked(juce::Button *buttonThatWasClicked);

  /** Implements the purely virtual comboBoxChanged()-method of the ComboBoxListener
  base-class.*/
  virtual void comboBoxChanged(ComboBox *comboBoxThatHasChanged);

  /** Implements the purely virtual labelTextChanged()-method of the LablelListener
  base-class. */
  virtual void labelTextChanged(Label *labelThatHasChanged);


protected:

  void saveNow();

  rsPlot *ownerSystem;

  Label        *targetFileLabel;
  Label        *targetFileEditLabel;
  TextButton   *browseButton;
  Label        *pixelWidthLabel;
  Label        *pixelWidthEditLabel;
  Label        *pixelHeightLabel;
  Label        *pixelHeightEditLabel;
  Label        *formatLabel;
  ComboBox     *formatComboBox;
  TextButton   *saveButton;

  juce_UseDebuggingNewOperator;
};

#endif
