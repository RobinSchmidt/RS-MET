#ifndef rojue_RHyperlinkButton_h
#define rojue_RHyperlinkButton_h

#include "rojue_RButton.h"

namespace rojue
{

  /** A button showing an underlined weblink, that will launch the link when it's clicked. */

  class RHyperlinkButton : public RButton
  {

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Creates a RHyperlinkButton with given text and refering to the given URL. */
    RHyperlinkButton(const juce::String& linkText, const URL& linkURL);

    /** Destructor. */
    ~RHyperlinkButton();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // setup:

    /** Changes the URL that the button will trigger. */
    void setURL(const URL& newURL) throw();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the URL that the button will trigger. */
    const URL& getURL() const throw() { return url; }

    //=====================================================================================================================================
    juce_UseDebuggingNewOperator

  protected:

    /** */
    void clicked();

    /**  */
    void paint(Graphics& g);

    URL url;


  private:

    RHyperlinkButton(const RHyperlinkButton&);
    const RHyperlinkButton& operator= (const RHyperlinkButton&);

  };

}

#endif  
