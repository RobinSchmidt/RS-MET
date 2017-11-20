#include "rojue_RHyperlinkButton.h"
using namespace rojue;

RHyperlinkButton::RHyperlinkButton(const String& linkText, const URL& linkURL) : RButton(linkText), url(linkURL)
{
  setMouseCursor(MouseCursor::PointingHandCursor);
}

RHyperlinkButton::~RHyperlinkButton()
{

}

void RHyperlinkButton::setURL (const URL& newURL) throw()
{
  url = newURL;
}

void RHyperlinkButton::clicked()
{
  if(url.isWellFormed())
    url.launchInDefaultBrowser();
}

void RHyperlinkButton::paint(Graphics& g)
{
  drawBitmapFontText(g, 0, 0, getButtonText(), &boldFont10px, getTextColour());
}
