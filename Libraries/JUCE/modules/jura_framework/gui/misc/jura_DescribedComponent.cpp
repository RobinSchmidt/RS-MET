//#include "rojue_DescribedComponent.h"
//#include "rojue_RTextField.h"
//using namespace rojue;

// construction/destruction:

DescribedItem::DescribedItem(const String& newDescription) 
{
  description       = newDescription;
  descriptionField  = NULL;
}

DescribedItem::~DescribedItem()
{

}

// setup:

void DescribedItem::setDescription(const String &newDescription)
{
  description = newDescription;
  if( descriptionField != NULL )
    descriptionField->setText(description);
}

void DescribedItem::setDescriptionField(RTextField *newDescriptionField)
{
  if( descriptionField != NULL )
    descriptionField->setText(String::empty); // clear the old field
  descriptionField = newDescriptionField;
}

// inquiry:

String DescribedItem::getDescription() const
{
  return description;
}

RTextField* DescribedItem::getDescriptionField() const
{
  return descriptionField;
}

//=================================================================================================

void DescribedComponent::mouseEnter(const juce::MouseEvent &e)
{
  if( descriptionField != NULL )
    descriptionField->setText(description);
}

void DescribedComponent::mouseExit(const MouseEvent &e)
{
  if( descriptionField != NULL )
    descriptionField->setText(String::empty);
}

void repaintComponent(juce::Component* c)
{
  c->repaint();
}

void testFunction()
{

  int dummy = 0;
}

void DescribedComponent::repaintOnMessageThread()
{
  // preliminary - repaint only if this is the message thread - what we actually should do is
  // somehow trigger an asynchronous repaint on on the message thread, if this isn't the message 
  // thread
  MessageManager* mm = MessageManager::getInstance();
  if(mm->isThisTheMessageThread())
    repaint();
  else
  {
    // https://juce.com/doc/classMessageManager
    // use: void * 	callFunctionOnMessageThread (MessageCallbackFunction *callback, void *userData)
    // ...but no - that blocks
    // or void MessageManager::callAsync	(	FunctionType 	functionToCall	)	

    // or this:
    // https://juce.com/doc/classCallbackMessage

    //MessageManager::callAsync(repaintComponent(this));

    //auto func = testFunction;

    //auto func = [this]{ repaint(); };
    //MessageManager::callAsync(func);

    MessageManager::callAsync([this]{ repaint(); });

    //jassertfalse;
  }
}