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
    descriptionField->setText(String()); // clear the old field
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
    descriptionField->setText(String());
}

void DescribedComponent::repaintOnMessageThread()
{
  if(MessageManager::getInstance()->isThisTheMessageThread())
    repaint();
  else
  {
      /*
       * Passing a safe ptr to avoid accidentally calling this while the editor is being
       * destroyed during an automation move.
       */
      Component::SafePointer<DescribedComponent> ptr = { this };
      MessageManager::callAsync([=] {if (ptr) ptr.getComponent()->repaint(); });
  }
    
}