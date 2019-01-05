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

void DescribedComponent::repaintOnMessageThread()
{
	triggerAsyncUpdate();
}