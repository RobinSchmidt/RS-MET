

#ifdef RS_BLAH  // for commenting out

//#include "rojue_WidgetSet.h"
//using namespace rojue;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

WidgetSetOld::WidgetSetOld() 
{

}

WidgetSetOld::~WidgetSetOld()
{

}

//-------------------------------------------------------------------------------------------------
// setup:

void WidgetSetOld::setColourScheme(const WidgetColourScheme& newColourScheme)
{
  widgets.getLock().enter();
  for(int w=0; w<widgets.size(); w++)
    widgets[w]->setColourScheme(newColourScheme);
  widgets.getLock().exit();
}

void WidgetSetOld::setColourSchemeFromXml(const XmlElement* widgetColours)
{
  if( widgetColours != NULL )
  {
    WidgetColourScheme tmpColourScheme;
    //tmpColourScheme.setColourSchemeFromXml(*widgetColours);
    setColourScheme(tmpColourScheme);
  }
}

void WidgetSetOld::setDescriptionField(RTextField *newDescriptionField)
{
  widgets.getLock().enter();
  for(int w=0; w<widgets.size(); w++)
  {
    widgets[w]->setDescriptionField(newDescriptionField);
  }
  widgets.getLock().exit();
}

//-------------------------------------------------------------------------------------------------
// others:

void WidgetSetOld::addWidget(RWidget *widgetToAdd, bool addAsChildComponent, bool makeVisible)
{
  widgets.getLock().enter();
  widgets.addIfNotAlreadyThere(widgetToAdd);
  widgets.getLock().exit();
  if( addAsChildComponent == true )
  {
    Component* comp = dynamic_cast<Component*> (widgetToAdd);
    if( comp != NULL )
    {
      addChildComponent(comp);
      if( makeVisible == true )
        comp->setVisible(true);
    }
  }

  widgetToAdd->setColourScheme(widgetColourScheme);

}
// -Compare to baseclass method 
//  -it has also a call to widgetToAdd->setColourScheme(widgetColourScheme);
//  -why do we override it at all?
// 

/*
void WidgetSet::addWidget(RLabel *labelToAdd, bool addAsChildComponent, bool makeVisible)
{
  labels.getLock().enter();
  labels.addIfNotAlreadyThere(labelToAdd);
  labels.getLock().exit();
  if( addAsChildComponent == true )
  {
    addChildComponent(labelToAdd);
    if( makeVisible == true )
      labelToAdd->setVisible(true);
  }
}
*/
#endif // end of commenting out