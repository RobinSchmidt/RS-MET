#include "rojue_RTimeGridComboBox.h"
using namespace rojue;

//-------------------------------------------------------------------------------------------------
// construction:

RTimeGridComboBox::RTimeGridComboBox(const String& componentName) : RComboBox(componentName)
{
  setDescription(String(T("Set time grid resolution (in seconds or beats).")));
  addItem(0, T("1/2")  );
  addItem(1, T("1/4")  );
  addItem(2, T("1/8")  );
  addItem(3, T("0.1")  );
  addItem(4, T("1/16") );
  addItem(5, T("1/32") );
  addItem(6, T("1/64") );
  addItem(7, T("0.01") );
  addItem(8, T("1/128"));
  selectItemByIndex(2, false);
}

//-------------------------------------------------------------------------------------------------
// setup:

/*
void RTimeGridComboBox::setValue(double newValue, bool sendMessage)
{
  if( isCloseTo(newValue, 1.0/16.0, 1.0/64.0) )
    setSelectedId(1, !sendMessage);
  else if( isCloseTo(newValue, 1.0/8.0, 1.0/64.0) )
    setSelectedId(2, !sendMessage);
  else if( isCloseTo(newValue, 1.0/4.0, 1.0/64.0) )
    setSelectedId(3, !sendMessage);
  else if( isCloseTo(newValue, 1.0/2.0, 1.0/64.0) )
    setSelectedId(4, !sendMessage);
  else if( isCloseTo(newValue, 1.0, 1.0/64.0) )
    setSelectedId(5, !sendMessage);
  else if( isCloseTo(newValue, 2.0, 1.0/64.0) )
    setSelectedId(6, !sendMessage);
  else if( isCloseTo(newValue, 4.0, 1.0/64.0) )
    setSelectedId(7, !sendMessage);
  else if( isCloseTo(newValue, 8.0, 1.0/64.0) )
    setSelectedId(8, !sendMessage);
  else if( isCloseTo(newValue, 16.0, 1.0/64.0) )
    setSelectedId(9, !sendMessage);
  else
    setSelectedId(5, !sendMessage);
}
*/

//-------------------------------------------------------------------------------------------------
// inquiry:
/*
double RTimeGridComboBox::getValueFromString(String stringToConvert)
{
  if( stringToConvert == String(T("1/16 beat")) )
    return 1.0 / 16.0;
  else if( stringToConvert == String(T("1/8 beat")) )
    return 1.0 / 8.0;
  else if( stringToConvert == String(T("1/4 beat")) )
    return 1.0 / 4.0;
  else if( stringToConvert == String(T("1/2 beat")) )
    return 1.0 / 2.0;
  else if( stringToConvert == String(T("1 beat")) )
    return 1.0;
  else if( stringToConvert == String(T("2 beats")) )
    return 2.0;
  else if( stringToConvert == String(T("4 beats")) )
    return 4.0;
  else if( stringToConvert == String(T("8 beats")) )
    return 8.0;
  else if( stringToConvert == String(T("16 beats")) )
    return 16.0;
  else 
    return 1.0;
}

String RTimeGridComboBox::getStringFromValue(double valueToConvert)
{
  if( isCloseTo(valueToConvert, 1.0/16.0, 1.0/64.0) )
    return String(T("1/16 beat"));
  else if( isCloseTo(valueToConvert, 1.0/8.0, 1.0/64.0) )
    return String(T("1/8 beat"));
  else if( isCloseTo(valueToConvert, 1.0/4.0, 1.0/64.0) )
    return String(T("1/4 beat"));
  else if( isCloseTo(valueToConvert, 1.0/2.0, 1.0/64.0) )
    return String(T("1/2 beat"));
  else if( isCloseTo(valueToConvert, 1.0, 1.0/64.0) )
    return String(T("1 beat"));
  else if( isCloseTo(valueToConvert, 2.0, 1.0/64.0) )
    return String(T("2 beats"));
  else if( isCloseTo(valueToConvert, 4.0, 1.0/64.0) )
    return String(T("4 beats"));
  else if( isCloseTo(valueToConvert, 8.0, 1.0/64.0) )
    return String(T("8 beats"));
  else if( isCloseTo(valueToConvert, 16.0, 1.0/64.0) )
    return String(T("16 beats"));
  else
    return String(T("1 beat"));
}
*/

double RTimeGridComboBox::getValue()
{
  switch( getSelectedItemIdentifier() )
  {
  case  0: return 1.0/2.0;
  case  1: return 1.0/4.0;
  case  2: return 1.0/8.0;
  case  3: return 0.1;
  case  4: return 1.0/16.0;
  case  5: return 1.0/32.0;
  case  6: return 1.0/64.0;
  case  7: return 0.01;
  case  8: return 1.0/128.0;

  default: return 1.0/4.0;
  }
}



