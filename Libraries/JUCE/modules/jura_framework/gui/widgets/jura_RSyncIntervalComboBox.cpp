
//-------------------------------------------------------------------------------------------------
// construction:

RSyncIntervalComboBox::RSyncIntervalComboBox(const String& componentName) : RComboBox(componentName)
{
  addItem(0, "1/16 beat");
  addItem(1, "1/8 beat" );
  addItem(2, "1/4 beat" );
  addItem(3, "1/2 beat" );
  addItem(4, "1 beat"   );
  addItem(5, "2 beats"  );
  addItem(6, "4 beats"  );
  addItem(7, "8 beats"  );
  addItem(8, "16 beats" );
  selectItemByIndex(4, false);
}

//-------------------------------------------------------------------------------------------------
// setup:

void RSyncIntervalComboBox::setValue(double newValue, bool sendMessage)
{
  if( isCloseTo(newValue, 1.0/16.0, 1.0/64.0) )
    selectItemByIndex(0, sendMessage);
  else if( isCloseTo(newValue, 1.0/8.0, 1.0/64.0) )
    selectItemByIndex(1, sendMessage);
  else if( isCloseTo(newValue, 1.0/4.0, 1.0/64.0) )
    selectItemByIndex(2, sendMessage);
  else if( isCloseTo(newValue, 1.0/2.0, 1.0/64.0) )
    selectItemByIndex(3, sendMessage);
  else if( isCloseTo(newValue, 1.0, 1.0/64.0) )
    selectItemByIndex(4, sendMessage);
  else if( isCloseTo(newValue, 2.0, 1.0/64.0) )
    selectItemByIndex(5, sendMessage);
  else if( isCloseTo(newValue, 4.0, 1.0/64.0) )
    selectItemByIndex(6, sendMessage);
  else if( isCloseTo(newValue, 8.0, 1.0/64.0) )
    selectItemByIndex(7, sendMessage);
  else if( isCloseTo(newValue, 16.0, 1.0/64.0) )
    selectItemByIndex(8, sendMessage);
  else
    selectItemByIndex(4, sendMessage);
}

//-------------------------------------------------------------------------------------------------
// inquiry:

double RSyncIntervalComboBox::getValueFromString(String stringToConvert)
{
  if( stringToConvert == String("1/16 beat") )
    return 1.0 / 16.0;
  else if( stringToConvert == String("1/8 beat") )
    return 1.0 / 8.0;
  else if( stringToConvert == String("1/4 beat") )
    return 1.0 / 4.0;
  else if( stringToConvert == String("1/2 beat") )
    return 1.0 / 2.0;
  else if( stringToConvert == String("1 beat") )
    return 1.0;
  else if( stringToConvert == String("2 beats") )
    return 2.0;
  else if( stringToConvert == String("4 beats") )
    return 4.0;
  else if( stringToConvert == String("8 beats") )
    return 8.0;
  else if( stringToConvert == String("16 beats") )
    return 16.0;
  else 
    return 1.0;
}

String RSyncIntervalComboBox::getStringFromValue(double valueToConvert)
{
  if( isCloseTo(valueToConvert, 1.0/16.0, 1.0/64.0) )
    return String("1/16 beat");
  else if( isCloseTo(valueToConvert, 1.0/8.0, 1.0/64.0) )
    return String("1/8 beat");
  else if( isCloseTo(valueToConvert, 1.0/4.0, 1.0/64.0) )
    return String("1/4 beat");
  else if( isCloseTo(valueToConvert, 1.0/2.0, 1.0/64.0) )
    return String("1/2 beat");
  else if( isCloseTo(valueToConvert, 1.0, 1.0/64.0) )
    return String("1 beat");
  else if( isCloseTo(valueToConvert, 2.0, 1.0/64.0) )
    return String("2 beats");
  else if( isCloseTo(valueToConvert, 4.0, 1.0/64.0) )
    return String("4 beats");
  else if( isCloseTo(valueToConvert, 8.0, 1.0/64.0) )
    return String("8 beats");
  else if( isCloseTo(valueToConvert, 16.0, 1.0/64.0) )
    return String("16 beats");
  else
    return String("1 beat");
}

double RSyncIntervalComboBox::getValue()
{
  switch( getSelectedItemIdentifier() )
  {
  case  0: return 1.0/16.0;
  case  1: return 1.0/8.0;
  case  2: return 1.0/4.0;
  case  3: return 1.0/2.0;
  case  4: return 1.0;
  case  5: return 2.0;
  case  6: return 4.0;
  case  7: return 8.0;
  case  8: return 16.0;

  default: return 1.0;
  }
}



