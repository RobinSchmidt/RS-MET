//#include "rosic_TuningTable.h"
//using namespace rosic;

#include "../_third_party/MarkHenning/TuningMap.cpp"

//-------------------------------------------------------------------------------------------------
// construction/destruction:

TuningTable::TuningTable()
{
  name = NULL;
  resetToDefaults();
}

TuningTable::~TuningTable()
{
  if( name != NULL )
    delete[] name;
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void TuningTable::setName(const char *newName)
{
  // free old and allocate new memory for the name:
  if( name != NULL )
  {
    delete[] name;
    name = NULL;
  }
  if( newName != NULL )
  {
    int newLength = (int) strlen(newName);
    name          = new char[newLength+1];
    for(int c=0; c<=newLength; c++) // the <= is valid here, because we have one more cell allocated
      name[c] = newName[c];
  }
  //defaultState = false;
}

void TuningTable::resetToDefaults()
{
  setName("12-TET");
  masterTuneA4 = 440.0;
  detuneFactor = masterTuneA4 / 440.0;
  for(int i = 0; i < 128; i++)
    table[i] = RAPT::rsPitchToFreq((double) i);
  //defaultState = true;
}

void TuningTable::assignFrequency(int note, double newFrequency)
{
  if( note >= 0 && note <= 127 )
    table[note] = newFrequency;
  //defaultState = false;
}

void TuningTable::assignFrequencies(double *newFrequencies)
{
  for(int i = 0; i < 128; i++)
    table[i] = newFrequencies[i];
  //defaultState = false;
}

void TuningTable::setMasterTuneA4(double newTuneA4)
{
  masterTuneA4 = newTuneA4;
  detuneFactor = masterTuneA4 / 440.0;
  //defaultState = false;
}

bool TuningTable::loadFromTunFile(char *path)
{
  // create a temporary object of Mark Henning's class CTuningMap and let it read in the tuning
  // from the file:
  CTuningMap cTuningMap;
  bool success = cTuningMap.ReadFromFile(path);

  // if everything went right, the object should now contain the desired tuning - we retrieve now
  // the note frequencies to set up ourselves:
  if( success == true )
  {
    //for(int n = 0; n < 127; n++)  // shouldn't is be "n <= 127" or "n < 128"?
    for(int n = 0; n < 128; n++)  // new - needs tests
      assignFrequency(n, cTuningMap.GetNoteFreq(n));
    return true;
  }
  else
    return false;
  // ToDo: maybe apply a sanity check to the note-values - maybe restrict the frequencies to
  // the range 0..20000 or soemthing
}

//-------------------------------------------------------------------------------------------------
// note to frequency conversion:

char* TuningTable::getName()
{
  return name;
}

bool TuningTable::isInDefaultState()
{
  if( masterTuneA4 != 440.0 )
    return false;

  for(int i = 0; i < 128; i++)
  {
    if( table[i] != RAPT::rsPitchToFreq((double) i) )
      return false;
  }

  return true;
}

double TuningTable::getFrequency(int note)
{
  if( note >= 0 && note <= 127 )
    return detuneFactor * table[note];
  else
    return table[127];  // new since 2024/01/03
    //return 440.0;     // old
}

double TuningTable::getFrequency(double note)
{
  if( note >= 0.0 && note < 127.0 )
  {
    int    iPart  = (int) note;
    double fPart  = note - (double) iPart;
    double pitch1 = RAPT::rsFreqToPitch(table[iPart]);
    double pitch2 = RAPT::rsFreqToPitch(table[iPart+1]);
    double pitch  = (1.0-fPart)*pitch1 + fPart*pitch2;
    return detuneFactor * RAPT::rsPitchToFreq(pitch);
  }
  else
  {
    //return 440.0;  // old
    return detuneFactor * table[127];
    // New since 2024/03/01. We now saturate at the highest frequency in the table
  }
}

double TuningTable::getMasterTuneA4()
{
  return masterTuneA4;
}

