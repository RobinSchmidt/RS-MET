#include "rosic_TuningTable.h"
using namespace rosic;

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
  for(int i=0; i<128; i++)
    table[i] = pitchToFreq((double) i);
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
  for(int i=0; i<128; i++)
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
    for(int n=0; n<127; n++)
      assignFrequency(n, cTuningMap.GetNoteFreq(n));
    return true;
  }
  else
    return false;
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

  for(int i=0; i<128; i++)
  {
    if( table[i] != pitchToFreq((double) i) )
      return false;
  }

  return true;
}

double TuningTable::getFrequency(int note)
{
  if( note >= 0 && note <= 127 )
    return detuneFactor * table[note];
  else
    return 440.0;
}

double TuningTable::getFrequency(double note)
{
  if( note >= 0.0 && note < 127.0 )
  {
    int    iPart  = (int) note;
    double fPart  = note - (double) iPart;
    double pitch1 = freqToPitch(table[iPart]);
    double pitch2 = freqToPitch(table[iPart+1]);
    double pitch  = (1.0-fPart)*pitch1 + fPart*pitch2;
    return detuneFactor * pitchToFreq(pitch);
  }
  else
    return 440.0;
}

double TuningTable::getMasterTuneA4()
{
  return masterTuneA4;
}

