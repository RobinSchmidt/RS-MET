#include "jura_StringTools.h"

char* toZeroTerminatedString(String stringToConvert)
{
  int length = stringToConvert.length();
  if( length <= 0 )
    return NULL;

  char* stringC = new char[length+1];
  //stringToConvert.copyToBuffer(stringC, length);
  //stringToConvert.copyToCString(stringC, length);
  size_t numWritten = stringToConvert.copyToUTF8(stringC, length+1);
  jassert((int) numWritten == length+1);
  return stringC;
}

String extractFileName(String fullPath)
{
  String fileName = fullPath;
  if( fileName.contains("\\") )
    fileName = fileName.fromLastOccurrenceOf("\\", false, false);
  if( fileName.contains("/") )
    fileName = fileName.fromLastOccurrenceOf("/", false, false);
  return fileName;
}

String createAudioFileInfoString(File fileToCreateStringFrom)
{
  // create a file input stream as needed for the reader
  FileInputStream *fileInputStream = new FileInputStream(fileToCreateStringFrom);
  if( fileInputStream == NULL )
    return String("Error opening file in 'createAudioFileInfoString'");

  bool   fileIsReadable = false;
  int64  numSamples  = 0;
  int    numChannels = 0;
  double sampleRate  = 0.0;
  int    bitDepth    = 0;
  String formatString;

  if( fileToCreateStringFrom.hasFileExtension(String("wav")) )
  {
    WavAudioFormat wavAudioFormat;
    fileIsReadable = wavAudioFormat.canHandleFile(fileToCreateStringFrom);
    formatString = String("wav");
    if( fileIsReadable )
    {
      AudioFormatReader* wavReader = wavAudioFormat.createReaderFor(fileInputStream, true);
      if( wavReader != NULL )
      {
        numSamples  = wavReader->lengthInSamples;
        numChannels = wavReader->numChannels;
        sampleRate  = wavReader->sampleRate;
        bitDepth    = wavReader->bitsPerSample;
        delete wavReader;
      }
      else
      {
        delete fileInputStream;
        return String("Error reading file in 'createAudioFileInfoString'");
      }
    }
    else
      return String("Error reading file in 'createAudioFileInfoString'");
  }
  else if( fileToCreateStringFrom.hasFileExtension(String("flac")) )
  {
    FlacAudioFormat flacAudioFormat;
    fileIsReadable = flacAudioFormat.canHandleFile(fileToCreateStringFrom);
    formatString = String("flac");
    if( fileIsReadable )
    {
      AudioFormatReader* flacReader = flacAudioFormat.createReaderFor(fileInputStream, true);
      if( flacReader != NULL )
      {
        numSamples  = flacReader->lengthInSamples;
        numChannels = flacReader->numChannels;
        sampleRate  = flacReader->sampleRate;
        bitDepth    = flacReader->bitsPerSample;
        delete flacReader;
      }
      else
      {
        delete fileInputStream;
        return String("Error reading file in 'createAudioFileInfoString'");
      }
    }
    else
      return String("Error reading file in 'createAudioFileInfoString'");
  }
  else
  {
    return String("unknown file format");
  }

  String channelString;
  if( numChannels == 1 )
    channelString = String("mono");
  else if( numChannels == 2 )
    channelString = String("stereo");
  else
    channelString = String(numChannels) + String(" channels");

  String infoString = formatString + String(", ") + String(sampleRate) + String(" Hz, ")
    + String(bitDepth) + String(" Bit, ") + String(numSamples) + String(" frames, ")
    + channelString;

  return infoString;
}

double toDouble(const juce::String& s)
{
  String st = s.trim();
  if(st == "-inf")
    return -std::numeric_limits<double>::infinity();
  return s.getDoubleValue();
}

inline double rsAmp2dB(double amp)
{
  // temporary - todo: remove and use function from RAPT instead
  return 8.6858896380650365530225783783321 * log(amp);
}
String amplitudeRawAndInDecibels(double amplitudeRaw)
{
  amplitudeRaw = fabs(amplitudeRaw);
  if( amplitudeRaw > 0.0 )
    return valueToStringTotal5(amplitudeRaw) + String(", ")
    + decibelsToStringWithUnit2(rsAmp2dB(amplitudeRaw));
  else
    return valueToStringTotal5(amplitudeRaw) + String(", -Inf dB");
}

String beatsToStringWithUnit4(double value)
{
  return String(value, 4) + String(" bts");
}

String centsToStringWithUnit2(double value)
{
  return String(value, 2) + String(" ct");
}

String decibelsToStringWithUnit(double value)
{
  return String(value) + String(" dB");
}

String decibelsToStringWithUnit1(double value)
{
  return String(value, 1) + String(" dB");
}

String decibelsToStringWithUnit2(double value)
{
  return String(value, 2) + String(" dB");
}

String degreesToStringWithUnit0(double value)
{
  return String(value, 2) + String(" deg");
  //return String(value, 0) + String(CharPointer_UTF8("°")); // ° not available in my pixelfont
  //return String(value, 0) + String("°"); // old - raises assertion
}

String frequencyInHzAndAsNote(double frequencyInHz)
{
  return hertzToStringWithUnitTotal5(frequencyInHz) + String(", ") + frequencyToNoteString(frequencyInHz);
}


inline double rsLog2(double x)
{
  return log(x) / log(2.0); // remove..
}
inline double rsFreqToPitch(double freq)
{
  // remove and use RAPT function instead...
  return 12.0 * rsLog2(freq / 440) + 69.0;
}
String frequencyToNoteString(double frequencyInHz)
{
  double pitch = rsFreqToPitch(frequencyInHz);
  int    note  = (int) round(pitch);
  double cents = 100.0*(pitch-note);
  return midiNoteToString(note) + String(" ") + valueToStringWithSign0(cents) + String(" cents");
  //return midiNoteToString(note) + String(T(", ")) + valueToString0(cents) + String(T(" cents"));
}

String decibelsPerOctaveToString(double value)
{
  return String( (int) round(value) ) + String(" dB/oct");
}

String decibelsPerOctaveToString2(double value)
{
  return String(value, 2) + String(" dB/oct");
}

String getSignAsString(double value)
{
  if( value > 0.0 )
    return String("+");
  else if( value < 0.0 )
    return String("-");
  else
    return String::empty;
}

String hertzToStringWithUnit1(double value)
{
  return String(value, 1) + String(" Hz");
}

String hertzToStringWithUnit2(double value)
{
  return String(value, 2) + String(" Hz");
}

String intToStringWithLeadingZeros(int value, int minNumDigits)
{
  String result = String(value);
  int big = (int) pow(10.0, minNumDigits-1.0);
  while( big > value )
  {
    big   /= 10;
    result = String("0") + result;
  }
  return result;
}

String hertzToStringWithUnitTotal5(double value)
{
  if( value >= 10000.0 )
    return String( (int) round(value) ) + String(" Hz");
  else if( value >= 1000.0 )
    return String(value, 1) + String(" Hz");
  else if( value >= 100.0 )
    return String(value, 2) + String(" Hz");
  else if( value >= 10.0 )
    return String(value, 3) + String(" Hz");
  else if( value >= 1.0 )
    return String(value, 4) + String(" Hz");
  else
    return String(value, 5) + String(" Hz");
}

String kiloToString0(double value)
{
  return String( (int) round(0.001*value) ) + String("k");
}

String midiMessageToString(MidiMessage message, bool addNewLine)
{
  String messageString = String::empty;

  int hours, minutes, seconds, frames;  // needed for the MMC goto message and time-code messages

  if( message.isNoteOn() )
  {
    messageString += String("Note On, ") + String("Channel: ") + String(message.getChannel())
      + String(", Key: ") + String(message.getNoteNumber())
      + String(" (") + message.getMidiNoteName(message.getNoteNumber(), true, true, 4)
      + String("), Velocity: ") + String(message.getVelocity());
  }
  else if( message.isNoteOff() )
  {
    messageString += String("Note Off, ") + String("Channel: ") + String(message.getChannel())
      + String(", Key: ") + String(message.getNoteNumber())
      + String(" (") + message.getMidiNoteName(message.getNoteNumber(), true, true, 4)
      + String("), Velocity: ") + String(message.getVelocity());
  }
  else if( message.isController() )
  {
    messageString += String("Controller, ") + String("Channel: ") + String(message.getChannel())
      + String(", Number: ") + String(message.getControllerNumber())
      + String(" (") + MidiMessage::getControllerName(message.getControllerNumber())
      + String("), Value: ") + String(message.getControllerValue());
  }
  else if( message.isPitchWheel() )
  {
    messageString += String("Pitch Wheel, ") + String("Channel: ") + String(message.getChannel())
      + String(", Value: ") + String(message.getPitchWheelValue());
  }
  else if( message.isProgramChange() )
  {
    messageString += String("Program Change, ") + String("Channel: ") + String(message.getChannel())
      + String(", Number: ") + String(message.getProgramChangeNumber());
  }
  else if( message.isAftertouch() )
  {
    messageString += String("Aftertouch, ") + String("Channel: ") + String(message.getChannel())
      + String(", Key: ") + String(message.getNoteNumber())
      + String(", Value: ") + String(message.getAfterTouchValue());
  }
  else if( message.isChannelPressure() )
  {
    messageString += String("Channel Pressure, ") + String("Channel: ") + String(message.getChannel())
      + String(", Value: ") + String(message.getChannelPressureValue());
  }
  else if( message.isSysEx() )
  {
    messageString += String("SysEx");
  }
  else if( message.isMetaEvent() )
  {
    messageString += String("Meta Event: ");
    if( message.isTrackMetaEvent() )
      messageString += String("track");
    else if( message.isEndOfTrackMetaEvent() )
      messageString += String("end of track");
    else if( message.isTrackNameEvent() )
      messageString += String("track name: ") + message.getTextFromTextMetaEvent();
    else if( message.isTextMetaEvent() )
      messageString += String("text: ") + message.getTextFromTextMetaEvent();
    else if( message.isTempoMetaEvent() )
    {
      double bpm = 60.0 / message.getTempoSecondsPerQuarterNote();
      messageString += String("tempo: ") + String(bpm,3) + String(" bpm");
    }
    else if( message.isTimeSignatureMetaEvent() )
    {
      int num, den;
      message.getTimeSignatureInfo(num, den);
      messageString += String("time signature: ") + String(num) + String("/") + String(den);
    }
    else if( message.isKeySignatureMetaEvent() )
    {
      messageString += String("key signature: ")
        + String(message.getKeySignatureNumberOfSharpsOrFlats()) + String(" sharps");
    }
    else if( message.isMidiChannelMetaEvent() )
    {
      messageString += String("midi channel: ")
        + String(message.getMidiChannelMetaEventChannel());
    }
  }
  else if( message.isMidiStart() )
    messageString += String("Start");
  else if( message.isMidiStop() )
    messageString += String("Stop");
  else if( message.isMidiContinue() )
    messageString += String("Continue");
  else if( message.isMidiClock() )
    messageString += String("Clock");
  else if( message.isSongPositionPointer() )
  {
    messageString += String("Song position: ") + String(message.getSongPositionPointerMidiBeat()/4)
      + String(" quarter notes");
  }
  else if( message.isQuarterFrame() )
  {
    messageString += String("MIDI Time Code (MTC), quarter frame - sequence number: ")
      + String(message.getQuarterFrameSequenceNumber()) + String(", value: ")
      + String(message.getQuarterFrameValue() );
  }
  else if( message.isFullFrame() )
  {
    MidiMessage::SmpteTimecodeType timecodeType;
    message.getFullFrameParameters(hours, minutes, seconds, frames, timecodeType);

    String timeCodeString = String(hours) + String(":") + String(minutes) + String(":")
      + String(seconds) + String(":") + String(frames);

    String typeString;
    switch( timecodeType )
    {
    case MidiMessage::fps24:     typeString = String("24 fps");                break;
    case MidiMessage::fps25:     typeString = String("25 fps");                break;
    case MidiMessage::fps30:     typeString = String("30 fps");                break;
    case MidiMessage::fps30drop: typeString = String("30 fps with dropframe"); break;
    default: typeString = String("unknown");
    }

    messageString += String("MIDI Time Code (MTC), full frame: ") + timeCodeString
      + String(" at framerate: ") + typeString;
  }
  else if( message.isMidiMachineControlMessage ()  )
  {
    messageString += String("MIDI Machine Control (MMC): ");
    MidiMessage::MidiMachineControlCommand  mmc = message.getMidiMachineControlCommand();
    if( mmc == MidiMessage::mmc_deferredplay )
      messageString += String("deferred play");
    else if( mmc == MidiMessage::mmc_play )
      messageString += String("play");
    else if( mmc == MidiMessage::mmc_stop )
      messageString += String("stop");
    else if( mmc == MidiMessage::mmc_fastforward )
      messageString += String("fast forward");
    else if( mmc == MidiMessage::mmc_rewind )
      messageString += String("rewind");
    else if( mmc == MidiMessage::mmc_recordStart )
      messageString += String("record start");
    else if( mmc == MidiMessage::mmc_recordStop )
      messageString += String("record stop");
    else if( mmc == MidiMessage::mmc_pause )
      messageString += String("pause");
  }
  else if( message.isMidiMachineControlGoto(hours, minutes, seconds, frames) )
  {
    messageString += String("MIDI Machine Control (MMC) goto: ");
    messageString += String(hours) + String(":") + String(minutes) + String(":") +
      String(seconds) + String(":") + String(frames);
  }
  else
  {
    messageString = String("unknown message");
  }

  if( addNewLine == true && messageString != String::empty )
    return messageString + String("\n");
  else
    return messageString;
}

/*
String midiControllerToString(int midiControllerNumber)
{
  switch( midiControllerNumber )
  {
  case   0: return String(T("Bank Select, coarse"));
  case   1: return String(T("Modulation Wheel, coarse"));
  case   2: return String(T("Breath Controller, coarse"));
  case   3: return String(T("unknown"));
  case   4: return String(T("Foot Pedal, coarse"));
  case   5: return String(T("Portamento Time, coarse"));
  case   6: return String(T("Data Entry, coarse"));
  case   7: return String(T("Volume, coarse"));
  case   8: return String(T("Balance, coarse"));
  case   9: return String(T("unknown"));
  case  10: return String(T("Pan Position, coarse"));
  case  11: return String(T("Expression, coarse"));

  case  12: return String(T("Effect Control 1, coarse"));
  case  13: return String(T("Effect Control 2, coarse"));
  case  14: return String(T("unknown"));
  case  15: return String(T("unknown"));
  case  16: return String(T("General Purpose Slider 1"));
  case  17: return String(T("General Purpose Slider 2"));
  case  18: return String(T("General Purpose Slider 3"));
  case  19: return String(T("General Purpose Slider 4"));
  case  20: return String(T("unknown"));
  case  21: return String(T("unknown"));
  case  22: return String(T("unknown"));
  case  23: return String(T("unknown"));

  case  24: return String(T("unknown"));
  case  25: return String(T("unknown"));
  case  26: return String(T("unknown"));
  case  27: return String(T("unknown"));
  case  28: return String(T("unknown"));
  case  29: return String(T("unknown"));
  case  30: return String(T("unknown"));
  case  31: return String(T("unknown"));
  case  32: return String(T("Bank Select, fine"));
  case  33: return String(T("Modulation Wheel, fine"));
  case  34: return String(T("Breath Controller, fine"));
  case  35: return String(T("unknown"));

  case  36: return String(T("Foot Pedal, fine"));
  case  37: return String(T("Portamento Time, fine"));
  case  38: return String(T("Data Entry, fine"));
  case  39: return String(T("Volume, fine"));
  case  40: return String(T("Balance, fine"));
  case  41: return String(T("unknown"));
  case  42: return String(T("Pan Position, fine"));
  case  43: return String(T("Expression, fine"));
  case  44: return String(T("Effect Control 1, fine"));
  case  45: return String(T("Effect Control 2, fine"));
  case  46: return String(T("unknown"));
  case  47: return String(T("unknown"));

  case  48: return String(T("unknown"));
  case  49: return String(T("unknown"));
  case  50: return String(T("unknown"));
  case  51: return String(T("unknown"));
  case  52: return String(T("unknown"));
  case  53: return String(T("unknown"));
  case  54: return String(T("unknown"));
  case  55: return String(T("unknown"));
  case  56: return String(T("unknown"));
  case  57: return String(T("unknown"));
  case  58: return String(T("unknown"));
  case  59: return String(T("unknown"));

  case  60: return String(T("unknown"));
  case  61: return String(T("unknown"));
  case  62: return String(T("unknown"));
  case  63: return String(T("unknown"));
  case  64: return String(T("Hold Pedal on/off"));
  case  65: return String(T("Portamento on/off"));
  case  66: return String(T("Sustenuto on/off"));
  case  67: return String(T("Soft Pedal on/off"));
  case  68: return String(T("Legato Pedal on/off"));
  case  69: return String(T("Hold 2 Pedal on/off"));
  case  70: return String(T("Variation"));
  case  71: return String(T("Timbre"));

  case  72: return String(T("Release Time"));
  case  73: return String(T("Attack Time"));
  case  74: return String(T("Brightness"));
  case  75: return String(T("Sound Control 6"));
  case  76: return String(T("Sound Control 7"));
  case  77: return String(T("Sound Control 8"));
  case  78: return String(T("Sound Control 9"));
  case  79: return String(T("Sound Control 10"));
  case  80: return String(T("General Purpose Button 1"));
  case  81: return String(T("General Purpose Button 2"));
  case  82: return String(T("General Purpose Button 3"));
  case  83: return String(T("General Purpose Button 4"));

  case  84: return String(T("unknown"));
  case  85: return String(T("unknown"));
  case  86: return String(T("unknown"));
  case  87: return String(T("unknown"));
  case  88: return String(T("unknown"));
  case  89: return String(T("unknown"));
  case  90: return String(T("unknown"));
  case  91: return String(T("Effects Level"));
  case  92: return String(T("Tremolo Level"));
  case  93: return String(T("Chorus Level"));
  case  94: return String(T("Celeste Level"));
  case  95: return String(T("Phaser Level"));

  case  96: return String(T("Data Button Increment"));
  case  97: return String(T("Data Button Decrement"));
  case  98: return String(T("Non-Registered Parameter, fine"));
  case  99: return String(T("Non-Registered Parameter, coarse"));
  case 100: return String(T("Registered Parameter, fine"));
  case 101: return String(T("Registered Parameter, coarse"));
  case 102: return String(T("unknown"));
  case 103: return String(T("unknown"));
  case 104: return String(T("unknown"));
  case 105: return String(T("unknown"));
  case 106: return String(T("unknown"));
  case 107: return String(T("unknown"));

  case 108: return String(T("unknown"));
  case 109: return String(T("unknown"));
  case 110: return String(T("unknown"));
  case 111: return String(T("unknown"));
  case 112: return String(T("unknown"));
  case 113: return String(T("unknown"));
  case 114: return String(T("unknown"));
  case 115: return String(T("unknown"));
  case 116: return String(T("unknown"));
  case 117: return String(T("unknown"));
  case 118: return String(T("unknown"));
  case 119: return String(T("unknown"));

  case 120: return String(T("All Sound Off"));
  case 121: return String(T("All Controllers Off"));
  case 122: return String(T("Local Keyboard on/off"));
  case 123: return String(T("All Notes Off"));
  case 124: return String(T("Omni Mode Off"));
  case 125: return String(T("Omni Mode On"));
  case 126: return String(T("Mono Operation"));
  case 127: return String(T("Poly Operation"));

  default: return String(T("?"));
  }
}
*/

String midiNoteToString(double midiNoteNumber)
{
  // This is stupid! We should use a StringArray!

  int intNoteNumber = (int) midiNoteNumber;
  switch( intNoteNumber )
  {
  case   0: return String("C-1");
  case   1: return String("C#-1");
  case   2: return String("D-1");
  case   3: return String("D#-1");
  case   4: return String("E-1");
  case   5: return String("F-1");
  case   6: return String("F#-1");
  case   7: return String("G-1");
  case   8: return String("G#-1");
  case   9: return String("A-1");
  case  10: return String("A#-1");
  case  11: return String("B-1");

  case  12: return String("C0");
  case  13: return String("C#0");
  case  14: return String("D0");
  case  15: return String("D#0");
  case  16: return String("E0");
  case  17: return String("F0");
  case  18: return String("F#0");
  case  19: return String("G0");
  case  20: return String("G#0");
  case  21: return String("A0");
  case  22: return String("A#0");
  case  23: return String("B0");

  case  24: return String("C1");
  case  25: return String("C#1");
  case  26: return String("D1");
  case  27: return String("D#1");
  case  28: return String("E1");
  case  29: return String("F1");
  case  30: return String("F#1");
  case  31: return String("G1");
  case  32: return String("G#1");
  case  33: return String("A1");
  case  34: return String("A#1");
  case  35: return String("B1");

  case  36: return String("C2");
  case  37: return String("C#2");
  case  38: return String("D2");
  case  39: return String("D#2");
  case  40: return String("E2");
  case  41: return String("F2");
  case  42: return String("F#2");
  case  43: return String("G2");
  case  44: return String("G#2");
  case  45: return String("A2");
  case  46: return String("A#2");
  case  47: return String("B2");

  case  48: return String("C3");
  case  49: return String("C#3");
  case  50: return String("D3");
  case  51: return String("D#3");
  case  52: return String("E3");
  case  53: return String("F3");
  case  54: return String("F#3");
  case  55: return String("G3");
  case  56: return String("G#3");
  case  57: return String("A3");
  case  58: return String("A#3");
  case  59: return String("B3");

  case  60: return String("C4");
  case  61: return String("C#4");
  case  62: return String("D4");
  case  63: return String("D#4");
  case  64: return String("E4");
  case  65: return String("F4");
  case  66: return String("F#4");
  case  67: return String("G4");
  case  68: return String("G#4");
  case  69: return String("A4");
  case  70: return String("A#4");
  case  71: return String("B4");

  case  72: return String("C5");
  case  73: return String("C#5");
  case  74: return String("D5");
  case  75: return String("D#5");
  case  76: return String("E5");
  case  77: return String("F5");
  case  78: return String("F#5");
  case  79: return String("G5");
  case  80: return String("G#5");
  case  81: return String("A5");
  case  82: return String("A#5");
  case  83: return String("B5");

  case  84: return String("C6");
  case  85: return String("C#6");
  case  86: return String("D6");
  case  87: return String("D#6");
  case  88: return String("E6");
  case  89: return String("F6");
  case  90: return String("F#6");
  case  91: return String("G6");
  case  92: return String("G#6");
  case  93: return String("A6");
  case  94: return String("A#6");
  case  95: return String("B6");

  case  96: return String("C7");
  case  97: return String("C#7");
  case  98: return String("D7");
  case  99: return String("D#7");
  case 100: return String("E7");
  case 101: return String("F7");
  case 102: return String("F#7");
  case 103: return String("G7");
  case 104: return String("G#7");
  case 105: return String("A7");
  case 106: return String("A#7");
  case 107: return String("B7");

  case 108: return String("C8");
  case 109: return String("C#8");
  case 110: return String("D8");
  case 111: return String("D#8");
  case 112: return String("E8");
  case 113: return String("F8");
  case 114: return String("F#8");
  case 115: return String("G8");
  case 116: return String("G#8");
  case 117: return String("A8");
  case 118: return String("A#8");
  case 119: return String("B8");

  case 120: return String("C9");
  case 121: return String("C#9");
  case 122: return String("D9");
  case 123: return String("D#9");
  case 124: return String("E9");
  case 125: return String("F9");
  case 126: return String("F#9");
  case 127: return String("G9");

  case 128: return String("C10");
  case 129: return String("C#10");
  case 130: return String("D10");
  case 131: return String("D#10");
  case 132: return String("E10");
  case 133: return String("F10");
  case 134: return String("F#10");
  case 135: return String("G10");

  case 136: return String("C11");
  case 137: return String("C#11");
  case 138: return String("D11");
  case 139: return String("D#11");
  case 140: return String("E11");
  case 141: return String("F11");
  case 142: return String("F#11");
  case 143: return String("G11");

  default: return String("?");
  }
}

String millisecondsToStringWithUnit0(double value)
{
  return String(round(value)) + String(" ms");
}

String millisecondsToStringWithUnit2(double value)
{
  return String(value, 2) + String(" ms");
}

String octavesToStringWithUnit2(double value)
{
  return String(value, 2) + String(" oct");
}

String percentToStringWithUnit0(double value)
{
  return String( (int) round(value) ) + String(" %");
}

String percentToStringWithUnit1(double value)
{
  return String(value, 1) + String(" %");
}

String percentToStringWithUnit2(double value)
{
  return String(value, 2) + String(" %");
}

String ratioToString0(double value)
{
  double percentage1 = round(100.0*(1.0-value));
  double percentage2 = round(100.0*value);
  return String((int) percentage1 ) + String("/") + String((int) percentage2 );
}

String ratioToString1(double value)
{
  double percentage1 = 100.0*(1.0-value);
  double percentage2 = 100.0*value;
  return String(percentage1, 1) + String("/") + String(percentage2, 1);
}

String ratioBothFullAtCenterToString0(double value)
{
  double p1, p2;
  if( value > 0.5 )
  {
    p2 = 1.0;
    p1 = 2.0 - 2.0*value;
  }
  else
  {
    p1 = 1.0;
    p2 = 2.0 * value;
  }
  p1 = round(100.0*p1);
  p2 = round(100.0*p2);
  return String((int) p1 ) + String("/") + String((int) p2 );
}

String secondsToStringWithUnit2(double value)
{
  return String(value, 2) + String(" s");
}

String secondsToStringWithUnit3(double value)
{
  return String(value, 3) + String(" s");
}

String secondsToStringWithUnit4(double value)
{
  return String(value, 4) + String(" s");
}

String secondsToStringWithUnitTotal4(double value)
{
  if( value >= 100.0 )
    return String(value, 2) + String(" s");
  else if( value >= 1.0 )
    return String(value, 3) + String(" s");
  else if( value >= 0.1 )
    return String(1000*value, 1) + String(" ms");
  else if( value >= 0.01 )
    return String(1000*value, 2) + String(" ms");
  else if( value >= 0.001 )
    return String(1000*value, 3) + String(" ms");
  else if( value >= 0.0001 )
    return String(1000*value, 3) + String(" ms");
  else
    return String(1000*value, 3) + String(" ms");
}

String semitonesToStringWithUnit2(double value)
{
  return String(value, 2) + String(" st");
}

String semitonesToStringWithUnit1(double value)
{
  return String(value, 1) + String(" st");
}

String valueToStringWithTotalNumDigits(double value, int totalNumDigits,
                                              const String& suffix)
{
  if( totalNumDigits <= 0 )
    return String( (int) round(value) ) + suffix;
  else
  {
    if( value >= pow(10.0, totalNumDigits) )
      return String( (int) round(value) ) + suffix;
    else
    {
      int tmp                  = (int) floor(log10(fabs(value)+0.00000001));
      int numDigitsBeforePoint = jmax(1, tmp+1);
      int numDigitsAfterPoint  = totalNumDigits-numDigitsBeforePoint;
      numDigitsAfterPoint      = jmax(0, numDigitsAfterPoint);
      if( numDigitsAfterPoint == 0 )
        return String( (int) round(value) ) + suffix;
      else
        return String(value, numDigitsAfterPoint) + suffix;
    }
  }
}

String valueToString(double value)
{
  return String( value ) ;
}

String valueToString0(double value)
{
  return String( (int) round(value) ) ;
}

String valueToString1(double value)
{
  return String(value, 1);
}

String valueToString2(double value)
{
  return String(value, 2);
}

String valueToString3(double value)
{
  return String(value, 3);
}

String valueToString4(double value)
{
  return String(value, 4);
}

String valueToString5(double value)
{
  return String(value, 5);
}

String valueToStringTotal5(double value)
{
  if( value >= 10000.0 )
    return String( (int) round(value) );
  else if( value >= 1000.0 )
    return String(value, 1);
  else if( value >= 100.0 )
    return String(value, 2);
  else if( value >= 10.0 )
    return String(value, 3);
  else if( value >= 1.0 )
    return String(value, 4);
  else
    return String(value, 5);
}

String valueToStringWithSign0(double value)
{
  return getSignAsString(value) + valueToString0(fabs(value));
}

String valueToStringWithSign1(double value)
{
  return getSignAsString(value) + valueToString1(fabs(value));
}
